#include "simulation.h"

#include <spark/collisions/mcc.h>
#include <spark/constants/constants.h>
#include <spark/core/matrix.h>
#include <spark/em/electric_field.h>
#include <spark/em/poisson.h>
#include <spark/interpolate/field.h>
#include <spark/interpolate/weight.h>
#include <spark/particle/boundary.h>
#include <spark/particle/pusher.h>
#include <spark/random/random.h>
#include <spark/spatial/grid.h>

#include "reactions.h"

#include <fstream>
#include <filesystem>

namespace {
    auto maxwellian_emitter(double t, double lx, double ly, double m) {
        return [t, lx, ly, m](spark::core::Vec<3>& v, spark::core::Vec<2>& x) {
            x.x = lx * spark::random::uniform();
            x.y = ly * spark::random::uniform();
            double vth = std::sqrt(spark::constants::kb * t / m);
            v = {spark::random::normal(0.0, vth), spark::random::normal(0.0, vth),
                 spark::random::normal(0.0, vth)};
        };
    }
} // namespace

namespace spark {

Simulation::Simulation(const Parameters& parameters, const std::string& data_path)
    : parameters_(parameters), data_path_(data_path), state_(StateInterface(*this)) {}

void Simulation::run() {
    set_initial_conditions();    
    
    auto electron_collisions = load_electron_reactions_iodine_1();
    auto electron_collisions_i2 = load_electron_reactions_iodine_2();
    auto ion_collisions = load_ion_collisions();
    auto ion_collisions_i2 = load_ion_collisions_i2();
    auto ion_collisions_molecular = load_ion_collisions_molecular();
    auto ion_collisions_negative = load_ion_collisions_negative();
    auto ion_collisions_negative_i2 = load_ion_collisions_negative_i2();


    em::StructPoissonSolver2D::DomainProp domain_prop;
    domain_prop.extents = {static_cast<int>(parameters_.nx), static_cast<int>(parameters_.ny)};
    domain_prop.dx = {parameters_.dx, parameters_.dy};

    events().notify(Event::Start, state_);

    std::vector<em::StructPoissonSolver2D::Region> regions;
    
    regions.push_back(em::StructPoissonSolver2D::Region{
        em::CellType::BoundaryNeumann,
        {1, 0},
        {static_cast<int>(parameters_.nx - 2), 0},
        []() { return 0.0; }
    });
    
    regions.push_back(em::StructPoissonSolver2D::Region{
        em::CellType::BoundaryNeumann,
        {1, static_cast<int>(parameters_.ny - 1)},
        {static_cast<int>(parameters_.nx - 2), static_cast<int>(parameters_.ny - 1)},
        []() { return 0.0; }
    });

    regions.push_back(em::StructPoissonSolver2D::Region{
        em::CellType::BoundaryDirichlet,
        {0, 0},
        {0, static_cast<int>(parameters_.ny - 1)},
        []() { return 0.0; }
    });

    double boundary_voltage = 0.0;
    regions.push_back(em::StructPoissonSolver2D::Region{
        em::CellType::BoundaryDirichlet,
        {static_cast<int>(parameters_.nx - 1), 0},
        {static_cast<int>(parameters_.nx - 1), static_cast<int>(parameters_.ny - 1)}, 
        [&boundary_voltage]() { return boundary_voltage; }
    });

    auto poisson_solver = em::StructPoissonSolver2D(domain_prop, regions);

    for (step = 0; step < parameters_.n_steps; ++step) {
        boundary_voltage = parameters_.volt * 
        std::sin(2.0 * spark::constants::pi * parameters_.f * parameters_.dt * static_cast<double>(step));

        spark::interpolate::weight_to_grid(electrons_, electron_density_);
        spark::interpolate::weight_to_grid(ions_, ion_density_);
        spark::interpolate::weight_to_grid(ions_i2_, ion_density_i2_);
        spark::interpolate::weight_to_grid(ions_im_slow_, ion_density_im_);
        
        reduce_rho();

        poisson_solver.solve(phi_field_.data(), rho_field_.data());

        spark::em::electric_field(phi_field_, electric_field_.data());

        spark::interpolate::field_at_particles(electric_field_, electrons_, electron_field);
        spark::interpolate::field_at_particles(electric_field_, ions_, ion_field);
        spark::interpolate::field_at_particles(electric_field_, ions_i2_, ion_field_i2);
        spark::interpolate::field_at_particles(electric_field_, ions_im_slow_, ion_field_im);

        spark::particle::move_particles(electrons_, electron_field, parameters_.dt);
        spark::particle::move_particles(ions_, ion_field, parameters_.dt);
        spark::particle::move_particles(ions_i2_, ion_field_i2, parameters_.dt);
        spark::particle::move_particles(ions_im_slow_, ion_field_im, parameters_.dt);

        tiled_boundary_.apply(electrons_);
        tiled_boundary_.apply(ions_);
        tiled_boundary_.apply(ions_i2_);
        tiled_boundary_.apply(ions_im_slow_);

        electron_collisions.react_all();
        electron_collisions_i2.react_all();
        ion_collisions.react_all();
        ion_collisions_i2.react_all();
        ion_collisions_molecular.react_all();
        ion_collisions_negative.react_all();
        ion_collisions_negative_i2.react_all();

        events().notify(Event::Step, state_);
    }
    events().notify(Event::End, state_);
}

void Simulation::reduce_rho() {
    const auto k = constants::e * parameters_.particle_weight / (parameters_.dx * parameters_.dx);

    auto* rho_ptr = rho_field_.data_ptr();
    auto* ne = electron_density_.data_ptr();
    auto* ni = ion_density_.data_ptr();
    auto* ni_i2 = ion_density_i2_.data_ptr();
    auto* ni_im = ion_density_im_.data_ptr();

    for (size_t i = 0; i < rho_field_.n_total(); ++i) {
        rho_ptr[i] = k * (ni[i] + ni_i2[i] + ni_im[i] - ne[i]);
    }
}

Events<Simulation::Event, Simulation::EventAction>& Simulation::events() {
    return events_;
}

void Simulation::set_initial_conditions() {
    electrons_ = spark::particle::ChargedSpecies<2, 3>(-spark::constants::e, spark::constants::m_e);
    electrons_.add(
        parameters_.n_initial,
        maxwellian_emitter(parameters_.te, parameters_.lx, parameters_.ly, spark::constants::m_e));

    ions_ = spark::particle::ChargedSpecies<2, 3>(spark::constants::e, parameters_.m_i);
    ions_.add(
        parameters_.n_initial,
        maxwellian_emitter(parameters_.ti, parameters_.lx, parameters_.ly, parameters_.m_i));

    ions_i2_ = spark::particle::ChargedSpecies<2, 3>(spark::constants::e, 2.0 * parameters_.m_i);
    ions_i2_.add(
        parameters_.n_initial,
        maxwellian_emitter(parameters_.ti, parameters_.lx, parameters_.ly, 2.0 * parameters_.m_i));
        
    ions_im_slow_ = spark::particle::ChargedSpecies<2, 3>(-spark::constants::e, parameters_.m_i);
    ions_im_slow_.add(
        parameters_.n_initial,
        maxwellian_emitter(parameters_.ti, parameters_.lx, parameters_.ly, parameters_.m_i));

    electron_density_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                                      {parameters_.nx, parameters_.ny});
    ion_density_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                                 {parameters_.nx, parameters_.ny});
    ion_density_i2_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                                 {parameters_.nx, parameters_.ny});
    ion_density_im_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                                 {parameters_.nx, parameters_.ny});                                                                                          
    rho_field_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                               {parameters_.nx, parameters_.ny});
    phi_field_ = spark::spatial::UniformGrid<2>({parameters_.lx, parameters_.ly},
                                               {parameters_.nx, parameters_.ny});

    electric_field_ = spark::spatial::TUniformGrid<core::TVec<double, 2>, 2>(
        {parameters_.lx, parameters_.ly}, {parameters_.nx, parameters_.ny});

    std::vector<spark::particle::TiledBoundary> boundaries = {
        {{-1, -1}, {static_cast<int>(parameters_.nx - 1), -1}, spark::particle::BoundaryType::Specular},
        {{0, static_cast<int>(parameters_.ny - 1)}, {static_cast<int>(parameters_.nx - 2), static_cast<int>(parameters_.ny - 1)}, spark::particle::BoundaryType::Specular},
        {{-1, 0}, {-1, static_cast<int>(parameters_.ny)}, spark::particle::BoundaryType::Absorbing},
        {{static_cast<int>(parameters_.nx - 1), -1}, {static_cast<int>(parameters_.nx - 1), static_cast<int>(parameters_.ny - 1)}, spark::particle::BoundaryType::Absorbing}
    };
    tiled_boundary_ = spark::particle::TiledBoundary2D(electric_field_.prop(), boundaries, parameters_.dt);
}

spark::collisions::MCCReactionSet<2, 3> Simulation::load_electron_reactions_iodine_1() {
    auto electron_reactions_iodine = cross_section::load_electron_reactions_iodine_1(data_path_, parameters_.tg, &ions_);
    spark::collisions::ReactionConfig<2, 3> electron_reaction_config{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ni, parameters_.tg),
        .reactions = electron_reactions_iodine,
        .dyn = spark::collisions::RelativeDynamics::FastProjectile
    };
    return spark::collisions::MCCReactionSet(&electrons_, std::move(electron_reaction_config));
}

spark::collisions::MCCReactionSet<2, 3> Simulation::load_electron_reactions_iodine_2() {
    auto electron_reactions_iodine2 = cross_section::load_electron_reactions_iodine_2(data_path_, parameters_.tg, &ions_i2_, &ions_, &ions_im_slow_);
    spark::collisions::ReactionConfig<2, 3> electron_reaction_config_i2{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ni2, parameters_.tg),
        .reactions = electron_reactions_iodine2,
        .dyn = spark::collisions::RelativeDynamics::FastProjectile
    };
    return spark::collisions::MCCReactionSet(&electrons_, std::move(electron_reaction_config_i2));
}

spark::collisions::MCCReactionSet<2, 3> Simulation::load_ion_collisions() {
    auto ion_reactions_iodine = cross_section::load_atomic_ion_reactions_iodine_1(data_path_, parameters_.tg, &ions_);
    spark::collisions::ReactionConfig<2, 3> ion_reaction_config{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ni, parameters_.tg),
        .reactions = ion_reactions_iodine,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_, std::move(ion_reaction_config));
}

spark::collisions::MCCReactionSet<2, 3> Simulation::load_ion_collisions_i2() {
    auto ion_reactions_iodine_i2 = cross_section::load_atomic_ion_reactions_iodine_2(data_path_, parameters_.tg, &ions_i2_);
    spark::collisions::ReactionConfig<2, 3> ion_reaction_config_i2{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ni2, parameters_.tg),
        .reactions = ion_reactions_iodine_i2,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_i2_, std::move(ion_reaction_config_i2));
}

spark::collisions::MCCReactionSet<2, 3> Simulation::load_ion_collisions_molecular() {
    auto ion_reactions_molecular = cross_section::load_molecular_ion_reactions_iodine_1(data_path_, parameters_.tg, &ions_i2_);
    spark::collisions::ReactionConfig<2, 3> ion_reaction_config_molecular{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ni, parameters_.tg),
        .reactions = ion_reactions_molecular,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_i2_, std::move(ion_reaction_config_molecular));
}

spark::collisions::MCCReactionSet<2, 3> Simulation::load_ion_collisions_negative() {
    auto ion_reactions_negative = cross_section::load_negative_ion_reactions_iodine_1(data_path_);
    spark::collisions::ReactionConfig<2, 3> ion_reaction_config_im{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ni, parameters_.tg),
        .reactions = ion_reactions_negative,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_im_slow_, std::move(ion_reaction_config_im));
}

spark::collisions::MCCReactionSet<2, 3> Simulation::load_ion_collisions_negative_i2() {
    auto ion_reactions_negative_i2 = cross_section::load_negative_ion_reactions_iodine_2(data_path_);
    spark::collisions::ReactionConfig<2, 3> ion_reaction_config{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<2, 3>>(parameters_.ni2, parameters_.tg),
        .reactions = ion_reactions_negative_i2,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_im_slow_, std::move(ion_reaction_config));
}


} // namespace spark