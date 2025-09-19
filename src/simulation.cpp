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
    auto maxwellian_emitter(double t, double lx, double m) {
        return [t, lx, m](spark::core::Vec<3>& v, spark::core::Vec<1>& x) {
            x.x = lx * spark::random::uniform();
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

    events().notify(Event::Start, state_);

    double boundary_voltage = 0.0;

    auto poisson_solver = em::ThomasPoissonSolver1D(parameters_.nx, parameters_.dx);

    for (step = 0; step < parameters_.n_steps; ++step) {
        boundary_voltage = parameters_.volt *
        std::sin(2.0 * constants::pi * parameters_.f * parameters_.dt * static_cast<double>(step));

        interpolate::weight_to_grid(electrons_, electron_density_);
        interpolate::weight_to_grid(ions_, ion_density_);
        interpolate::weight_to_grid(ions_i2_, ion_density_i2_);
        interpolate::weight_to_grid(ions_im_slow_, ion_density_im_);
        
        reduce_rho();

    poisson_solver.solve(rho_field_.data().data(), phi_field_.data().data(), 0.0, boundary_voltage);

        em::electric_field(phi_field_, electric_field_.data());

        interpolate::field_at_particles(electric_field_, electrons_, electron_field);
        interpolate::field_at_particles(electric_field_, ions_, ion_field);
        interpolate::field_at_particles(electric_field_, ions_i2_, ion_field_i2);
        interpolate::field_at_particles(electric_field_, ions_im_slow_, ion_field_im);

        particle::move_particles(electrons_, electron_field, parameters_.dt);
        particle::move_particles(ions_, ion_field, parameters_.dt);
        particle::move_particles(ions_i2_, ion_field_i2, parameters_.dt);
        particle::move_particles(ions_im_slow_, ion_field_im, parameters_.dt);

        particle::apply_absorbing_boundary(electrons_, 0, parameters_.lx);
        particle::apply_absorbing_boundary(ions_, 0, parameters_.lx);
        particle::apply_absorbing_boundary(ions_i2_, 0, parameters_.lx);
        particle::apply_absorbing_boundary(ions_im_slow_, 0, parameters_.lx);

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
    const auto k = constants::e * parameters_.particle_weight / (parameters_.dx);

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
    electrons_ = spark::particle::ChargedSpecies<1, 3>(-spark::constants::e, spark::constants::m_e);
    electrons_.add(
        parameters_.n_initial,
        maxwellian_emitter(parameters_.te, parameters_.lx, spark::constants::m_e));

    ions_ = spark::particle::ChargedSpecies<1, 3>(spark::constants::e, parameters_.m_i);
    ions_.add(
        parameters_.n_initial,
        maxwellian_emitter(parameters_.ti, parameters_.lx, parameters_.m_i));

    ions_i2_ = spark::particle::ChargedSpecies<1, 3>(spark::constants::e, 2.0 * parameters_.m_i);
    ions_i2_.add(
        parameters_.n_initial,
        maxwellian_emitter(parameters_.ti, parameters_.lx, 2.0 * parameters_.m_i));
        
    ions_im_slow_ = spark::particle::ChargedSpecies<1, 3>(-spark::constants::e, parameters_.m_i);
    ions_im_slow_.add(
        parameters_.n_initial,
        maxwellian_emitter(parameters_.ti, parameters_.lx, parameters_.m_i));

    // Construct grids using core::Vec and core::ULongVec as required by TUniformGrid constructors
    spark::core::Vec<1> l{parameters_.lx};
    spark::core::ULongVec<1> n{parameters_.nx};

    electron_density_ = spark::spatial::UniformGrid<1>(l, n);
    ion_density_ = spark::spatial::UniformGrid<1>(l, n);
    ion_density_i2_ = spark::spatial::UniformGrid<1>(l, n);
    ion_density_im_ = spark::spatial::UniformGrid<1>(l, n);
    rho_field_ = spark::spatial::UniformGrid<1>(l, n);
    phi_field_ = spark::spatial::UniformGrid<1>(l, n);

    electric_field_ = spark::spatial::TUniformGrid<core::TVec<double, 1>, 1>(l, n);
}

spark::collisions::MCCReactionSet<1, 3> Simulation::load_electron_reactions_iodine_1() {
    auto electron_reactions_iodine = cross_section::load_electron_reactions_iodine_1(data_path_, parameters_.tg, &ions_);
    spark::collisions::ReactionConfig<1, 3> electron_reaction_config{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<1, 3>>(parameters_.ni, parameters_.tg),
        .reactions = electron_reactions_iodine,
        .dyn = spark::collisions::RelativeDynamics::FastProjectile
    };
    return spark::collisions::MCCReactionSet(&electrons_, std::move(electron_reaction_config));
}

spark::collisions::MCCReactionSet<1, 3> Simulation::load_electron_reactions_iodine_2() {
    auto electron_reactions_iodine2 = cross_section::load_electron_reactions_iodine_2(data_path_, parameters_.tg, &ions_i2_, &ions_, &ions_im_slow_);
    spark::collisions::ReactionConfig<1, 3> electron_reaction_config_i2{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<1, 3>>(parameters_.ni2, parameters_.tg),
        .reactions = electron_reactions_iodine2,
        .dyn = spark::collisions::RelativeDynamics::FastProjectile
    };
    return spark::collisions::MCCReactionSet(&electrons_, std::move(electron_reaction_config_i2));
}

spark::collisions::MCCReactionSet<1, 3> Simulation::load_ion_collisions() {
    auto ion_reactions_iodine = cross_section::load_atomic_ion_reactions_iodine_1(data_path_, parameters_.tg, &ions_);
    spark::collisions::ReactionConfig<1, 3> ion_reaction_config{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<1, 3>>(parameters_.ni, parameters_.tg),
        .reactions = ion_reactions_iodine,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_, std::move(ion_reaction_config));
}

spark::collisions::MCCReactionSet<1, 3> Simulation::load_ion_collisions_i2() {
    auto ion_reactions_iodine_i2 = cross_section::load_atomic_ion_reactions_iodine_2(data_path_, parameters_.tg, &ions_i2_);
    spark::collisions::ReactionConfig<1, 3> ion_reaction_config_i2{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<1, 3>>(parameters_.ni2, parameters_.tg),
        .reactions = ion_reactions_iodine_i2,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_i2_, std::move(ion_reaction_config_i2));
}

spark::collisions::MCCReactionSet<1, 3> Simulation::load_ion_collisions_molecular() {
    auto ion_reactions_molecular = cross_section::load_molecular_ion_reactions_iodine_1(data_path_, parameters_.tg, &ions_i2_);
    spark::collisions::ReactionConfig<1, 3> ion_reaction_config_molecular{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<1, 3>>(parameters_.ni, parameters_.tg),
        .reactions = ion_reactions_molecular,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_i2_, std::move(ion_reaction_config_molecular));
}

spark::collisions::MCCReactionSet<1, 3> Simulation::load_ion_collisions_negative() {
    auto ion_reactions_negative = cross_section::load_negative_ion_reactions_iodine_1(data_path_);
    spark::collisions::ReactionConfig<1, 3> ion_reaction_config_im{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<1, 3>>(parameters_.ni, parameters_.tg),
        .reactions = ion_reactions_negative,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_im_slow_, std::move(ion_reaction_config_im));
}

spark::collisions::MCCReactionSet<1, 3> Simulation::load_ion_collisions_negative_i2() {
    auto ion_reactions_negative_i2 = cross_section::load_negative_ion_reactions_iodine_2(data_path_);
    spark::collisions::ReactionConfig<1, 3> ion_reaction_config{
        .dt = parameters_.dt,
        .target = std::make_shared<spark::collisions::StaticUniformTarget<1, 3>>(parameters_.ni2, parameters_.tg),
        .reactions = ion_reactions_negative_i2,
        .dyn = spark::collisions::RelativeDynamics::SlowProjectile
    };
    return spark::collisions::MCCReactionSet(&ions_im_slow_, std::move(ion_reaction_config));
}


} // namespace spark