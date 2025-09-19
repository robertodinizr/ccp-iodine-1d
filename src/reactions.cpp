#include "reactions.h"

#include <spark/collisions/basic_reactions.h>

#include <memory>

#include "rapidcsv.h"
#include "spark/collisions/reaction.h"
#include "spark/constants/constants.h"

using namespace spark::collisions;
using namespace spark::collisions::reactions;

namespace {
CrossSection load_cross_section(const std::filesystem::path& path, double energy_threshold) {
    CrossSection cs;
    const rapidcsv::Document doc(path.string(), rapidcsv::LabelParams(-1, -1),
                                 rapidcsv::SeparatorParams(';'));
    cs.energy = doc.GetColumn<double>(0);
    cs.cross_section = doc.GetColumn<double>(1);
    cs.threshold = energy_threshold;
    return cs;
}
}  // namespace

namespace cross_section {
std::shared_ptr<spark::collisions::Reactions<2, 3>> load_electron_reactions_iodine_1(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i_ions) {
    spark::collisions::Reactions<2, 3> ei_reactions;
    auto cfg = BasicCollisionConfig{.atomic_mass = spark::constants::m_i};

    // Reactions
    // 1. I_elastic
    ei_reactions.push_back(std::make_unique<ElectronElasticCollision<2, 3>>(
        cfg, load_cross_section(dir / "I_elastic.txt", 0.0)));

    // 2. I_excitation{1, 2, 3-25}
    ei_reactions.push_back(std::make_unique<ExcitationCollision<2, 3>>(
        cfg, load_cross_section(dir / "I_excitation_1.txt", 0.9529)));

    ei_reactions.push_back(std::make_unique<ExcitationCollision<2, 3>>(
        cfg, load_cross_section(dir / "I_excitation_2.txt", 7.2015)));

    // 3. I_ionization
    ei_reactions.push_back(std::make_unique<IonizationCollision<2, 3>>(
        i_ions, neutral_temperature, cfg, load_cross_section(dir / "I_ionization.txt", 12.13)));

    return std::make_shared<spark::collisions::Reactions<2, 3>>(std::move(ei_reactions));
}

class DissociativeAttachmentCollision final : public BasicCollision<2, 3> {
public:
    using BasicCollision<2, 3>::BasicCollision;

    DissociativeAttachmentCollision(spark::particle::ChargedSpecies<2, 3>* ions,
                                    const double t_neutral,
                                    const BasicCollisionConfig config,
                                    CrossSection&& cs)
        : BasicCollision<2, 3>(config, std::forward<CrossSection>(cs)),
          ions_(ions),
          t_neutral_(t_neutral) {};

    ReactionOutcome react(spark::particle::ChargedSpecies<2, 3>& projectile,
                          size_t id,
                          double kinetic_energy) override {
        if (kinetic_energy < this->m_cross_section.threshold)
            return ReactionOutcome::NotCollided;

        auto event_pos = projectile.x()[id];

        double ion_mass = ions_->m();
        double neutral_temperature = t_neutral_;

        // Generated ion
        ions_->add(1, [event_pos, ion_mass, neutral_temperature](spark::core::Vec<3>& v,
                                                                 spark::core::Vec<2>& x) {
            x = event_pos;
            const double v_th = std::sqrt(spark::constants::kb * neutral_temperature / ion_mass);
            v = {spark::random::normal(0.0, v_th), spark::random::normal(0.0, v_th),
                 spark::random::normal(0.0, v_th)};
        });

        return ReactionOutcome::Collided | ReactionOutcome::ProjectileToBeRemoved;
    }

private:
    spark::particle::ChargedSpecies<2, 3>* ions_ = nullptr;
    double t_neutral_ = 0.0;
};

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_electron_reactions_iodine_2(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i2_ions,
    spark::particle::ChargedSpecies<2, 3>* i_ions,
    spark::particle::ChargedSpecies<2, 3>* in_ions) {
    Reactions<2, 3> ei2_reactions;
    auto cfg = BasicCollisionConfig{.atomic_mass = 2 * spark::constants::m_i};

    // Reactions

    // 1. e+I2_elastic
    ei2_reactions.push_back(std::make_unique<ElectronElasticCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_elastic.txt", 0.0)));
    ei2_reactions.push_back(std::make_unique<ElectronElasticCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_elastic_momentum.txt", 0.0)));

    // 2. e+I2_excitation{1-4}
    ei2_reactions.push_back(std::make_unique<ExcitationCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_excitation_1.txt", 2.18)));
    ei2_reactions.push_back(std::make_unique<ExcitationCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_excitation_2.txt", 3.0)));
    ei2_reactions.push_back(std::make_unique<ExcitationCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_excitation_3.txt", 5.18)));
    ei2_reactions.push_back(std::make_unique<ExcitationCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_excitation_4.txt", 5.70)));
    ei2_reactions.push_back(std::make_unique<ElectronElasticCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_excitation_vibrational.txt", 0.0)));

    // 3. I2_ionization
    ei2_reactions.push_back(std::make_unique<IonizationCollision<2, 3>>(
        i2_ions, neutral_temperature, cfg, load_cross_section(dir / "I2_ionization.txt", 9.3074)));

    // 4. Dissociative ionization
    ei2_reactions.push_back(std::make_unique<IonizationCollision<2, 3>>(
        i_ions, neutral_temperature, cfg,
        load_cross_section(dir / "I2_dissociative_ionization.txt", 10.8745)));

    // 5. Dissociative attachment
    ei2_reactions.push_back(std::make_unique<DissociativeAttachmentCollision>(
        in_ions, neutral_temperature, cfg,
        load_cross_section(dir / "I2_dissociative_attachment.txt", 0.0)));

    // 6. Dissociation
    ei2_reactions.push_back(std::make_unique<ExcitationCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_dissociation.txt", 1.567)));

    return std::make_shared<spark::collisions::Reactions<2, 3>>(std::move(ei2_reactions));
}

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_atomic_ion_reactions_iodine_1(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i_ions) {
    auto r = std::make_unique<Reactions<2, 3>>();
    auto cfg = BasicCollisionConfig{.atomic_mass = spark::constants::m_i};

    r->push_back(std::make_unique<IonElasticCollision<2, 3>>(
        cfg, load_cross_section(dir / "I_I+_scattering.txt", 0.0)));

    r->push_back(std::make_unique<ChargeExchangeWithRemovalCollision<2, 3>>(
        i_ions, neutral_temperature, cfg,
        load_cross_section(dir / "I_I+_charge_exchange.txt", 0.0)));

    return r;
}

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_atomic_ion_reactions_iodine_2(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i2_ions) {
    auto r = std::make_unique<Reactions<2, 3>>();
    auto cfg = BasicCollisionConfig{};

    r->push_back(std::make_unique<IonDiffMassesElasticCollision<2, 3>>(
        spark::constants::m_i, 2 * spark::constants::m_i,
        load_cross_section(dir / "I2_I+_scattering.txt", 0.0)));

    r->push_back(std::make_unique<ChargeExchangeDiffCollision<2, 3>>(
        i2_ions, neutral_temperature, cfg,
        load_cross_section(dir / "I2_I+_charge_exchange.txt", 0.0)));

    return r;
}

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_molecular_ion_reactions_iodine_1(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i2_ions) {
    auto r = std::make_unique<Reactions<2, 3>>();
    auto cfg = BasicCollisionConfig{.atomic_mass = 2 * spark::constants::m_i};

    r->push_back(std::make_unique<IonElasticCollision<2, 3>>(
        cfg, load_cross_section(dir / "I2_I2+_scattering.txt", 0.0)));

    r->push_back(std::make_unique<ChargeExchangeWithRemovalCollision<2, 3>>(
        i2_ions, neutral_temperature, cfg,
        load_cross_section(dir / "I2_I2+_charge_exchange.txt", 0.0)));

    return r;
}

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_negative_ion_reactions_iodine_1(
    const std::filesystem::path& dir) {
    auto r = std::make_unique<Reactions<2, 3>>();
    auto cfg = BasicCollisionConfig{};

    r->push_back(std::make_unique<SinkCollision<2, 3>>(
        cfg, load_cross_section(dir / "I-_I+_volume_recombination.txt", 0.0)));

    return r;
}

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_negative_ion_reactions_iodine_2(
    const std::filesystem::path& dir) {
    auto r = std::make_unique<Reactions<2, 3>>();
    auto cfg = BasicCollisionConfig{};

    r->push_back(std::make_unique<SinkCollision<2, 3>>(
        cfg, load_cross_section(dir / "I-_I2+_volume_recombination.txt", 0.0)));

    return r;
}

}  // namespace cross_sections