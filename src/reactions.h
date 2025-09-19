#include <spark/collisions/mcc.h>

#include <filesystem>

namespace cross_section {
std::shared_ptr<spark::collisions::Reactions<2, 3>> load_electron_reactions_iodine_1(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i_ions);
std::shared_ptr<spark::collisions::Reactions<2, 3>> load_electron_reactions_iodine_2(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i2_ions,
    spark::particle::ChargedSpecies<2, 3>* i_ions,
    spark::particle::ChargedSpecies<2, 3>* in_ions);

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_atomic_ion_reactions_iodine_1(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i_ions);

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_atomic_ion_reactions_iodine_2(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i2_ions);

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_molecular_ion_reactions_iodine_1(
    const std::filesystem::path& dir,
    double neutral_temperature,
    spark::particle::ChargedSpecies<2, 3>* i2_ions);

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_negative_ion_reactions_iodine_1(
    const std::filesystem::path& dir);

std::shared_ptr<spark::collisions::Reactions<2, 3>> load_negative_ion_reactions_iodine_2(
    const std::filesystem::path& dir);

}  // namespace cross_section