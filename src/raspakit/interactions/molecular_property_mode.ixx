module;

export module interactions_molecular_property_mode;

import std;

export namespace Interactions
{
enum class MolecularPropertyMode : std::uint8_t
{
  EnergyOnly,
  EnergyAndPolarizationField,
  EnergyVirialAndPolarizationFieldStrain
};

[[nodiscard]] constexpr bool computesVirial(MolecularPropertyMode mode) noexcept
{
  return mode == MolecularPropertyMode::EnergyVirialAndPolarizationFieldStrain;
}

[[nodiscard]] constexpr bool gathersPolarizationField(MolecularPropertyMode mode) noexcept
{
  return mode != MolecularPropertyMode::EnergyOnly;
}

[[nodiscard]] constexpr bool gathersPolarizationFieldStrain(MolecularPropertyMode mode) noexcept
{
  return mode == MolecularPropertyMode::EnergyVirialAndPolarizationFieldStrain;
}
}  // namespace Interactions
