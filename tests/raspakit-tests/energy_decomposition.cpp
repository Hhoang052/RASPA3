#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "irmof_fixtures.hpp"

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import energy_status;
import energy_status_inter;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import integrators;
import integrators_update;
import integrators_compute;
import interpolation_energy_grid;
import cif_reader;
import mc_moves_probabilities;

namespace
{

void printEnergyDecomposition(const RunningEnergy &energy, const std::string &label)
{
  const double conv = Units::EnergyToKelvin;
  const auto toKelvin = [&](double internal) { return conv * internal; };
  const auto toKJmol = [&](double internal) { return conv * internal / 1.2027242847; };
  std::print("=== {} ===\n", label);
  std::print("  Total:                 {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.potentialEnergy()),
             toKelvin(energy.potentialEnergy()));
  std::print("  Host/ads VDW:          {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.frameworkMoleculeVDW),
             toKelvin(energy.frameworkMoleculeVDW));
  std::print("  Host/ads Real:         {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.frameworkMoleculeCharge),
             toKelvin(energy.frameworkMoleculeCharge));
  std::print("  Ewald Fourier:         {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.ewald_fourier),
             toKelvin(energy.ewald_fourier));
  std::print("  Ewald self:            {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.ewald_self),
             toKelvin(energy.ewald_self));
  std::print("  Ewald exclusion:       {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.ewald_exclusion),
             toKelvin(energy.ewald_exclusion));
  std::print("  Total VDW:             {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.VanDerWaalsEnergy()),
             toKelvin(energy.VanDerWaalsEnergy()));
  std::print("  Total Coulomb:         {: .2f} kJ/mol  ({: .2f} K)\n", toKJmol(energy.CoulombEnergy()),
             toKelvin(energy.CoulombEnergy()));
  std::print("\n");
}

void expectEnergyInterNear(const EnergyInter& lhs, const EnergyInter& rhs, double tolerance = 1e-10)
{
  EXPECT_NEAR(lhs.VanDerWaals.energy, rhs.VanDerWaals.energy, tolerance);
  EXPECT_NEAR(lhs.VanDerWaals.dUdlambda, rhs.VanDerWaals.dUdlambda, tolerance);
  EXPECT_NEAR(lhs.VanDerWaalsTailCorrection.energy, rhs.VanDerWaalsTailCorrection.energy, tolerance);
  EXPECT_NEAR(lhs.VanDerWaalsTailCorrection.dUdlambda, rhs.VanDerWaalsTailCorrection.dUdlambda, tolerance);
  EXPECT_NEAR(lhs.CoulombicReal.energy, rhs.CoulombicReal.energy, tolerance);
  EXPECT_NEAR(lhs.CoulombicReal.dUdlambda, rhs.CoulombicReal.dUdlambda, tolerance);
  EXPECT_NEAR(lhs.CoulombicFourier.energy, rhs.CoulombicFourier.energy, tolerance);
  EXPECT_NEAR(lhs.CoulombicFourier.dUdlambda, rhs.CoulombicFourier.dUdlambda, tolerance);
  EXPECT_NEAR(lhs.totalInter.energy, rhs.totalInter.energy, tolerance);
  EXPECT_NEAR(lhs.totalInter.dUdlambda, rhs.totalInter.dUdlambda, tolerance);
}

void expectDetailedEnergyNear(const EnergyStatus& lhs, const EnergyStatus& rhs, double tolerance = 1e-10)
{
  EXPECT_NEAR(lhs.totalEnergy.energy, rhs.totalEnergy.energy, tolerance);
  EXPECT_NEAR(lhs.totalEnergy.dUdlambda, rhs.totalEnergy.dUdlambda, tolerance);
  EXPECT_NEAR(lhs.polarizationEnergy.energy, rhs.polarizationEnergy.energy, tolerance);
  EXPECT_NEAR(lhs.polarizationEnergy.dUdlambda, rhs.polarizationEnergy.dUdlambda, tolerance);
  EXPECT_NEAR(lhs.dUdlambda, rhs.dUdlambda, tolerance);
  EXPECT_NEAR(lhs.translationalKineticEnergy, rhs.translationalKineticEnergy, tolerance);
  EXPECT_NEAR(lhs.rotationalKineticEnergy, rhs.rotationalKineticEnergy, tolerance);
  EXPECT_NEAR(lhs.noseHooverEnergy, rhs.noseHooverEnergy, tolerance);
  ASSERT_EQ(lhs.intraComponentEnergies.size(), rhs.intraComponentEnergies.size());
  ASSERT_EQ(lhs.externalFieldComponentEnergies.size(), rhs.externalFieldComponentEnergies.size());
  ASSERT_EQ(lhs.frameworkComponentEnergies.size(), rhs.frameworkComponentEnergies.size());
  ASSERT_EQ(lhs.interComponentEnergies.size(), rhs.interComponentEnergies.size());
  for (std::size_t i = 0; i < lhs.intraComponentEnergies.size(); ++i)
    EXPECT_NEAR(lhs.intraComponentEnergies[i].total().energy, rhs.intraComponentEnergies[i].total().energy, tolerance);
  for (std::size_t i = 0; i < lhs.externalFieldComponentEnergies.size(); ++i)
    expectEnergyInterNear(lhs.externalFieldComponentEnergies[i], rhs.externalFieldComponentEnergies[i], tolerance);
  for (std::size_t i = 0; i < lhs.frameworkComponentEnergies.size(); ++i)
    expectEnergyInterNear(lhs.frameworkComponentEnergies[i], rhs.frameworkComponentEnergies[i], tolerance);
  for (std::size_t i = 0; i < lhs.interComponentEnergies.size(); ++i)
    expectEnergyInterNear(lhs.interComponentEnergies[i], rhs.interComponentEnergies[i], tolerance);
}
}  // namespace

TEST(energy_decomposition, CO2_in_IRMOF1_OMS_geometry)
{
  TemporaryDirectory fixtureDir;
  fixtureDir.write("force_field.json", irmof_fixtures::kMinimizationForceFieldJson);
  fixtureDir.write("CO2.json", irmof_fixtures::kCO2Json);

  ForceField forceField = ForceField::readForceField(fixtureDir.path().string(), "force_field.json").value();
  forceField.chargeMethod = ForceField::ChargeMethod::Ewald;
  forceField.useCharge = true;

  const auto cif = CIFReader::readCIFString(std::string(irmof_fixtures::kIrmof1Cif), forceField,
                                            CIFReader::UseChargesFrom::PseudoAtoms);
  ASSERT_TRUE(cif.has_value());
  auto [simulationBox, spaceGroupHallNumber, definedAtoms, fractionalAtomsUnitCell] = cif.value();
  Framework framework = Framework(forceField, "IRMOF-1", simulationBox, spaceGroupHallNumber, definedAtoms,
                                  fractionalAtomsUnitCell, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  Component co2(Component::Type::Adsorbate, 0, forceField, "CO2", (fixtureDir.path() / "CO2").string(), 5, 21,
                probabilities, std::nullopt, false);

  auto evaluateAtPositions = [&](const std::string &label, const std::array<double3, 3> &positions)
  {
    System system = System(forceField, std::nullopt, false, 300.0, 0.0, 0.81, framework, {co2}, {}, {1}, 5);
    std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
    for (std::size_t atom = 0; atom < positions.size(); ++atom)
    {
      moleculeAtoms[atom].position = positions[atom];
    }
    RunningEnergy energy = system.computeTotalEnergies();
    printEnergyDecomposition(energy, label);
    return energy;
  };

  // These positions were recorded against the old P1 IRMOF-1.cif. The current symmetry-based
  // CIF describes the same crystal in a setting shifted by c/2, so the z coordinates carry a
  // +12.916 Angstrom (half unit cell) translation relative to the original values.
  const std::array<double3, 3> raspa2Oms = {double3(10.136569193677, 10.136569193677, 20.335907533962),
                                            double3(9.417951022032, 9.417951022032, 20.871982945414),
                                            double3(8.699332850387, 8.699332850387, 21.408058356865)};
  const std::array<double3, 3> raspa3Pore = {double3(9.4276655408314216, 17.808185239505558, 3.488334459168581),
                                             double3(9.919740710490089, 18.722474610990751, 2.996259289509912),
                                             double3(10.411815880148756, 19.636763982475944, 2.504184119851244)};

  const RunningEnergy omsEnergy = evaluateAtPositions("RASPA3 at RASPA2 OMS geometry", raspa2Oms);
  evaluateAtPositions("RASPA3 at RASPA3 pore geometry (seed 12345)", raspa3Pore);
  EXPECT_NEAR(omsEnergy.potentialEnergy() * Units::EnergyToKelvin, -2701.54, 5.0);
}

TEST(energy_decomposition, CO2_Methane_in_Box)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);

  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane, co2}, {}, {15, 30}, 5);

  RunningEnergy energy = system.computeTotalEnergies();

  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, system.atomData, system.atomDynamics);
  strainDerivative.first.sumTotal();

  EXPECT_NEAR(energy.moleculeMoleculeVDW + energy.moleculeMoleculeCharge, strainDerivative.first.totalEnergy.energy,
              1e-6);
}

TEST(energy_decomposition, CO2_Methane_in_Box_Ewald)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);
  Component co2_2 = Component::makeCO2(forceField, 2, true);

  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2, co2_2}, {}, {15, 30}, 5);

  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);

  RunningEnergy energy = system.computeTotalEnergies();

  Interactions::computeEwaldFourierEnergySingleIon(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                   system.forceField, system.simulationBox, double3(0.0, 0.0, 0.0),
                                                   1.0);

  system.precomputeTotalRigidEnergy();
  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.framework, system.components,
      system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
      system.netChargeFramework, system.netChargePerComponent);

  strainDerivative.first.sumTotal();

  EXPECT_NEAR(energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion,
              strainDerivative.first.totalEnergy.energy, 1e-6);
}

inline std::pair<EnergyStatus, double3x3> pair_acc(const std::pair<EnergyStatus, double3x3> &lhs,
                                                   const std::pair<EnergyStatus, double3x3> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

TEST(energy_decomposition, CO2_Methane_in_Framework)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));
  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {methane, co2}, {}, {10, 15}, 5);

  system.precomputeTotalRigidEnergy();

  RunningEnergy energy = system.computeTotalEnergies();

  RunningEnergy energyForces = Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox,
      system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik,
      system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent);

  std::pair<EnergyStatus, double3x3> strainDerivative = system.computeMolecularPressure();
  EnergyStatus energyOnly = system.computeMolecularEnergyStatus();

  EXPECT_NEAR(energy.potentialEnergy(), strainDerivative.first.totalEnergy.energy, 1e-6);
  EXPECT_NEAR(energy.potentialEnergy(), energyForces.potentialEnergy(), 1e-6);
  expectDetailedEnergyNear(strainDerivative.first, energyOnly, 1e-8);
}

TEST(energy_decomposition, ComputePressureFalseKeepsDetailedEnergyAndReturnsZeroTensor)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  Component co2 = Component::makeCO2(forceField, 1, true);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane, co2}, {}, {8, 12}, 5);

  const std::pair<EnergyStatus, double3x3> pressure = system.computeMolecularPressure();
  const auto eikXBefore = system.eik_x;
  const auto eikYBefore = system.eik_y;
  const auto eikZBefore = system.eik_z;
  const auto eikXYBefore = system.eik_xy;
  const auto fixedEikBefore = system.fixedFrameworkStoredEik;
  const auto storedEikBefore = system.storedEik;
  std::vector<double3> gradientsBefore;
  gradientsBefore.reserve(system.atomDynamics.size());
  for (const AtomDynamics& dynamics : system.atomDynamics) gradientsBefore.push_back(dynamics.gradient);

  system.computePressure = false;
  const std::pair<EnergyStatus, double3x3> energyOnly = system.computeMolecularPropertiesForSampling();

  expectDetailedEnergyNear(pressure.first, energyOnly.first, 1e-8);
  for (double value : energyOnly.second.m) EXPECT_DOUBLE_EQ(value, 0.0);
  EXPECT_EQ(system.eik_x, eikXBefore);
  EXPECT_EQ(system.eik_y, eikYBefore);
  EXPECT_EQ(system.eik_z, eikZBefore);
  EXPECT_EQ(system.eik_xy, eikXYBefore);
  EXPECT_EQ(system.fixedFrameworkStoredEik, fixedEikBefore);
  EXPECT_EQ(system.storedEik, storedEikBefore);
  ASSERT_EQ(system.atomDynamics.size(), gradientsBefore.size());
  for (std::size_t i = 0; i < gradientsBefore.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(system.atomDynamics[i].gradient.x, gradientsBefore[i].x);
    EXPECT_DOUBLE_EQ(system.atomDynamics[i].gradient.y, gradientsBefore[i].y);
    EXPECT_DOUBLE_EQ(system.atomDynamics[i].gradient.z, gradientsBefore[i].z);
  }
}

TEST(energy_decomposition, EnergyOnlyMatchesPressurePathForFiniteCutoffChargeMethods)
{
  for (const ForceField::ChargeMethod method :
       {ForceField::ChargeMethod::Wolf, ForceField::ChargeMethod::DampedShiftedForce,
        ForceField::ChargeMethod::ModifiedShiftedForce, ForceField::ChargeMethod::ZeroDipole})
  {
    ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
    forceField.chargeMethod = method;
    forceField.omitEwaldFourier = true;
    forceField.EwaldAlpha = 0.25;
    Component co2 = Component::makeCO2(forceField, 0, true);
    System system =
        System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {4}, 5);

    const EnergyStatus pressure = system.computeMolecularPressure().first;
    const EnergyStatus energyOnly = system.computeMolecularEnergyStatus();
    expectDetailedEnergyNear(pressure, energyOnly, 1e-8);
  }
}

TEST(energy_decomposition, EnergyOnlyHandlesZeroMolecules)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, true);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {0}, 5);

  const EnergyStatus pressure = system.computeMolecularPressure().first;
  system.computePressure = false;
  const auto energyOnly = system.computeMolecularPropertiesForSampling();
  expectDetailedEnergyNear(pressure, energyOnly.first, 0.0);
  for (double value : energyOnly.second.m) EXPECT_DOUBLE_EQ(value, 0.0);
}

TEST(energy_decomposition, ExternalFieldEnergyIsIncludedInDetailedEnergy)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction;
  Component methane = Component::makeMethane(forceField, 0);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), true, 300.0, 1e4, 1.0, {}, {methane}, {}, {1}, 5);
  system.spanOfMoleculeAtoms().front().position = double3(2.0, 3.0, 4.0);

  const RunningEnergy running = system.computeTotalEnergies();
  EnergyStatus detailed = system.computeMolecularEnergyStatus();
  const EnergyStatus pressureDetailed = system.computeMolecularPressure().first;

  EXPECT_NEAR(running.externalFieldVDW, 24.0, 1e-12);
  EXPECT_NEAR(detailed.externalFieldComponentEnergy(0, 0).VanDerWaals.energy, running.externalFieldVDW, 1e-12);
  EXPECT_NEAR(detailed.externalFieldMoleculeEnergy.VanDerWaals.energy, running.externalFieldVDW, 1e-12);
  EXPECT_DOUBLE_EQ(detailed.externalFieldComponentEnergy(0, 0).VanDerWaals.dUdlambda, 0.0);
  EXPECT_NEAR(detailed.totalEnergy.energy, running.potentialEnergy(), 1e-8);
  expectDetailedEnergyNear(pressureDetailed, detailed, 1e-12);
}
