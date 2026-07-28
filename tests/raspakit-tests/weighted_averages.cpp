#include <gtest/gtest.h>

import std;

import int3;
import double3x3;
import forcefield;
import framework;
import component;
import system;
import mc_moves_move_types;
import property_block_average;
import property_lambda_probability_histogram;
import property_pressure;

TEST(weighted_averages, weights_smaller_than_one_are_normalized_by_the_actual_weight)
{
  BlockAverage<double> average(2);
  average.addSample(0, 10.0, 0.01);
  average.addSample(0, 20.0, 0.01);
  average.addSample(1, 30.0, 0.02);

  EXPECT_DOUBLE_EQ(average.averaged(0), 15.0);
  EXPECT_DOUBLE_EQ(average.averaged(1), 30.0);
  EXPECT_DOUBLE_EQ(average.averaged(), 22.5);
}

TEST(weighted_averages, block_average_is_invariant_under_global_weight_rescaling)
{
  BlockAverage<double> reference(3);
  BlockAverage<double> rescaled(3);

  for (std::size_t block = 0; block < 3; ++block)
  {
    const double sample = 2.0 + static_cast<double>(block);
    const double weight = 0.25 * static_cast<double>(block + 1);
    reference.addSample(block, sample, weight);
    rescaled.addSample(block, sample, 1.0e-12 * weight);
  }

  const auto [referenceMean, referenceError] = reference.average();
  const auto [rescaledMean, rescaledError] = rescaled.average();
  EXPECT_NEAR(rescaledMean, referenceMean, 1.0e-14);
  EXPECT_NEAR(rescaledError, referenceError, 1.0e-14);
}

TEST(weighted_averages, histogram_is_invariant_under_global_weight_rescaling)
{
  BlockHistogram<double> reference(1, 1, 2);
  BlockHistogram<double> rescaled(1, 1, 2);

  reference(0, 0, 0) += 0.2 * 3.0;
  reference(0, 0, 1) += 0.2 * 7.0;
  reference.addCount(0, 0.2);

  rescaled(0, 0, 0) += 1.0e-12 * 0.2 * 3.0;
  rescaled(0, 0, 1) += 1.0e-12 * 0.2 * 7.0;
  rescaled.addCount(0, 1.0e-12 * 0.2);

  const std::vector<double> referenceAverage = reference.averaged(0);
  const std::vector<double> rescaledAverage = rescaled.averaged(0);
  ASSERT_EQ(referenceAverage.size(), rescaledAverage.size());
  EXPECT_NEAR(rescaledAverage[0], referenceAverage[0], 1.0e-14);
  EXPECT_NEAR(rescaledAverage[1], referenceAverage[1], 1.0e-14);
}

TEST(weighted_averages, non_empty_zero_weight_block_is_reported_as_invalid)
{
  BlockAverage<double> average(1);
  average.addSample(0, 12.0, 0.0);
  EXPECT_THROW(average.averaged(), std::runtime_error);

  BlockHistogram<double> histogram(1, 1, 1);
  histogram.addCount(0, 0.0);
  EXPECT_THROW(histogram.averaged(0), std::runtime_error);
}

TEST(lambda_reweighting, inactive_component_bias_does_not_enter_system_weight)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework framework = Framework::makeMFI(forceField, int3(1, 1, 1));
  Component inactive = Component::makeMethane(forceField, 0);
  Component active = Component::makeCO2(forceField, 1, true);

  System system =
      System(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, {framework}, {inactive, active}, {}, {0, 0}, 5);
  system.components[1].mc_moves_probabilities.setProbability(Move::Types::SwapCFCMC, 0.2);
  system.determineFractionalComponents();

  ASSERT_FALSE(system.gcLambdaActive(0));
  ASSERT_TRUE(system.gcLambdaActive(1));
  ASSERT_FALSE(system.components[0].lambdaGC.biasFactor.empty());
  ASSERT_FALSE(system.components[1].lambdaGC.biasFactor.empty());

  system.components[0].lambdaGC.biasFactor[system.components[0].lambdaGC.currentBin] = 1000.0;
  system.components[1].lambdaGC.biasFactor[system.components[1].lambdaGC.currentBin] = 2.0;

  EXPECT_NEAR(system.lambdaLogWeight(), -2.0, 1.0e-14);
  EXPECT_NEAR(system.weight(), std::exp(-2.0), 1.0e-14);
}

TEST(lambda_reweighting, multi_component_system_without_cfc_has_unit_weight)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework framework = Framework::makeMFI(forceField, int3(1, 1, 1));
  Component methane = Component::makeMethane(forceField, 0);
  Component carbonDioxide = Component::makeCO2(forceField, 1, true);
  System system = System(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, {framework},
                         {methane, carbonDioxide}, {}, {0, 0}, 5);

  system.components[0].lambdaGC.biasFactor[system.components[0].lambdaGC.currentBin] = 1000.0;
  system.components[1].lambdaGC.biasFactor[system.components[1].lambdaGC.currentBin] = 2000.0;

  EXPECT_FALSE(system.hasAnyActiveLambdaCoordinate());
  EXPECT_DOUBLE_EQ(system.weight(), 1.0);
}

TEST(lambda_reweighting, fixed_lambda_is_reportable_but_not_adaptively_reweighted)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework framework = Framework::makeMFI(forceField, int3(1, 1, 1));
  Component component = Component::makeMethane(forceField, 0);
  component.fixedLambdaBin = 10;
  component.fixedLambdaCoordinate = Component::FixedLambdaCoordinate::GC;

  System system =
      System(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, {framework}, {component}, {}, {0}, 5);

  ASSERT_TRUE(system.gcLambdaActive(0));
  ASSERT_FALSE(system.gcLambdaMovesEnabled(0));
  ASSERT_FALSE(system.gcLambdaAdaptiveBiasEnabled(0));
  system.components[0].lambdaGC.biasFactor[system.components[0].lambdaGC.currentBin] = 1000.0;
  EXPECT_DOUBLE_EQ(system.weight(), 1.0);
}

TEST(lambda_reweighting, each_independent_histogram_normalizes_to_its_own_minimum)
{
  PropertyLambdaProbabilityHistogram first(5, 3);
  PropertyLambdaProbabilityHistogram second(5, 3);
  first.biasFactor = {2417.0, 2418.0, 2419.0};
  second.biasFactor = {0.0, 3.0, 7.0};

  first.normalizeToMinimum();
  second.normalizeToMinimum();

  EXPECT_DOUBLE_EQ(*std::ranges::min_element(first.biasFactor), 0.0);
  EXPECT_DOUBLE_EQ(*std::ranges::min_element(second.biasFactor), 0.0);
  EXPECT_DOUBLE_EQ(first.biasFactor[1], 1.0);
  EXPECT_DOUBLE_EQ(second.biasFactor[1], 3.0);
}

TEST(pressure_sampling, samples_the_freshly_computed_excess_tensor)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework framework = Framework::makeMFI(forceField, int3(1, 1, 1));
  Component component = Component::makeMethane(forceField, 0);
  System system =
      System(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, {framework}, {component}, {}, {0}, 5);

  const double3x3 excessPressure(4.0, 4.0, 4.0);
  system.sampleEnergyAndPressure(0, system.currentEnergyStatus, excessPressure);

  const PressureData sampled = system.averagePressure.averaged();
  EXPECT_DOUBLE_EQ(sampled.excessPressure, 4.0);
  EXPECT_DOUBLE_EQ(sampled.excessPressureTensor.trace(), 12.0);
}
