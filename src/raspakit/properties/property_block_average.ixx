module;

export module property_block_average;

import std;

import archive;
import averages;

/**
 * \brief Requirements for a sample type that can be block-averaged directly.
 *
 * The per-block accumulators require weighted accumulation (\c += and \c scalar*sample) and
 * normalisation by the accumulated weight (\c sample/scalar). The value itself must be
 * \ref BlockAverageable so its confidence-interval error can be estimated.
 */
export template <typename Sample>
concept DirectlyAverageable = BlockAverageable<Sample> && requires(Sample sample, double scalar) {
  { sample / scalar } -> std::convertible_to<Sample>;
};

template <typename Sample>
struct BlockAverageDataType
{
  using type = Sample;
};

template <CompositePropertyTerms Sample>
struct BlockAverageDataType<Sample>
{
  using type = std::remove_cvref_t<decltype(std::declval<Sample>().compositeProperty())>;
};

/**
 * \brief Generic block-averaging accumulator for simulation properties.
 *
 * Owns the per-block book-keeping (accumulated weighted samples plus accumulated weights) and
 * provides the complete statistical pipeline: per-block averages, the overall average, and the
 * confidence-interval error estimated from the spread of the block averages.
 *
 * Two families of properties are supported through the same interface:
 * - **Direct** properties (\ref DirectlyAverageable): the block average is the normalised
 *   accumulator itself and the overall average pools all blocks (weighted mean).
 * - **Composite** properties (\ref CompositePropertyTerms): raw linear fluctuation terms are
 *   accumulated; each block is reduced through \c compositeProperty() and the overall average is
 *   the mean of the sufficiently-filled block composites.
 *
 * Non-linear derived quantities (e.g. chemical potential as \f$-k_BT\ln\langle W\rangle\f$) are
 * supported through statistics(), which propagates any transformation of the block averages
 * through the same error estimate.
 *
 * \tparam Sample The accumulated per-sample type (a data struct, \c double, or fluctuation terms).
 */
export template <typename Sample>
  requires DirectlyAverageable<Sample> || CompositePropertyTerms<Sample>
struct BlockAverage
{
  using sample_type = Sample;

  /// The type averages and errors are expressed in: Sample itself for direct properties,
  /// the result of Sample::compositeProperty() for composite properties.
  using Data = typename BlockAverageDataType<Sample>::type;

  std::uint64_t versionNumber{2};
  std::size_t numberOfBlocks{};
  Sample zeroSample{};
  std::vector<std::pair<Sample, double>> bookKeeping;
  std::vector<std::size_t> rawSampleCounts;

  BlockAverage() = default;

  explicit BlockAverage(std::size_t numberOfBlocks, Sample zero = Sample{})
      : numberOfBlocks(numberOfBlocks),
        zeroSample(zero),
        bookKeeping(numberOfBlocks, std::make_pair(zero, 0.0)),
        rawSampleCounts(numberOfBlocks)
  {
  }

  bool operator==(BlockAverage const &) const = default;

  /// Rebuild the book-keeping for a new correctly-shaped zero sample.
  void resize(Sample zero)
  {
    zeroSample = zero;
    bookKeeping = std::vector<std::pair<Sample, double>>(numberOfBlocks, std::make_pair(zero, 0.0));
    rawSampleCounts.assign(numberOfBlocks, 0);
  }

  void addSample(std::size_t blockIndex, const Sample &sample, const double &weight)
  {
    if (!std::isfinite(weight) || weight < 0.0)
    {
      throw std::invalid_argument("BlockAverage sample weight must be finite and non-negative");
    }
    ++rawSampleCounts[blockIndex];
    if (weight == 0.0)
    {
      return;
    }

    bookKeeping[blockIndex].first += weight * sample;
    bookKeeping[blockIndex].second += weight;
  }

  /// Average of a single block.
  Data averaged(std::size_t blockIndex) const
  {
    const double accumulatedWeight = bookKeeping[blockIndex].second;
    if (accumulatedWeight == 0.0)
    {
      if (rawSampleCounts[blockIndex] > 0)
      {
        throw std::runtime_error("All BlockAverage sample weights in a non-empty block are zero");
      }
      if constexpr (CompositePropertyTerms<Sample>)
      {
        return zeroSample.zeroCompositeProperty();
      }
      else
      {
        return zeroSample;
      }
    }

    Sample normalized = bookKeeping[blockIndex].first / accumulatedWeight;
    if constexpr (CompositePropertyTerms<Sample>)
    {
      return normalized.compositeProperty();
    }
    else
    {
      return normalized;
    }
  }

  /// Overall average: pooled weighted mean (direct) or mean of the sufficiently-filled block
  /// composites (composite).
  Data averaged() const
  {
    if constexpr (CompositePropertyTerms<Sample>)
    {
      Data accumulated = zeroSample.zeroCompositeProperty();
      std::size_t numberOfSamples = 0;
      const std::size_t reference =
          rawSampleCounts.empty() ? 0 : *std::ranges::max_element(rawSampleCounts);
      for (std::size_t blockIndex = 0; blockIndex != bookKeeping.size(); ++blockIndex)
      {
        if (reference > 0 && rawSampleCounts[blockIndex] > reference / 2)
        {
          accumulated += averaged(blockIndex);
          ++numberOfSamples;
        }
      }
      return (1.0 / static_cast<double>(std::max<std::size_t>(1, numberOfSamples))) * accumulated;
    }
    else
    {
      std::pair<Sample, double> summedBlocks =
          std::accumulate(bookKeeping.begin(), bookKeeping.end(), std::make_pair(zeroSample, 0.0),
                          pairSum<Sample, double>);
      if (summedBlocks.second > 0.0)
      {
        return summedBlocks.first / summedBlocks.second;
      }
      if (std::ranges::any_of(rawSampleCounts, [](std::size_t count) { return count > 0; }))
      {
        throw std::runtime_error("All BlockAverage sample weights are zero");
      }
      return zeroSample;
    }
  }

  /// Mean and confidence-interval error of any (possibly non-linear) function of the averages.
  /// The mean is the transform of the overall average; the error is estimated from the spread of
  /// the transformed block averages.
  template <typename Transform>
  auto statistics(Transform transform) const
  {
    auto mean = transform(averaged());
    auto confidenceIntervalError = blockErrorEstimate(
        bookKeeping, rawSampleCounts, mean, [&](std::size_t i) { return transform(averaged(i)); });
    return std::make_pair(mean, confidenceIntervalError);
  }

  /// Overall average together with its confidence-interval error.
  std::pair<Data, Data> average() const
  {
    return statistics([](const Data &data) { return data; });
  }
};

/**
 * \brief Generic block-averaged histogram for binned properties.
 *
 * Owns per-block, per-channel bin accumulators together with a per-block sample count that is
 * shared by all channels. Samples are accumulated directly into individual bins through
 * operator() (an O(1) increment; histogram sampling touches single bins rather than whole
 * vectors), and addCount() registers the sample weight once per sample. Provides the complete
 * statistical pipeline: per-block averages (bins normalised by the block count), the pooled
 * overall average, and the per-bin confidence-interval error estimated from the spread of the
 * block averages (blocks without samples are skipped).
 *
 * Channels distinguish independent histograms that share the same sampling (per component,
 * per pseudo-atom pair, ...); single-histogram properties use one channel and index 0.
 *
 * \tparam Value The per-bin value type (\c double or a data struct such as AverageEnergyType).
 */
export template <DirectlyAverageable Value>
struct BlockHistogram
{
  std::uint64_t versionNumber{2};
  std::size_t numberOfBlocks{};
  std::size_t numberOfChannels{};
  std::size_t numberOfBins{};
  /// Accumulated bins indexed as [block][channel][bin].
  std::vector<std::vector<std::vector<Value>>> bookKeeping;
  /// Accumulated sample weight per block, shared by all channels.
  std::vector<double> numberOfCounts;
  std::vector<std::size_t> rawSampleCounts;
  double totalNumberOfCounts{};

  BlockHistogram() = default;

  BlockHistogram(std::size_t numberOfBlocks, std::size_t numberOfChannels, std::size_t numberOfBins)
      : numberOfBlocks(numberOfBlocks),
        numberOfChannels(numberOfChannels),
        numberOfBins(numberOfBins),
        bookKeeping(numberOfBlocks,
                    std::vector<std::vector<Value>>(numberOfChannels, std::vector<Value>(numberOfBins))),
        numberOfCounts(numberOfBlocks),
        rawSampleCounts(numberOfBlocks)
  {
  }

  bool operator==(BlockHistogram const &) const = default;

  /// Direct accumulation access to a single bin.
  Value &operator()(std::size_t blockIndex, std::size_t channelIndex, std::size_t binIndex)
  {
    return bookKeeping[blockIndex][channelIndex][binIndex];
  }

  /// Register a sample weight (once per sample, not per channel).
  void addCount(std::size_t blockIndex, double weight)
  {
    if (!std::isfinite(weight) || weight < 0.0)
    {
      throw std::invalid_argument("BlockHistogram sample weight must be finite and non-negative");
    }
    ++rawSampleCounts[blockIndex];
    if (weight == 0.0)
    {
      return;
    }

    numberOfCounts[blockIndex] += weight;
    totalNumberOfCounts += weight;
  }

  /// Average histogram of a single block.
  std::vector<Value> averaged(std::size_t blockIndex, std::size_t channelIndex) const
  {
    std::vector<Value> averagedData(numberOfBins);
    if (numberOfCounts[blockIndex] == 0.0)
    {
      if (rawSampleCounts[blockIndex] > 0)
      {
        throw std::runtime_error("All BlockHistogram sample weights in a non-empty block are zero");
      }
      return averagedData;
    }

    std::transform(bookKeeping[blockIndex][channelIndex].begin(), bookKeeping[blockIndex][channelIndex].end(),
                   averagedData.begin(),
                   [&](const Value &sample) { return sample / numberOfCounts[blockIndex]; });
    return averagedData;
  }

  /// Overall average histogram, pooled over all blocks.
  std::vector<Value> averaged(std::size_t channelIndex) const
  {
    std::vector<Value> summedBlocks(numberOfBins);
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      for (std::size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
      {
        summedBlocks[binIndex] += bookKeeping[blockIndex][channelIndex][binIndex];
      }
    }

    std::vector<Value> average(numberOfBins);
    if (totalNumberOfCounts == 0.0)
    {
      if (std::ranges::any_of(rawSampleCounts, [](std::size_t count) { return count > 0; }))
      {
        throw std::runtime_error("All BlockHistogram sample weights are zero");
      }
      return average;
    }

    std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                   [&](const Value &sample) { return sample / totalNumberOfCounts; });
    return average;
  }

  /// Overall average histogram together with its per-bin confidence-interval error.
  std::pair<std::vector<Value>, std::vector<Value>> average(std::size_t channelIndex) const
  {
    std::vector<Value> mean = averaged(channelIndex);

    std::vector<std::vector<Value>> blockAverages;
    blockAverages.reserve(numberOfBlocks);
    const std::size_t reference =
        rawSampleCounts.empty() ? 0 : *std::ranges::max_element(rawSampleCounts);
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (reference > 0 && rawSampleCounts[blockIndex] > reference / 2)
      {
        blockAverages.push_back(averaged(blockIndex, channelIndex));
      }
    }

    std::vector<Value> confidenceIntervalError = blockErrorEstimate(blockAverages, mean);

    return std::make_pair(mean, confidenceIntervalError);
  }
};

export template <typename Value>
Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BlockHistogram<Value> &h)
{
  archive << h.versionNumber;
  archive << h.numberOfBlocks;
  archive << h.numberOfChannels;
  archive << h.numberOfBins;
  archive << h.bookKeeping;
  archive << h.numberOfCounts;
  archive << h.rawSampleCounts;
  archive << h.totalNumberOfCounts;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

export template <typename Value>
Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BlockHistogram<Value> &h)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > h.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BlockHistogram' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> h.numberOfBlocks;
  archive >> h.numberOfChannels;
  archive >> h.numberOfBins;
  archive >> h.bookKeeping;
  archive >> h.numberOfCounts;
  if (versionNumber >= 2)
  {
    archive >> h.rawSampleCounts;
  }
  else
  {
    h.rawSampleCounts.resize(h.numberOfBlocks);
    std::ranges::transform(h.numberOfCounts, h.rawSampleCounts.begin(),
                           [](double weight) { return weight > 0.0 ? 1uz : 0uz; });
  }
  archive >> h.totalNumberOfCounts;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("BlockHistogram: Error in binary restart\n"));
  }
#endif

  return archive;
}

export template <typename Property>
  requires std::derived_from<Property, BlockAverage<typename Property::sample_type>>
Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Property &p)
{
  archive << p.versionNumber;
  archive << p.numberOfBlocks;
  archive << p.zeroSample;
  archive << p.bookKeeping;
  archive << p.rawSampleCounts;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

export template <typename Property>
  requires std::derived_from<Property, BlockAverage<typename Property::sample_type>>
Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Property &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'BlockAverage' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.numberOfBlocks;
  archive >> p.zeroSample;
  archive >> p.bookKeeping;
  if (versionNumber >= 2)
  {
    archive >> p.rawSampleCounts;
  }
  else
  {
    p.rawSampleCounts.resize(p.numberOfBlocks);
    std::ranges::transform(p.bookKeeping, p.rawSampleCounts.begin(),
                           [](const auto &entry) { return entry.second > 0.0 ? 1uz : 0uz; });
  }

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("BlockAverage: Error in binary restart\n"));
  }
#endif

  return archive;
}
