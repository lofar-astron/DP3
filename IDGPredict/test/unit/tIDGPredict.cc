#include "../../IDGPredict.h"
#include "../../../Common/ParameterSet.h"
#include "../../../DPPP/DPBuffer.h"
#include "../../../DPPP/DPInput.h"
#include "../../../DPPP/test/unit/mock/MockInput.h"

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using DP3::DPPP::DPBuffer;
using DP3::DPPP::MultiResultStep;

namespace {
const unsigned int kNCorr = 4;
const unsigned int kNChan = 5;
const std::vector<std::size_t> kChannelCounts(kNChan, 1);
const unsigned int kStartChan = 0;
const unsigned int kNTime = 10;
const double kStartTime = 0.0;
const double kInterval = 1.0;
const std::size_t kTimeSteps = 5;
const std::string kMsName{};
const std::string kAntennaSet{};
const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2", "ant3"};
const std::vector<casacore::MPosition> kAntPos{
    casacore::MVPosition{0, 0, 0}, casacore::MVPosition{300, 0, 0},
    casacore::MVPosition{200, 0, 0}, casacore::MVPosition{2400, 0, 0}};
const std::vector<double> kAntDiam(4, 1.0);
const std::vector<int> kAnt1_1Bl{0, 0, 0};
const std::vector<int> kAnt2_1Bl{3, 1, 2};
const std::size_t kNBaselines = 3;

void InitInfo(DP3::DPPP::DPInfo& info, const std::vector<int>& ant1,
              const std::vector<int>& ant2, std::size_t n_chan = kNChan) {
  BOOST_REQUIRE_EQUAL(ant1.size(), ant2.size());
  std::vector<double> chan_freqs(n_chan);
  std::vector<double> chan_widths(n_chan, 5000.0);
  for (std::size_t i = 0; i < n_chan; i++) {
    chan_freqs[i] = i * 10000.0;
  }

  info.init(kNCorr, kStartChan, n_chan, kNTime, kStartTime, kInterval, kMsName,
            kAntennaSet);
  info.set(kAntNames, kAntDiam, kAntPos, ant1, ant2);
  info.set(std::move(chan_freqs), std::move(chan_widths));
}

/**
 * Create a buffer with artifical data values.
 * @param time Start time for the buffer.
 * @param interval Interval duration for the buffer.
 * @param n_baselines Number of baselines in the buffer.
 * @param base_value Base value for the data values, for distinguising buffers.
 *        For distinguishing baselines, this function adds baseline_nr * 100.0.
 *        When the buffer represents averaged data, the base_value should be
 *        the total of the base values of the original buffers.
 *        This function divides the base_value by the supplied weight so the
 *        caller does not have to do that division.
 * @param channel_counts List for generating channel data.
 *        For input buffers, this list should contain a 1 for each channel.
 *        When generating expected output data, this list should contain the
 *        number of averaged input buffers for each output buffer.
 * @param weight Weight value for the data values in the buffer.
 */
std::unique_ptr<DPBuffer> CreateBuffer(
    const double time, const double interval, std::size_t n_baselines,
    const std::vector<std::size_t>& channel_counts, const float base_value,
    const float weight = 1.0) {
  casacore::Cube<casacore::Complex> data(kNCorr, channel_counts.size(),
                                         n_baselines);
  casacore::Cube<bool> flags(data.shape(), false);
  casacore::Cube<float> weights(data.shape(), weight);
  casacore::Cube<bool> full_res_flags(channel_counts.size(), 1, n_baselines,
                                      false);
  casacore::Matrix<double> uvw(3, n_baselines);
  for (std::size_t bl = 0; bl < n_baselines; ++bl) {
    // Base value for this baseline.
    const float bl_value = (bl * 100.0) + (base_value / weight);

    std::size_t chan = 0;
    float chan_value = bl_value;  // Base value for a group of channels.
    for (std::size_t ch_count : channel_counts) {
      // For each channel, increase chan_base by 10.0.
      // When ch_count == 1, 'value' should equal chan_base.
      // When ch_count > 1, 'value' should be the average for multiple channels.
      const float value = chan_value + 5.0 * (ch_count - 1);
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        data(corr, chan, bl) = value + corr;
        weights(corr, chan, bl) *= ch_count;
      }
      ++chan;
      chan_value += ch_count * 10.0;
    }
    uvw(0, bl) = bl_value + 0.0;
    uvw(1, bl) = bl_value + 1.0;
    uvw(2, bl) = bl_value + 2.0;
  }

  auto buffer = boost::make_unique<DPBuffer>();
  buffer->setTime(time);
  buffer->setExposure(interval);
  buffer->setData(data);
  buffer->setWeights(weights);
  buffer->setFlags(flags);
  buffer->setFullResFlags(full_res_flags);
  buffer->setUVW(uvw);

  return buffer;
}

DP3::DPPP::MockInput mock_input;
}  // namespace

using aocommon::UVector;
using DP3::DPPP::IDGPredict;

BOOST_AUTO_TEST_SUITE(idgpredict)

BOOST_AUTO_TEST_CASE(getreaders, *boost::unit_test::tolerance(0.000001)) {
  std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<double>>>
      readers = IDGPredict::GetReaders({"sources-model.fits"});

  BOOST_TEST(readers.first.size() == 1u);
  BOOST_TEST(readers.first[0].PhaseCentreRA() == 0.426245723);
  BOOST_TEST(readers.first[0].PhaseCentreDec() == 0.578746973);
  BOOST_TEST(readers.first[0].ImageWidth() == 512u);
  BOOST_TEST(readers.first[0].ImageHeight() == 512u);
  BOOST_TEST(readers.first[0].PixelSizeX() == 0.000174533);
  BOOST_TEST(readers.first[0].PixelSizeY() == 0.000174533);
}

BOOST_AUTO_TEST_CASE(getfacets) {
  std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<double>>>
      readers = IDGPredict::GetReaders({"sources-model.fits"});

  std::vector<Facet> facets =
      IDGPredict::GetFacets("sources.reg", readers.first.front());

  BOOST_TEST(facets.size() == 4u);
  BOOST_TEST(facets[0].RA() == 0.);
  BOOST_TEST(facets[0].Dec() == 0.);

  auto it = facets[0].begin();
  BOOST_TEST(it->x == 379);
  BOOST_TEST(it->y == -9);
  BOOST_TEST((it + 1)->x == 377);
  BOOST_TEST((it + 1)->y == 91);
  BOOST_TEST((it + 2)->x == 548);
  BOOST_TEST((it + 2)->y == 95);
  BOOST_TEST((it + 3)->x == 551);
  BOOST_TEST((it + 3)->y == -5);

  it = facets[1].begin();
  BOOST_TEST(it->x == 364);
  BOOST_TEST(it->y == 129);
  BOOST_TEST((it + 1)->x == 366);
  BOOST_TEST((it + 1)->y == 301);
  BOOST_TEST((it + 2)->x == 541);
  BOOST_TEST((it + 2)->y == 305);
  BOOST_TEST((it + 3)->x == 547);
  BOOST_TEST((it + 3)->y == 128);
}

BOOST_AUTO_TEST_CASE(constructor) {
  DP3::ParameterSet parset;
  std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<double>>>
      readers = IDGPredict::GetReaders({"sources-model.fits"});

  std::vector<Facet> facets =
      IDGPredict::GetFacets("sources.reg", readers.first.front());

  IDGPredict predict(mock_input, parset, "", readers, std::move(facets));
  predict.SetBufferSize(42);

  BOOST_TEST(predict.IsStarted() == false);
  BOOST_TEST(predict.GetBufferSize() == 42u);

  // Check that the print methods do not throw any error
  std::ostream nullout(nullptr);
  BOOST_CHECK_NO_THROW(predict.show(nullout));
  BOOST_CHECK_NO_THROW(predict.showTimings(nullout, 10));
}

BOOST_AUTO_TEST_CASE(lean_constructor) {
  DP3::ParameterSet parset;
  parset.add("regions", "sources.reg");
  parset.add("images", {"sources-model.fits"});

  IDGPredict predict(mock_input, parset, "");

  BOOST_TEST(predict.IsStarted() == false);
  BOOST_TEST(predict.GetBufferSize() == 0u);
}

BOOST_AUTO_TEST_CASE(no_regions) {
  DP3::ParameterSet parset;
  parset.add("regions", "im-updside-down.reg");
  parset.add("images", {"sources-model.fits"});

  BOOST_CHECK_THROW(new IDGPredict(mock_input, parset, ""), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(no_models) {
  DP3::ParameterSet parset;
  parset.add("regions", "sources.reg");
  parset.add("images", {});

  BOOST_CHECK_THROW(new IDGPredict(mock_input, parset, ""), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(update_info_wrong) {
  DP3::ParameterSet parset;
  parset.add("regions", "sources.reg");
  parset.add("images", {"sources-model.fits"});

  DP3::DPPP::DPInfo info;

  IDGPredict predict(mock_input, parset, "");

  BOOST_CHECK_THROW(predict.setInfo(info), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(update_info) {
  DP3::ParameterSet parset;
  parset.add("regions", "sources.reg");
  parset.add("images", {"sources-model.fits"});

  DP3::DPPP::DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);

  IDGPredict predict(mock_input, parset, "");

  predict.setInfo(info);

  BOOST_TEST(predict.IsStarted() == true);
  BOOST_TEST(predict.GetBufferSize() > 0u);
}

BOOST_AUTO_TEST_CASE(process, *boost::unit_test::tolerance(0.1f)) {
  DP3::ParameterSet parset;
  parset.add("regions", "sources.reg");
  parset.add("images", {"sources-model.fits"});

  DP3::DPPP::DPInfo info;
  InitInfo(info, kAnt1_1Bl, kAnt2_1Bl);

  std::vector<std::unique_ptr<DPBuffer>> buffers;
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    buffers.push_back(CreateBuffer(kStartTime + i * kInterval, kInterval,
                                   kNBaselines, kChannelCounts, i * 1000.0));
  }

  auto result_step = std::make_shared<MultiResultStep>(kTimeSteps);

  IDGPredict predict(mock_input, parset, "");
  predict.setNextStep(result_step);
  predict.SetBufferSize(kTimeSteps);
  predict.setInfo(info);
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    predict.process(*buffers[i]);
  }
  predict.finish();

  BOOST_TEST(result_step.get()->size() == kTimeSteps);
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(result_step->get()[i].getData().data()->imag() == 0.f);
    BOOST_TEST(result_step->get()[i].getData().data()->real() == 60.f);
  }
}

BOOST_AUTO_TEST_SUITE_END()
