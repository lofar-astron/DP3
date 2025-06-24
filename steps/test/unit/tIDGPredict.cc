// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../../base/test/LoggerFixture.h"
#include "../../../common/ParameterSet.h"
#include <dp3/base/DPBuffer.h>
#include "../../../steps/InputStep.h"
#include "../../../steps/IDGPredict.h"
#include "../../../steps/MultiResultStep.h"

#include <schaapcommon/facets/facet.h>

#include <xtensor/xcomplex.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <memory>

using dp3::base::DPBuffer;
using dp3::steps::InputStep;
using dp3::steps::MultiResultStep;
using schaapcommon::facets::Facet;

namespace {
const unsigned int kNCorr = 4;
const unsigned int kNChan = 5;
const std::vector<std::size_t> kChannelCounts(kNChan, 1);
const unsigned int kStartChan = 0;
const double kFirstTime = 0.5;
const double kLastTime = 9.5;
const double kInterval = 1.0;
const std::size_t kTimeSteps = 5;
const std::string kMsName{};
const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2", "ant3"};
const std::vector<casacore::MPosition> kAntPos{
    casacore::MVPosition{0, 0, 0}, casacore::MVPosition{300, 0, 0},
    casacore::MVPosition{200, 0, 0}, casacore::MVPosition{2400, 0, 0}};
const std::vector<double> kAntDiam(4, 1.0);
const std::vector<int> kAnt1_1Bl{0, 0, 0};
const std::vector<int> kAnt2_1Bl{3, 1, 2};
const std::size_t kNBaselines = 3;

dp3::base::DPInfo InitInfo(const std::vector<int>& ant1,
                           const std::vector<int>& ant2,
                           std::size_t n_chan = kNChan) {
  BOOST_REQUIRE_EQUAL(ant1.size(), ant2.size());
  std::vector<double> chan_freqs(n_chan);
  std::vector<double> chan_widths(n_chan, 5000.0);
  for (std::size_t i = 0; i < n_chan; i++) {
    chan_freqs[i] = i * 10000.0;
  }

  dp3::base::DPInfo info(kNCorr, kNChan);
  info.setTimes(kFirstTime, kLastTime, kInterval);
  info.setAntennas(kAntNames, kAntDiam, kAntPos, ant1, ant2);
  info.setChannels(std::move(chan_freqs), std::move(chan_widths));
  return info;
}

/**
 * Create a buffer with artificial data values.
 * @param time Start time for the buffer.
 * @param interval Interval duration for the buffer.
 * @param n_baselines Number of baselines in the buffer.
 * @param base_value Base value for the data values, for distinguishing buffers.
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
std::unique_ptr<dp3::base::DPBuffer> CreateBuffer(
    const double time, const double interval, std::size_t n_baselines,
    const std::vector<std::size_t>& channel_counts, const float base_value,
    const float weight = 1.0) {
  const std::array<std::size_t, 3> kShape{n_baselines, channel_counts.size(),
                                          kNCorr};

  auto buffer = std::make_unique<dp3::base::DPBuffer>(time, interval);
  buffer->GetData().resize(kShape);
  buffer->GetWeights().resize(kShape);
  buffer->GetFlags().resize(kShape);
  buffer->GetUvw().resize({n_baselines, 3});

  buffer->GetFlags().fill(false);
  buffer->GetWeights().fill(weight);

  for (std::size_t baseline = 0; baseline < n_baselines; ++baseline) {
    // Base value for this baseline.
    const float baseline_value = (baseline * 100.0) + (base_value / weight);

    std::size_t channel = 0;
    float channel_value = baseline_value;  // Base value for a group of channels
    for (std::size_t channel_count : channel_counts) {
      // For each channel, increase channel_value by 10.0.
      // When channel_count == 1, 'value' should equal channel_value.
      // When channel_count > 1, 'value' should be the average for multiple
      // channels.
      const float value = channel_value + 5.0 * (channel_count - 1);
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        buffer->GetData()(baseline, channel, corr) = value + corr;
        buffer->GetWeights()(baseline, channel, corr) *= channel_count;
      }
      ++channel;
      channel_value += channel_count * 10.0;
    }
    buffer->GetUvw()(baseline, 0) = baseline_value + 0.0;
    buffer->GetUvw()(baseline, 1) = baseline_value + 1.0;
    buffer->GetUvw()(baseline, 2) = baseline_value + 2.0;
  }

  return buffer;
}

dp3::common::ParameterSet CreateParset() {
  dp3::common::ParameterSet parset;
  parset.add("regions", "sources.reg");
  parset.add("images", "sources-model.fits");
  return parset;
}

}  // namespace

using dp3::steps::IDGPredict;

BOOST_AUTO_TEST_SUITE(
    idgpredict, *boost::unit_test::fixture<dp3::base::test::LoggerFixture>())

BOOST_AUTO_TEST_CASE(getreaders, *boost::unit_test::tolerance(0.000001)) {
  const std::vector<aocommon::FitsReader> readers =
      IDGPredict::GetReaders({"sources-model.fits"});

  BOOST_TEST(readers.size() == 1u);
  BOOST_TEST(readers[0].PhaseCentreRA() == 0.426245723);
  BOOST_TEST(readers[0].PhaseCentreDec() == 0.578746973);
  BOOST_TEST(readers[0].ImageWidth() == 512u);
  BOOST_TEST(readers[0].ImageHeight() == 512u);
  BOOST_TEST(readers[0].PixelSizeX() == 0.000174533);
  BOOST_TEST(readers[0].PixelSizeY() == 0.000174533);
}

BOOST_AUTO_TEST_CASE(getfacets, *boost::unit_test::tolerance(1e-6)) {
  const std::vector<aocommon::FitsReader> readers =
      IDGPredict::GetReaders({"sources-model.fits"});

  const std::vector<Facet> facets =
      IDGPredict::GetFacets("sources.reg", readers.front());

  BOOST_TEST(facets.size() == 4u);
  BOOST_TEST(facets[0].RA() == 0.3877368);
  BOOST_TEST(facets[0].Dec() == 0.5417391);

  const std::vector<schaapcommon::facets::PixelPosition>& p0 =
      facets[0].GetPixels();
  BOOST_TEST_REQUIRE(p0.size() == 4u);
  BOOST_TEST(p0[0].x == 378);
  BOOST_TEST(p0[0].y == 0);
  BOOST_TEST(p0[1].x == 377);
  BOOST_TEST(p0[1].y == 91);
  BOOST_TEST(p0[2].x == 512);
  BOOST_TEST(p0[2].y == 94);
  BOOST_TEST(p0[3].x == 512);
  BOOST_TEST(p0[3].y == 0);

  const std::vector<schaapcommon::facets::PixelPosition>& p1 =
      facets[1].GetPixels();
  BOOST_TEST_REQUIRE(p1.size() == 4u);
  BOOST_TEST(p1[0].x == 512);
  BOOST_TEST(p1[0].y == 304);
  BOOST_TEST(p1[1].x == 512);
  BOOST_TEST(p1[1].y == 128);
  BOOST_TEST(p1[2].x == 364);
  BOOST_TEST(p1[2].y == 129);
  BOOST_TEST(p1[3].x == 366);
  BOOST_TEST(p1[3].y == 301);
}

BOOST_AUTO_TEST_CASE(constructor) {
  dp3::common::ParameterSet parset;
  const std::vector<aocommon::FitsReader> readers =
      IDGPredict::GetReaders({"sources-model.fits"});

  std::vector<Facet> facets =
      IDGPredict::GetFacets("sources.reg", readers.front());

  IDGPredict predict(parset, "", readers, std::move(facets));
  predict.SetBufferSize(42);

  BOOST_TEST(predict.IsStarted() == false);
  BOOST_TEST(predict.GetBufferSize() == 42u);

  // Check that the print methods do not throw any error
  std::ostream nullout(nullptr);
  BOOST_CHECK_NO_THROW(predict.show(nullout));
  BOOST_CHECK_NO_THROW(predict.showTimings(nullout, 10));
}

BOOST_AUTO_TEST_CASE(lean_constructor) {
  const dp3::common::ParameterSet parset = CreateParset();

  IDGPredict predict(parset, "");

  BOOST_TEST(predict.IsStarted() == false);
  BOOST_TEST(predict.GetBufferSize() == 0u);
}

BOOST_AUTO_TEST_CASE(no_regions) {
  dp3::common::ParameterSet parset;
  parset.add("regions", "im-updside-down.reg");
  parset.add("images", "sources-model.fits");

  BOOST_CHECK_THROW(std::make_unique<IDGPredict>(parset, ""),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(no_models) {
  dp3::common::ParameterSet parset;
  parset.add("regions", "sources.reg");

  BOOST_CHECK_THROW(std::make_unique<IDGPredict>(parset, ""),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(update_info_wrong) {
  const dp3::common::ParameterSet parset = CreateParset();

  dp3::base::DPInfo info;

  IDGPredict predict(parset, "");

  BOOST_CHECK_THROW(predict.setInfo(info), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(update_info) {
  const dp3::common::ParameterSet parset = CreateParset();

  const dp3::base::DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);

  IDGPredict predict(parset, "");

  predict.setInfo(info);

  BOOST_TEST(predict.IsStarted() == true);
  BOOST_TEST(predict.GetBufferSize() > 0u);
}

BOOST_AUTO_TEST_CASE(process, *boost::unit_test::tolerance(0.1f) *
                                  boost::unit_test::label("slow")) {
  const dp3::common::ParameterSet parset = CreateParset();

  IDGPredict predict(parset, "");
  predict.SetBufferSize(kTimeSteps);

  auto result_step = std::make_shared<MultiResultStep>(kTimeSteps);

  predict.setNextStep(result_step);

  const dp3::base::DPInfo info = InitInfo(kAnt1_1Bl, kAnt2_1Bl);
  predict.setInfo(info);

  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    std::unique_ptr<DPBuffer> buffer =
        CreateBuffer(kFirstTime + i * kInterval, kInterval, kNBaselines,
                     kChannelCounts, i * 1000.0);
    predict.process(std::move(buffer));
  }
  predict.finish();

  BOOST_TEST(result_step->size() == kTimeSteps);
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    BOOST_TEST(xt::real(result_step->get()[i]->GetData()(0, 0, 0)) == 60.f);
    BOOST_TEST(xt::imag(result_step->get()[i]->GetData()(0, 0, 0)) == 0.f);
  }
}

BOOST_AUTO_TEST_CASE(process_beam, *boost::unit_test::tolerance(0.0001f) *
                                       boost::unit_test::label("slow")) {
  dp3::common::ParameterSet parset;
  parset.add("msin", "tNDPPP-generic.MS");
  parset.add("aterms", "beam");
  parset.add("beam.element_response_model", "hamaker");

  const std::vector<aocommon::FitsReader> fits_readers =
      IDGPredict::GetReaders({"sources-model.fits"});
  std::vector<Facet> facets =
      IDGPredict::GetFacets("sources.reg", fits_readers.front());

  // Set RA and Dec pointing attached to the facets to 0, since kExpectedData
  // was generated for these facet pointings - rather than the facet centroids.
  for (auto& facet : facets) {
    facet.SetRA(0.0);
    facet.SetDec(0.0);
  }

  // This test needs a real reader, since IDGPredict passes a MeasurementSet to
  // Everybeam::Load when using aterms. MockInput does not suffice.
  std::shared_ptr<InputStep> reader = InputStep::CreateReader(parset);
  auto predict =
      std::make_shared<IDGPredict>(parset, "", fits_readers, std::move(facets));
  predict->SetBufferSize(kTimeSteps);
  auto result_step = std::make_shared<MultiResultStep>(kTimeSteps);

  reader->setNextStep(predict);
  predict->setNextStep(result_step);

  reader->setInfo(dp3::base::DPInfo());
  reader->setFieldsToRead(predict->getRequiredFields());

  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    reader->process(std::make_unique<DPBuffer>());
  }
  reader->finish();

  // These samples were gathered while running this test initially.
  // Ideally this test should use a mocked aterm, and check that IDGPredict
  // calls the proper functions on the mock. For now, this approach will do.
  const std::complex<float> kExpectedData[kTimeSteps][3] = {
      {{0.00020892531, 0.000115041898},
       {-0.000436737319, -0.000309092226},
       {0.00273089739, 0.00102706277}},
      {{0.000207837016, 0.000117528049},
       {-0.000436606992, -0.000309304567},
       {0.00272684987, 0.00102799805}},
      {{0.000206714481, 0.000119998658},
       {-0.000436475908, -0.000309517025},
       {0.00272280071, 0.00102892378}},
      {{0.000205557706, 0.000122453159},
       {-0.000436344679, -0.000309729861},
       {0.00271875248, 0.00102983904}},
      {{0.000204367054, 0.00012489123},
       {-0.000436212751, -0.00030994293},
       {0.00271470915, 0.00103074359}}};

  BOOST_TEST(result_step->size() == kTimeSteps);
  for (std::size_t i = 0; i < kTimeSteps; ++i) {
    const auto& data = result_step->get()[i]->GetData();
    const size_t n_bl = data.shape(0);
    const size_t n_chan = data.shape(1);
    const size_t n_corr = data.shape(2);

    BOOST_REQUIRE_EQUAL(n_bl, reader->getInfoOut().nbaselines());
    BOOST_REQUIRE_EQUAL(n_chan, reader->getInfoOut().nchan());
    BOOST_REQUIRE_EQUAL(n_corr, reader->getInfoOut().ncorr());

    // Take samples of the result values and compare those.
    const std::complex<float>& data0 = data(0, 0, 0);
    const std::complex<float>& data1 = data(1, 6, 11);
    const std::complex<float>& data2 = data(n_corr - 1, n_chan - 1, n_bl - 1);
    BOOST_TEST(data0.real() == kExpectedData[i][0].real());
    BOOST_TEST(data0.imag() == kExpectedData[i][0].imag());
    BOOST_TEST(data1.real() == kExpectedData[i][1].real());
    BOOST_TEST(data1.imag() == kExpectedData[i][1].imag());
    BOOST_TEST(data2.real() == kExpectedData[i][2].real());
    BOOST_TEST(data2.imag() == kExpectedData[i][2].imag());
  }
}

BOOST_AUTO_TEST_SUITE_END()
