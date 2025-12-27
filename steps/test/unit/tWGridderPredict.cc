// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <memory>

#include <aocommon/fits/fitsreader.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <schaapcommon/facets/facet.h>
#include <xtensor/xcomplex.hpp>

#include "base/DPBuffer.h"

#include "../../../base/test/LoggerFixture.h"
#include "../../../common/ParameterSet.h"
#include "../../../steps/InputStep.h"
#include "../../../steps/WGridderPredict.h"
#include "../../../steps/MultiResultStep.h"
#include "../../../steps/OnePredict.h"

using dp3::base::DPBuffer;
using dp3::steps::InputStep;
using dp3::steps::MultiResultStep;
using dp3::steps::OnePredict;
using dp3::steps::WGridderPredict;
using schaapcommon::facets::Facet;

namespace {
const unsigned int kNCorrelations = 4;
const unsigned int kNChannels = 5;
const double kStartFrequency = 150.0e6;
const double kFrequencyStep = 100.0e3;
const double kFirstTime = 0.5;
const double kLastTime = 9.5;
const double kInterval = 1.0;
const std::size_t kNTimeSteps = 5;
const std::string kMsName{};
const std::vector<std::string> kAntennaNames{"ant0", "ant1", "ant2", "ant3"};
const std::vector<casacore::MPosition> kAntennaPositions{
    casacore::MVPosition{0, 0, 0}, casacore::MVPosition{300, 0, 0},
    casacore::MVPosition{200, 0, 0}, casacore::MVPosition{2400, 0, 0}};
const std::vector<double> kAntennaDiameters(4, 1.0);
const std::vector<int> kAntennas1{0, 0, 0};
const std::vector<int> kAntennas2{3, 1, 2};
const std::size_t kNBaselines = 3;
// Distance between generated uvw points is 1.0 km
const double kBaselineStepSize = 1.0e3;
// Do not use relative tolerance for visibilities
const double kVisibilitiesRelativeTolerance = 0.0;
// Gridded visibilities are an approximation
// Tolerance needs to be a bit higher than in the case of an exact precit
const double kVisibilitiesAbsoluteTolerance = 1.0e-3;

dp3::base::DPInfo InitInfo() {
  std::vector<double> channel_frequencies(kNChannels);
  std::vector<double> channel_widths(kNChannels, 5000.0);
  for (std::size_t i = 0; i < kNChannels; i++) {
    channel_frequencies[i] = kStartFrequency + i * kFrequencyStep;
  }

  dp3::base::DPInfo info(kNCorrelations, kNChannels);
  info.setTimes(kFirstTime, kLastTime, kInterval);
  info.setAntennas(kAntennaNames, kAntennaDiameters, kAntennaPositions,
                   kAntennas1, kAntennas2);
  info.setChannels(std::move(channel_frequencies), std::move(channel_widths));
  return info;
}

/**
 * Create a buffer with artificial data values.
 * @param time Start time for the buffer.
 * @param interval Interval duration for the buffer.
 */
std::unique_ptr<dp3::base::DPBuffer> CreateBuffer(const double time,
                                                  const double interval) {
  auto buffer = std::make_unique<dp3::base::DPBuffer>(time, interval);
  buffer->GetUvw().resize({kNBaselines, 3});

  buffer->GetFlags().fill(false);

  // Create different arbitrary uvw coordinates
  for (std::size_t baseline = 0; baseline < kNBaselines; ++baseline) {
    buffer->GetUvw()(baseline, 0) = kBaselineStepSize * (baseline % 10);
    buffer->GetUvw()(baseline, 1) = kBaselineStepSize * (baseline / 10);
    buffer->GetUvw()(baseline, 2) = kBaselineStepSize;
  }

  return buffer;
}

dp3::common::ParameterSet CreateParset() {
  dp3::common::ParameterSet parset;
  parset.add("regions", "sources.reg");
  parset.add("images", "sources-model.fits");
  parset.add("sumfacets", "true");
  return parset;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(
    wgridderpredict,
    *boost::unit_test::fixture<dp3::base::test::LoggerFixture>())

BOOST_AUTO_TEST_CASE(getreaders, *boost::unit_test::tolerance(1.0e-6)) {
  const std::vector<aocommon::FitsReader> readers =
      WGridderPredict::GetReaders({"sources-model.fits"});

  BOOST_TEST_REQUIRE(readers.size() == 1u);
  BOOST_TEST(readers[0].PhaseCentreRA() == 0.426245723);
  BOOST_TEST(readers[0].PhaseCentreDec() == 0.578746973);
  BOOST_TEST(readers[0].ImageWidth() == 512u);
  BOOST_TEST(readers[0].ImageHeight() == 512u);
  BOOST_TEST(readers[0].PixelSizeX() == 0.000174533);
  BOOST_TEST(readers[0].PixelSizeY() == 0.000174533);
}

BOOST_AUTO_TEST_CASE(getfacets, *boost::unit_test::tolerance(1.0e-6)) {
  const std::vector<aocommon::FitsReader> readers =
      WGridderPredict::GetReaders({"sources-model.fits"});

  std::vector<Facet> facets =
      WGridderPredict::GetFacets("sources.reg", readers.front());

  BOOST_TEST_REQUIRE(facets.size() == 4u);
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
  std::vector<aocommon::FitsReader> readers =
      WGridderPredict::GetReaders({"sources-model.fits"});

  std::vector<Facet> facets =
      WGridderPredict::GetFacets("sources.reg", readers.front());

  WGridderPredict predict(parset, "", std::move(readers), std::move(facets));

  BOOST_TEST(predict.GetBufferSize() == 0u);
  predict.SetBufferSize(42);
  BOOST_TEST(predict.GetBufferSize() == 42u);

  // Check that the show methods do not throw any error.
  std::ostream nullout(nullptr);
  BOOST_CHECK_NO_THROW(predict.show(nullout));
  BOOST_CHECK_NO_THROW(predict.showTimings(nullout, 10));
}

BOOST_AUTO_TEST_CASE(lean_constructor) {
  const dp3::common::ParameterSet parset = CreateParset();

  const WGridderPredict predict(parset, "");

  BOOST_TEST(predict.GetBufferSize() == 0u);
}

BOOST_AUTO_TEST_CASE(no_regions) {
  dp3::common::ParameterSet parset;
  parset.add("regions", "im-updside-down.reg");
  parset.add("images", "sources-model.fits");

  BOOST_CHECK_THROW(std::make_unique<WGridderPredict>(parset, ""),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(no_models) {
  dp3::common::ParameterSet parset;
  parset.add("regions", "sources.reg");

  BOOST_CHECK_THROW(std::make_unique<WGridderPredict>(parset, ""),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(update_info_wrong) {
  const dp3::common::ParameterSet parset = CreateParset();

  dp3::base::DPInfo info;

  WGridderPredict predict(parset, "");

  BOOST_CHECK_THROW(predict.updateInfo(info), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(update_info) {
  const dp3::common::ParameterSet parset = CreateParset();

  const dp3::base::DPInfo info = InitInfo();

  WGridderPredict predict(parset, "");

  predict.updateInfo(info);

  BOOST_TEST(predict.GetBufferSize() > 0u);
}

BOOST_AUTO_TEST_CASE(process) {
  // This test compares the visibilities generated bu WGridderPredict to those
  // generated by OnePredict. WGridderPredict uses the model image
  // "sources-model.fits" as input. OnePredict uses the catalog
  // "sources.skymodel" as input. The image and the catalog contain the same
  // sources, so the visibilities should be the same.

  const dp3::common::ParameterSet parset = CreateParset();

  dp3::common::ParameterSet one_predict_parset = parset;
  one_predict_parset.add("sourcedb", "sources.skymodel");
  OnePredict one_predict(one_predict_parset, "", std::vector<std::string>());
  WGridderPredict wgridder_predict(parset, "");

  wgridder_predict.SetBufferSize(kNTimeSteps);

  auto one_predict_result_step = std::make_shared<MultiResultStep>(kNTimeSteps);
  auto wgridder_predict_result_step =
      std::make_shared<MultiResultStep>(kNTimeSteps);

  one_predict.setNextStep(one_predict_result_step);
  wgridder_predict.setNextStep(wgridder_predict_result_step);

  aocommon::FitsReader fits_reader("sources-model.fits");
  casacore::MDirection phase_center(
      casacore::Quantity(fits_reader.PhaseCentreRA(), "rad"),
      casacore::Quantity(fits_reader.PhaseCentreDec(), "rad"));
  dp3::base::DPInfo info = InitInfo();
  info.setPhaseCenter(phase_center);
  wgridder_predict.setInfo(info);
  one_predict.setInfo(info);

  for (std::size_t i = 0; i < kNTimeSteps; ++i) {
    std::unique_ptr<DPBuffer> buffer =
        CreateBuffer(kFirstTime + i * kInterval, kInterval);
    std::unique_ptr<DPBuffer> one_predict_buffer =
        std::make_unique<DPBuffer>(*buffer);
    wgridder_predict.process(std::move(buffer));
    one_predict.process(std::move(one_predict_buffer));
  }
  wgridder_predict.finish();
  one_predict.finish();

  BOOST_TEST(wgridder_predict_result_step->size() == kNTimeSteps);
  BOOST_TEST(one_predict_result_step->size() == kNTimeSteps);
  for (std::size_t i = 0; i < kNTimeSteps; ++i) {
    BOOST_CHECK(xt::allclose(wgridder_predict_result_step->get()[i]->GetData(),
                             one_predict_result_step->get()[i]->GetData(),
                             kVisibilitiesRelativeTolerance,
                             kVisibilitiesAbsoluteTolerance));
  }
}

BOOST_AUTO_TEST_SUITE_END()
