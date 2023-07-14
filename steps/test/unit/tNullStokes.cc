// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cstddef>

#include <boost/test/unit_test.hpp>
#include <xtensor/xcomplex.hpp>

#include "../../NullStokes.h"
#include "../../ResultStep.h"
#include "tStepCommon.h"
#include "mock/ThrowStep.h"
#include "../../../common/ParameterSet.h"

BOOST_AUTO_TEST_SUITE(null_stokes)

namespace {
// Constants with arbitrary values used solely to ensure unexpected data isn't
// changed by the step
constexpr double kTime{4.87128e+09};
constexpr double kExposure{10.0139};
constexpr double kWeightStepSize{0.7};
constexpr double kUvwStepSize{0.1};
}  // namespace

bool TestNullStokes(std::size_t n_baselines, std::size_t n_channels,
                    std::size_t n_correlations, bool modify_q, bool modify_u) {
  dp3::common::ParameterSet parset;
  if (modify_q) {
    parset.add("modify_q", "true");
  }
  if (modify_u) {
    parset.add("modify_u", "true");
  }

  auto null_stokes_step = std::make_shared<dp3::steps::NullStokes>(parset, "");
  auto result_step = std::make_shared<dp3::steps::ResultStep>();
  null_stokes_step->setNextStep(result_step);

  auto buffer = std::make_unique<dp3::base::DPBuffer>(kTime, kExposure);
  const std::array<std::size_t, 3> data_shape = {n_baselines, n_channels,
                                                 n_correlations};
  const std::array<std::size_t, 2> uvw_shape = {n_baselines, std::size_t{3}};
  // Fill the data with arbitrary values in which we can detect if the correct
  // changes are applied
  buffer->GetData().resize(data_shape);
  const std::size_t data_size = buffer->GetData().size();
  for (std::size_t i = 0; i < data_size; ++i) {
    buffer->GetData().data()[i] = std::complex<float>(
        i + 1 * 10, static_cast<int>(i) - 1000 + static_cast<int>(1) * 6);
  }
  // Fill the flags, weights and uvw with arbitrary values so that we can detect
  // any unwanted changes
  xt::xtensor<std::complex<float>, 3> input_data = buffer->GetData();
  buffer->ResizeFlags(data_shape);
  buffer->GetFlags().fill(false);
  xt::xtensor<bool, 3> expected_flags = buffer->GetFlags();
  buffer->ResizeWeights(data_shape);
  buffer->GetWeights() =
      xt::arange(0.0,
                 kWeightStepSize * n_baselines * n_channels * n_correlations,
                 kWeightStepSize)
          .reshape(data_shape);
  xt::xtensor<float, 3> expected_weights = buffer->GetWeights();
  buffer->ResizeUvw(n_baselines);
  buffer->GetUvw() =
      xt::arange(0.0, kUvwStepSize * n_baselines * n_channels * n_correlations,
                 kUvwStepSize)
          .reshape(uvw_shape);
  xt::xtensor<double, 2> expected_uvw = buffer->GetUvw();

  null_stokes_step->process(std::move(buffer));

  // Ensure nothing other than the data has changed
  BOOST_CHECK(result_step->get().getTime() == kTime);
  BOOST_CHECK(result_step->get().getExposure() == kExposure);
  BOOST_TEST(result_step->get().GetWeights() == expected_weights);
  BOOST_TEST(result_step->get().GetFlags() == expected_flags);
  BOOST_TEST(result_step->get().GetUvw() == expected_uvw);

  // Test that the data changes are correct
  if (!(modify_q || modify_u)) {
    BOOST_TEST(result_step->get().GetData() == input_data);
  } else {
    BOOST_TEST(result_step->get().GetData() != input_data);
    const std::complex<float>* data = result_step->get().GetData().data();
    for (std::size_t i = 0; i < data_size; i += 4) {
      std::complex<float> xx = data[i];
      std::complex<float> xy = data[i + 1];
      std::complex<float> yx = data[i + 2];
      std::complex<float> yy = data[i + 3];
      if (modify_q) {
        BOOST_CHECK(yy == xx);
      }
      if (modify_u) {
        BOOST_CHECK(yx == -xy);
      }
    }
  }
  return true;
}

// Perform a very small/minimal run of the ArrayBeam step and compare against
// expected result to help catch any regressions
BOOST_AUTO_TEST_CASE(test_basic_nullstokes) {
  BOOST_CHECK(TestNullStokes(28, 8, 4, false, false));
  BOOST_CHECK(TestNullStokes(28, 8, 4, true, true));
  BOOST_CHECK(TestNullStokes(28, 8, 4, true, false));
  BOOST_CHECK(TestNullStokes(28, 8, 4, false, true));
}

BOOST_AUTO_TEST_SUITE_END()
