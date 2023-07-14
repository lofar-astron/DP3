// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/SolverTools.h"

#include <limits>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <xtensor/xcomplex.hpp>

using dp3::base::DPBuffer;
using dp3::ddecal::AssignAndWeight;

namespace {

// Values in the first unweighted buffer.
// If a test uses multiple buffers, successive buffers have the buffer index
// added to the real component of the weighted (model) visibility data.
const std::complex<float> kBaseDataValue{1, 1};
const std::complex<float> kBaseModelDataValue{42, 42};
// Make AssignAndWeight double all values.
const float kWeightFactor = 2.0f;
const float kWeightValue = kWeightFactor * kWeightFactor;

const std::size_t kNBaselines = 1;
const std::size_t kNChannels = 2;
const std::size_t kNCorrelations = 4;
const std::array<std::size_t, 3> kShape{kNBaselines, kNChannels,
                                        kNCorrelations};
const std::vector<std::string> kDirectionNames{"foodir", "bardir"};

template <std::size_t NBuffers>
class AssignAndWeightFixture {
 public:
  AssignAndWeightFixture() {
    for (std::size_t i = 0; i < NBuffers; ++i) {
      unweighted_buffers_.emplace_back(std::make_unique<DPBuffer>());
      weighted_buffers_.emplace_back();

      DPBuffer& buffer = *unweighted_buffers_.back();
      buffer.GetData().resize(kShape);
      buffer.GetData().fill(kBaseDataValue + static_cast<float>(i));
      buffer.ResizeWeights(kShape);
      buffer.GetWeights().fill(kWeightValue);
      buffer.ResizeFlags(kShape);
      buffer.GetFlags().fill(false);
      for (const std::string& name : kDirectionNames) {
        buffer.AddData(name);
        buffer.GetData(name).fill(kBaseModelDataValue + static_cast<float>(i));
      }
    }
  }

  /// @return The expected value for the main data in a buffer.
  static std::complex<float> WeightedValue(std::size_t buffer_index) {
    return (kBaseDataValue + static_cast<float>(buffer_index)) * kWeightFactor;
  }

  /// @return The expected value for the model data in a buffer.
  static std::complex<float> WeightedModelValue(std::size_t buffer_index) {
    return (kBaseModelDataValue + static_cast<float>(buffer_index)) *
           kWeightFactor;
  }

  std::vector<std::unique_ptr<DPBuffer>>& UnweightedBuffers() {
    return unweighted_buffers_;
  }
  const std::vector<DPBuffer>& WeightedBuffers() const {
    return weighted_buffers_;
  }
  const DPBuffer::DataType& WeightedData(std::size_t buffer_index,
                                         const std::string& name) const {
    return weighted_buffers_[buffer_index].GetData(name);
  }

  void DoAssignAndWeight(const bool keep_original_model_data) {
    dp3::ddecal::AssignAndWeight(unweighted_buffers_, kDirectionNames,
                                 weighted_buffers_, keep_original_model_data);
  }

 private:
  std::vector<std::unique_ptr<DPBuffer>> unweighted_buffers_;
  std::vector<DPBuffer> weighted_buffers_;
};

template <typename T>
void CheckTensor(const T& tensor,
                 const typename T::value_type& expected_value) {
  BOOST_TEST(tensor.shape() == kShape);
  BOOST_TEST(xt::allclose(tensor, expected_value));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(solvertools)

BOOST_DATA_TEST_CASE_F(AssignAndWeightFixture<1>,
                       assign_and_weight_move_or_keep_data,
                       boost::unit_test::data::make({false, true}),
                       keep_original_model_data) {
  DoAssignAndWeight(keep_original_model_data);

  const DPBuffer& unweighted_buffer = *UnweightedBuffers().front();
  const DPBuffer& weighted_buffer = WeightedBuffers().front();

  // Check that the main data and weights did not change.
  CheckTensor(unweighted_buffer.GetData(), kBaseDataValue);
  CheckTensor(unweighted_buffer.GetWeights(), kWeightValue);

  // Check that the main weighted data is correct.
  CheckTensor(weighted_buffer.GetData(), WeightedValue(0));

  for (const std::string& name : kDirectionNames) {
    BOOST_TEST(unweighted_buffer.HasData(name) == keep_original_model_data);
    BOOST_TEST(weighted_buffer.HasData(name));

    if (keep_original_model_data) {
      // Check that the model data did not change.
      CheckTensor(unweighted_buffer.GetData(name), kBaseModelDataValue);
    }

    // Check that the weighted model data is correct.
    CheckTensor(weighted_buffer.GetData(name), WeightedModelValue(0));
  }
}

BOOST_FIXTURE_TEST_CASE(assign_and_weight_nan, AssignAndWeightFixture<3>) {
  // This test case passes three DPBuffers to AssignAndWeight:
  // buffers[0]: A DPBuffer that has a single NaN value in the main data buffer,
  //             in the first channel.
  // buffers[1]: A DPBuffer that has a single NaN value in a model data buffer,
  //             in the second/last channel.
  // buffers[2]: A DPBuffer that has NaN values in all data buffers, for
  //             all channels and correllations.
  // (assign_and_weight_move_or_keep_data tests a DPBuffer without NaN's.)

  const float kNan = std::numeric_limits<float>::quiet_NaN();
  const std::complex<float> kNanReal{kNan, 0.0f};
  const std::complex<float> kNanImag{0.0f, kNan};
  const std::complex<float> kNanRealImag{kNan, kNan};

  UnweightedBuffers()[0]->GetData()(0, 0, kNCorrelations - 1) = kNanReal;
  UnweightedBuffers()[1]->GetData(kDirectionNames.front())(0, kNChannels - 1,
                                                           0) = kNanImag;
  UnweightedBuffers()[2]->GetData().fill(kNanRealImag);
  for (const std::string& name : kDirectionNames) {
    UnweightedBuffers()[2]->GetData(name).fill(kNanRealImag);
  }

  DoAssignAndWeight(false);

  // Check all correlations have the same result: Either zero or the expected
  // weighted value.
  const std::complex<float> kZero{0.0f, 0.0f};
  for (std::size_t correlation = 0; correlation < kNCorrelations;
       ++correlation) {
    // Check main data buffers.
    BOOST_TEST(WeightedData(0, "")(0, 0, correlation) == kZero);
    BOOST_TEST(WeightedData(0, "")(0, 1, correlation) == WeightedValue(0));

    BOOST_TEST(WeightedData(1, "")(0, 0, correlation) == WeightedValue(1));
    BOOST_TEST(WeightedData(1, "")(0, 1, correlation) == kZero);

    BOOST_TEST(WeightedData(2, "")(0, 0, correlation) == kZero);
    BOOST_TEST(WeightedData(2, "")(0, 1, correlation) == kZero);

    // Check model data buffers.
    for (const std::string& name : kDirectionNames) {
      BOOST_TEST(WeightedData(0, name)(0, 0, correlation) == kZero);
      BOOST_TEST(WeightedData(0, name)(0, 1, correlation) ==
                 WeightedModelValue(0));

      BOOST_TEST(WeightedData(1, name)(0, 0, correlation) ==
                 WeightedModelValue(1));
      BOOST_TEST(WeightedData(1, name)(0, 1, correlation) == kZero);

      BOOST_TEST(WeightedData(2, name)(0, 0, correlation) == kZero);
      BOOST_TEST(WeightedData(2, name)(0, 1, correlation) == kZero);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(assign_and_weight_flags, AssignAndWeightFixture<3>) {
  // This test case passes three DPBuffers to AssignAndWeight:
  // buffers[0]: A DPBuffer where all flags are set.
  // buffers[1]: A DPBuffer where the first flag is set.
  // buffers[2]: A DPBuffer where the last flag is set.
  // (assign_and_weight_move_or_keep_data tests a DPBuffer where all flags are
  // unset.)

  UnweightedBuffers()[0]->GetFlags().fill(true);
  UnweightedBuffers()[1]->GetFlags()(0, 0, 0) = true;
  UnweightedBuffers()[2]->GetFlags()(kNBaselines - 1, kNChannels - 1,
                                     kNCorrelations - 1) = true;

  DoAssignAndWeight(false);

  // Check all correlations have the same result: Either zero or the expected
  // weighted value.
  const std::complex<float> kZero{0.0f, 0.0f};
  for (std::size_t correlation = 0; correlation < kNCorrelations;
       ++correlation) {
    // Check main data buffers.
    BOOST_TEST(WeightedData(0, "")(0, 0, correlation) == kZero);
    BOOST_TEST(WeightedData(0, "")(0, 1, correlation) == kZero);

    BOOST_TEST(WeightedData(1, "")(0, 0, correlation) == kZero);
    BOOST_TEST(WeightedData(1, "")(0, 1, correlation) == WeightedValue(1));

    BOOST_TEST(WeightedData(2, "")(0, 0, correlation) == WeightedValue(2));
    BOOST_TEST(WeightedData(2, "")(0, 1, correlation) == kZero);

    // Check model data buffers.
    for (const std::string& name : kDirectionNames) {
      BOOST_TEST(WeightedData(0, name)(0, 0, correlation) == kZero);
      BOOST_TEST(WeightedData(0, name)(0, 1, correlation) == kZero);

      BOOST_TEST(WeightedData(1, name)(0, 0, correlation) == kZero);
      BOOST_TEST(WeightedData(1, name)(0, 1, correlation) ==
                 WeightedModelValue(1));

      BOOST_TEST(WeightedData(2, name)(0, 0, correlation) ==
                 WeightedModelValue(2));
      BOOST_TEST(WeightedData(2, name)(0, 1, correlation) == kZero);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
