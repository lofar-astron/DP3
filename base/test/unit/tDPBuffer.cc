// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../DPBuffer.h"

using dp3::base::DPBuffer;

namespace {
const double kTime = 42.0;
const double kExposure = 4.0;
const dp3::common::rownr_t kRowNr = 42;
const std::complex<float> kDataValue{42.0, -42.0};
const float kWeightValue = 0.5;
const double KUVWValue = 2.0;
const size_t kRowNrs = 10;
const size_t kNCorrelations = 4;
const size_t kNChannels = 3;
const size_t kNBaselines = 5;
const size_t kNTimeAvg = 6;
const casacore::IPosition kDimensions(3, kNCorrelations, kNChannels,
                                      kNBaselines);
const casacore::IPosition kFRFDimensions(3, kNChannels, kNTimeAvg, kNBaselines);

template <class T>
void CompareArray(const casacore::Array<T>& left,
                  const casacore::Array<T>& right) {
  BOOST_REQUIRE(left.shape() == right.shape());
  BOOST_CHECK_EQUAL_COLLECTIONS(left.begin(), left.end(), right.begin(),
                                right.end());
}

/// Verify that two DPBuffers do not share the same data, e.g., using
/// casacore references.
void CheckIndependent(const DPBuffer& left, const DPBuffer& right) {
  BOOST_CHECK_NE(left.getRowNrs().data(), right.getRowNrs().data());
  BOOST_CHECK_NE(left.getData().data(), right.getData().data());
  BOOST_CHECK_NE(left.getFlags().data(), right.getFlags().data());
  BOOST_CHECK_NE(left.getUVW().data(), right.getUVW().data());
  BOOST_CHECK_NE(left.getWeights().data(), right.getWeights().data());
  BOOST_CHECK_NE(left.getFullResFlags().data(), right.getFullResFlags().data());
}

DPBuffer CreateFilledBuffer() {
  DPBuffer buffer;
  buffer.setTime(kTime);
  buffer.setExposure(kExposure);
  buffer.setRowNrs(casacore::Vector<dp3::common::rownr_t>(kRowNrs, kRowNr));
  buffer.setData(casacore::Cube<std::complex<float>>(kDimensions, kDataValue));
  buffer.setFlags(casacore::Cube<bool>(kDimensions, false));
  buffer.setUVW(casacore::Matrix<double>(3, kNBaselines, KUVWValue));
  buffer.setWeights(casacore::Cube<float>(kDimensions, kWeightValue));
  buffer.setFullResFlags(casacore::Cube<bool>(kFRFDimensions, true));
  return buffer;
}

void CheckFilledBuffer(const DPBuffer& buffer) {
  BOOST_CHECK(buffer.getTime() == kTime);
  BOOST_CHECK(buffer.getExposure() == kExposure);
  CompareArray(buffer.getRowNrs(),
               casacore::Vector<dp3::common::rownr_t>(kRowNrs, kRowNr));
  CompareArray(buffer.getData(),
               casacore::Cube<std::complex<float>>(kDimensions, kDataValue));
  CompareArray(buffer.getFlags(), casacore::Cube<bool>(kDimensions, false));
  CompareArray(buffer.getUVW(),
               casacore::Matrix<double>(3, kNBaselines, KUVWValue));
  CompareArray(buffer.getWeights(),
               casacore::Cube<float>(kDimensions, kWeightValue));
  CompareArray(buffer.getFullResFlags(),
               casacore::Cube<bool>(kFRFDimensions, true));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(dpbuffer)

BOOST_AUTO_TEST_CASE(constructor) {
  const DPBuffer buffer;
  BOOST_CHECK(buffer.getTime() == 0.0);
  BOOST_CHECK(buffer.getExposure() == 0.0);
  BOOST_CHECK(buffer.getRowNrs().empty());
  BOOST_CHECK(buffer.getData().empty());
  BOOST_CHECK(buffer.getFlags().empty());
  BOOST_CHECK(buffer.getUVW().empty());
  BOOST_CHECK(buffer.getWeights().empty());
  BOOST_CHECK(buffer.getFullResFlags().empty());
}

BOOST_AUTO_TEST_CASE(move_constructor) {
  DPBuffer source = CreateFilledBuffer();
  CheckFilledBuffer(source);

  DPBuffer move_constructed(std::move(source));
  CheckFilledBuffer(move_constructed);
  CheckIndependent(source, move_constructed);
}

BOOST_AUTO_TEST_CASE(move_assignment) {
  DPBuffer source = CreateFilledBuffer();
  DPBuffer move_assigned;
  CheckFilledBuffer(source);

  move_assigned = std::move(source);
  CheckFilledBuffer(move_assigned);
  CheckIndependent(source, move_assigned);
}

BOOST_AUTO_TEST_SUITE_END()