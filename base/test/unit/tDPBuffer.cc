// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include <xtensor/xio.hpp>

#include <dp3/base/DPBuffer.h>

using dp3::base::DPBuffer;

namespace {
const double kTime = 42.0;
const double kExposure = 4.0;
const dp3::common::rownr_t kRowNr = 42;
const std::complex<float> kDataValue{42.0, -42.0};
const std::string kFooDataName = "foo";
const std::string kBarDataName = "bar";
const std::complex<float> kFooDataValue{43.0, -43.0};
const std::complex<float> kBarDataValue{44.0, -44.0};
const float kWeightValue = 0.5;
const double kUVWValue = 2.0;
const size_t kRowNrs = 10;
const size_t kNCorrelations = 4;
const size_t kNChannels = 3;
const size_t kNBaselines = 5;
const size_t kNTimeAvg = 6;
const std::array<std::size_t, 3> kShape{kNBaselines, kNChannels,
                                        kNCorrelations};

template <class T>
void CompareArray(const casacore::Array<T>& left,
                  const casacore::Array<T>& right) {
  BOOST_REQUIRE(left.shape() == right.shape());
  BOOST_CHECK_EQUAL_COLLECTIONS(left.begin(), left.end(), right.begin(),
                                right.end());
}

/// Verify that two DPBuffers do share the same data, e.g., using
/// casacore references.
void CheckDependent(const DPBuffer& left, const DPBuffer& right) {
  BOOST_CHECK_EQUAL(left.getRowNrs().data(), right.getRowNrs().data());
  BOOST_CHECK_EQUAL(left.GetData().data(), right.GetData().data());
  BOOST_CHECK_EQUAL(left.GetFlags().data(), right.GetFlags().data());
  BOOST_CHECK_EQUAL(left.GetUvw().data(), right.GetUvw().data());
  BOOST_CHECK_EQUAL(left.GetWeights().data(), right.GetWeights().data());
}

/// Verify that two DPBuffers do not share the same data, e.g., using
/// casacore references.
void CheckIndependent(const DPBuffer& left, const DPBuffer& right) {
  BOOST_CHECK_NE(left.getRowNrs().data(), right.getRowNrs().data());
  BOOST_CHECK_NE(left.GetData().data(), right.GetData().data());
  BOOST_CHECK_NE(left.GetFlags().data(), right.GetFlags().data());
  BOOST_CHECK_NE(left.GetUvw().data(), right.GetUvw().data());
  BOOST_CHECK_NE(left.GetWeights().data(), right.GetWeights().data());
}

DPBuffer CreateFilledBuffer() {
  DPBuffer buffer;
  buffer.setTime(kTime);
  buffer.setExposure(kExposure);
  buffer.setRowNrs(casacore::Vector<dp3::common::rownr_t>(kRowNrs, kRowNr));
  buffer.ResizeData(kNBaselines, kNChannels, kNCorrelations);
  buffer.AddData(kFooDataName);
  buffer.AddData(kBarDataName);
  buffer.ResizeFlags(kNBaselines, kNChannels, kNCorrelations);
  buffer.ResizeWeights(kNBaselines, kNChannels, kNCorrelations);
  buffer.ResizeUvw(kNBaselines);
  buffer.GetData().fill(kDataValue);
  buffer.GetData(kFooDataName).fill(kFooDataValue);
  buffer.GetData(kBarDataName).fill(kBarDataValue);
  buffer.GetFlags().fill(false);
  buffer.GetUvw().fill(kUVWValue);
  buffer.GetWeights().fill(kWeightValue);
  return buffer;
}

void CheckFilledBuffer(const DPBuffer& buffer) {
  BOOST_CHECK(buffer.getTime() == kTime);
  BOOST_CHECK(buffer.getExposure() == kExposure);
  CompareArray(buffer.getRowNrs(),
               casacore::Vector<dp3::common::rownr_t>(kRowNrs, kRowNr));

  const xt::xtensor<std::complex<float>, 3> data(kShape, kDataValue);
  const xt::xtensor<std::complex<float>, 3> foo_data(kShape, kFooDataValue);
  const xt::xtensor<std::complex<float>, 3> bar_data(kShape, kBarDataValue);
  const xt::xtensor<bool, 3> flags(kShape, false);
  const xt::xtensor<double, 2> uvw({kNBaselines, 3}, kUVWValue);
  const xt::xtensor<float, 3> weights(kShape, kWeightValue);
  BOOST_CHECK_EQUAL(buffer.GetData(), data);
  BOOST_CHECK_EQUAL(buffer.GetData(kFooDataName), foo_data);
  BOOST_CHECK_EQUAL(buffer.GetData(kBarDataName), bar_data);
  BOOST_CHECK_EQUAL(buffer.GetFlags(), flags);
  BOOST_CHECK_EQUAL(buffer.GetUvw(), uvw);
  BOOST_CHECK_EQUAL(buffer.GetWeights(), weights);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(dpbuffer)

BOOST_AUTO_TEST_CASE(constructor) {
  const DPBuffer buffer;
  BOOST_CHECK(buffer.getTime() == 0.0);
  BOOST_CHECK(buffer.getExposure() == 0.0);
  BOOST_CHECK(buffer.getRowNrs().empty());
  BOOST_CHECK(buffer.GetData().size() == 0);
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(!buffer.HasData(kFooDataName));
  BOOST_CHECK(buffer.GetFlags().size() == 0);
  BOOST_CHECK(buffer.GetUvw().size() == 0);
  BOOST_CHECK(buffer.GetWeights().size() == 0);
}

BOOST_AUTO_TEST_CASE(copy_constructor) {
  const DPBuffer source = CreateFilledBuffer();
  CheckFilledBuffer(source);

  const DPBuffer copy_constructed(source);
  CheckFilledBuffer(copy_constructed);
  CheckDependent(source, copy_constructed);
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

BOOST_AUTO_TEST_CASE(copy_to_empty) {
  const DPBuffer source = CreateFilledBuffer();
  DPBuffer copy;
  copy.copy(source);
  CheckFilledBuffer(copy);
  CheckIndependent(source, copy);
}

BOOST_AUTO_TEST_CASE(copy_to_reference_copy) {
  const DPBuffer source = CreateFilledBuffer();
  DPBuffer copy(source);  // 'copy' gets references (see copy_constructor test).
  copy.copy(source);      // 'copy' becomes indedependent.
  CheckFilledBuffer(copy);
  CheckIndependent(source, copy);
}

BOOST_AUTO_TEST_CASE(make_independent) {
  using dp3::common::Fields;
  const DPBuffer source = CreateFilledBuffer();
  DPBuffer copy(source);  // 'copy' gets references (see copy_constructor test).
  CheckDependent(source, copy);
  copy.getRowNrs().unique();  // MakeIndependent does not support row numbers.
  copy.MakeIndependent(
      Fields(Fields::Single::kData) | Fields(Fields::Single::kFlags) |
      Fields(Fields::Single::kWeights) | Fields(Fields::Single::kUvw));
  CheckFilledBuffer(copy);
  CheckIndependent(source, copy);
}

BOOST_AUTO_TEST_CASE(resize_data) {
  const size_t kNewNBaselines = kNBaselines + 1;
  const size_t kNewNChannels = kNChannels + 1;
  const std::array<size_t, 3> kNewShape{kNewNBaselines, kNewNChannels,
                                        kNCorrelations};
  DPBuffer buffer = CreateFilledBuffer();
  buffer.ResizeData(kNewNBaselines, kNewNChannels, kNCorrelations);
  BOOST_CHECK(buffer.GetData().shape() == kNewShape);
  BOOST_CHECK(buffer.GetData(kFooDataName).shape() == kNewShape);
  BOOST_CHECK(buffer.GetData(kBarDataName).shape() == kNewShape);
}

BOOST_AUTO_TEST_CASE(remove_data_one_by_one) {
  DPBuffer buffer = CreateFilledBuffer();
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(buffer.HasData(kFooDataName));
  BOOST_CHECK(buffer.HasData(kBarDataName));

  buffer.RemoveData(kFooDataName);
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(!buffer.HasData(kFooDataName));
  BOOST_CHECK(buffer.HasData(kBarDataName));

  buffer.RemoveData(kBarDataName);
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(!buffer.HasData(kFooDataName));
  BOOST_CHECK(!buffer.HasData(kBarDataName));
}

BOOST_AUTO_TEST_CASE(remove_data_all_at_once) {
  DPBuffer buffer = CreateFilledBuffer();
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(buffer.HasData(kFooDataName));
  BOOST_CHECK(buffer.HasData(kBarDataName));

  buffer.RemoveData();
  BOOST_CHECK(buffer.HasData());
  BOOST_CHECK(!buffer.HasData(kFooDataName));
  BOOST_CHECK(!buffer.HasData(kBarDataName));
}

BOOST_AUTO_TEST_SUITE_END()
