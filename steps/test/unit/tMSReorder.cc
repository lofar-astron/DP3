// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MSReorder.h"
#include "aocommon/polarization.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <tuple>
#include <vector>
#include <complex>
#include <cstddef>

using namespace dp3::reorder;
using namespace aocommon;

typedef std::vector<std::complex<float>> ComplexVector;

const ComplexVector test_data{10, 11, 12, 13, 20, 21, 22, 23,
                              30, 31, 32, 33, 40, 41, 42, 43};
const std::vector<float> test_weights{0.1f, 0.1f, 0.1f, 0.1f, 0.2f, 0.2f,
                                      0.2f, 0.2f, 0.3f, 0.3f, 0.3f, 0.3f,
                                      0.4f, 0.4f, 0.4f, 0.4f};
// std::vector<bool> does have a .data() method
const bool test_flags[] = {false, false, false, false, false, false,
                           false, false, true,  true,  true,  true,
                           false, false, false, false};

BOOST_AUTO_TEST_SUITE(msreorder)

BOOST_AUTO_TEST_CASE(filename_prefix) {
  const std::string filename_prefix = GetFilenamePrefix("test.ms", "");
  BOOST_CHECK_EQUAL(filename_prefix, "test.ms");
}

BOOST_AUTO_TEST_CASE(filename_prefix_tmp_dir) {
  const std::string filename_prefix = GetFilenamePrefix("tmp/test.ms", "tmp");
  BOOST_CHECK_EQUAL(filename_prefix, "tmp/test.ms");
}

BOOST_AUTO_TEST_CASE(filenameprefix_remove_trailing_separator) {
  const std::string filename_prefix = GetFilenamePrefix("test.ms/", "");
  BOOST_CHECK_EQUAL(filename_prefix, "test.ms");
}

BOOST_AUTO_TEST_CASE(metafilename) {
  const std::string meta_filename = GetMetaFilename("test.ms", "", 0);
  BOOST_CHECK_EQUAL(meta_filename, "test.ms-spw0-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(metafilename_with_tmp) {
  const std::string meta_filename = GetMetaFilename("test.ms", "tmp", 0);
  BOOST_CHECK_EQUAL(meta_filename, "tmp/test.ms-spw0-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(metafilename_ddi) {
  const std::string meta_filename = GetMetaFilename("test.ms", "", 1);
  BOOST_CHECK_EQUAL(meta_filename, "test.ms-spw1-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(partprefix) {
  const std::string partprefix =
      GetPartPrefix("test.ms", 0, Polarization::StokesI, 0, "");
  BOOST_CHECK_EQUAL(partprefix, "test.ms-part0000-I-b0");
}

BOOST_DATA_TEST_CASE(extractdata_linear_pol_to_stokes,
                     boost::unit_test::data::make({Polarization::StokesI,
                                                   Polarization::StokesQ,
                                                   Polarization::StokesU,
                                                   Polarization::StokesV}),
                     pol_out) {
  size_t nchan = 4;
  std::vector<std::complex<float>> actual(nchan);
  // Boost auto test case doesn't suport STL inside data::make
  std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::StokesI, ComplexVector{11.5f, 21.5f, 31.5f, 41.5f}},
      {Polarization::StokesQ, ComplexVector{-1.5f, -1.5f, -1.5f, -1.5f}},
      {Polarization::StokesU, ComplexVector{11.5f, 21.5f, 31.5f, 41.5f}},
      {Polarization::StokesV, ComplexVector{std::complex<float>(0.0f, 0.5f),
                                            std::complex<float>(0.0f, 0.5f),
                                            std::complex<float>(0.0f, 0.5f),
                                            std::complex<float>(0.0f, 0.5f)}},
  };

  ExtractData(actual.data(), 0, nchan,
              std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                         Polarization::YX, Polarization::YY},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected[pol_out].begin(),
                                expected[pol_out].end());
}

BOOST_DATA_TEST_CASE(extractdata_circular_pol_to_stokes,
                     boost::unit_test::data::make({Polarization::StokesI,
                                                   Polarization::StokesQ,
                                                   Polarization::StokesU,
                                                   Polarization::StokesV}),
                     pol_out) {
  size_t nchan = 4;
  std::vector<std::complex<float>> actual(nchan);
  // Boost auto test case doesn't suport STL inside data::make
  std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::StokesI, ComplexVector{11.5f, 21.5f, 31.5f, 41.5f}},
      {Polarization::StokesQ, ComplexVector{11.5f, 21.5f, 31.5f, 41.5f}},
      {Polarization::StokesU, ComplexVector{std::complex<float>(0.0f, 0.5f),
                                            std::complex<float>(0.0f, 0.5f),
                                            std::complex<float>(0.0f, 0.5f),
                                            std::complex<float>(0.0f, 0.5f)}},
      {Polarization::StokesV, ComplexVector{-1.5f, -1.5f, -1.5f, -1.5f}},
  };

  ExtractData(actual.data(), 0, nchan,
              std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                         Polarization::LR, Polarization::LL},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected[pol_out].begin(),
                                expected[pol_out].end());
}

BOOST_DATA_TEST_CASE(
    extractdata_linear_pol_to_linear,
    boost::unit_test::data::make({Polarization::XX, Polarization::XY,
                                  Polarization::YX, Polarization::YY}),
    pol_out) {
  size_t nchan = 4;
  std::vector<std::complex<float>> actual(nchan);
  // Boost auto test case doesn't suport STL inside data::make
  std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::XX, ComplexVector{10, 20, 30, 40}},
      {Polarization::XY, ComplexVector{11, 21, 31, 41}},
      {Polarization::YX, ComplexVector{12, 22, 32, 42}},
      {Polarization::YY, ComplexVector{13, 23, 33, 43}},
  };

  ExtractData(actual.data(), 0, nchan,
              std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                         Polarization::YX, Polarization::YY},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected[pol_out].begin(),
                                expected[pol_out].end());
}

BOOST_DATA_TEST_CASE(
    extractdata_circular_pol_to_circular,
    boost::unit_test::data::make({Polarization::RR, Polarization::RL,
                                  Polarization::LR, Polarization::LL}),
    pol_out) {
  size_t nchan = 4;
  std::vector<std::complex<float>> actual(nchan);
  // Boost auto test case doesn't suport STL inside data::make
  std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::RR, ComplexVector{10, 20, 30, 40}},
      {Polarization::RL, ComplexVector{11, 21, 31, 41}},
      {Polarization::LR, ComplexVector{12, 22, 32, 42}},
      {Polarization::LL, ComplexVector{13, 23, 33, 43}},
  };

  ExtractData(actual.data(), 0, nchan,
              std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                         Polarization::LR, Polarization::LL},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected[pol_out].begin(),
                                expected[pol_out].end());
}

BOOST_AUTO_TEST_CASE(extractdata_linear_pol_to_instrumental) {
  size_t nchan = 4;
  size_t pol_per_file = 4;
  std::vector<std::complex<float>> actual(nchan * pol_per_file);
  ComplexVector expected(test_data);
  PolarizationEnum pol_out = Polarization::Instrumental;

  ExtractData(actual.data(), 0, nchan,
              std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                         Polarization::YX, Polarization::YY},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extractdata_circular_pol_to_instrumental) {
  size_t nchan = 4;
  size_t pol_per_file = 4;
  std::vector<std::complex<float>> actual(nchan * pol_per_file);
  ComplexVector expected(test_data);
  PolarizationEnum pol_out = Polarization::Instrumental;

  ExtractData(actual.data(), 0, nchan,
              std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                         Polarization::LR, Polarization::LL},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extractdata_linear_pol_to_diag_instrumental) {
  size_t nchan = 4;
  size_t pol_per_file = 2;
  std::vector<std::complex<float>> actual(nchan * pol_per_file);
  ComplexVector expected{10, 13, 20, 23, 30, 33, 40, 43};
  PolarizationEnum pol_out = Polarization::DiagonalInstrumental;

  ExtractData(actual.data(), 0, nchan,
              std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                         Polarization::YX, Polarization::YY},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extractdata_circular_pol_to_diag_instrumental) {
  size_t nchan = 4;
  size_t pol_per_file = 2;
  std::vector<std::complex<float>> actual(nchan * pol_per_file);
  ComplexVector expected{10, 13, 20, 23, 30, 33, 40, 43};
  PolarizationEnum pol_out = Polarization::DiagonalInstrumental;

  ExtractData(actual.data(), 0, nchan,
              std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                         Polarization::LR, Polarization::LL},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_DATA_TEST_CASE(extractweights_linear_pol_to_stokes,
                     boost::unit_test::data::make({Polarization::StokesI,
                                                   Polarization::StokesQ,
                                                   Polarization::StokesU,
                                                   Polarization::StokesV}),
                     pol_out) {
  size_t nchan = 4;
  std::vector<float> actual(nchan);
  // Boost auto test case doesn't suport STL inside data::make
  std::map<PolarizationEnum, std::vector<float>> expected{
      {Polarization::StokesI, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesQ, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesU, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesV, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
  };

  ExtractWeights(actual.data(), 0, nchan,
                 std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected[pol_out].begin(),
                                expected[pol_out].end());
}

BOOST_DATA_TEST_CASE(extractweights_circular_pol_to_stokes,
                     boost::unit_test::data::make({Polarization::StokesI,
                                                   Polarization::StokesQ,
                                                   Polarization::StokesU,
                                                   Polarization::StokesV}),
                     pol_out) {
  size_t nchan = 4;
  std::vector<float> actual(nchan);
  // Boost auto test case doesn't suport STL inside data::make
  std::map<PolarizationEnum, std::vector<float>> expected{
      {Polarization::StokesI, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesQ, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesU, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesV, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
  };

  ExtractWeights(actual.data(), 0, nchan,
                 std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                            Polarization::LR, Polarization::LL},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected[pol_out].begin(),
                                expected[pol_out].end());
}

BOOST_DATA_TEST_CASE(
    extractweights_linear_pol_to_linear,
    boost::unit_test::data::make({Polarization::XX, Polarization::XY,
                                  Polarization::YX, Polarization::YY}),
    pol_out) {
  size_t nchan = 4;
  std::vector<float> actual(nchan);
  // Boost auto test case doesn't suport STL inside data::make
  std::map<PolarizationEnum, std::vector<float>> expected{
      {Polarization::XX, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::XY, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::YX, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::YY, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
  };

  ExtractWeights(actual.data(), 0, nchan,
                 std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected[pol_out].begin(),
                                expected[pol_out].end());
}

BOOST_DATA_TEST_CASE(
    extractweights_circular_pol_to_circular,
    boost::unit_test::data::make({Polarization::RR, Polarization::RL,
                                  Polarization::LR, Polarization::LL}),
    pol_out) {
  size_t nchan = 4;
  std::vector<float> actual(nchan);
  // Boost auto test case doesn't suport STL inside data::make
  std::map<PolarizationEnum, std::vector<float>> expected{
      {Polarization::RR, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::RL, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::LR, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::LL, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
  };

  ExtractWeights(actual.data(), 0, nchan,
                 std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                            Polarization::LR, Polarization::LL},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected[pol_out].begin(),
                                expected[pol_out].end());
}

BOOST_AUTO_TEST_CASE(extractweights_linear_pol_to_instrumental) {
  size_t nchan = 4;
  size_t pol_per_file = 4;
  std::vector<float> actual(nchan * pol_per_file);
  std::vector<float> expected(test_weights);
  size_t idx = 0;
  for (float& weight : expected) {
    weight *= 4;
    // Channel 3 is flagged
    weight *= idx >= 8 && idx < 12 ? 0.0f : 1.0f;
    idx++;
  }

  PolarizationEnum pol_out = Polarization::Instrumental;

  ExtractWeights(actual.data(), 0, nchan,
                 std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extractweights_circular_pol_to_instrumental) {
  size_t nchan = 4;
  size_t pol_per_file = 4;
  std::vector<float> actual(nchan * pol_per_file);
  std::vector<float> expected(test_weights);
  size_t idx = 0;
  for (float& weight : expected) {
    weight *= 4;
    // Channel 3 is flagged
    weight *= idx >= 8 && idx < 12 ? 0.0f : 1.0f;
    idx++;
  }

  PolarizationEnum pol_out = Polarization::Instrumental;

  ExtractWeights(actual.data(), 0, nchan,
                 std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                            Polarization::LR, Polarization::LL},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extractweight_linear_pol_to_diag_instrumental) {
  size_t nchan = 4;
  size_t pol_per_file = 2;
  std::vector<float> actual(nchan * pol_per_file);
  std::vector<float> expected{0.4f, 0.4f, 0.8f, 0.8f, 0.0f, 0.0f, 1.6f, 1.6f};
  PolarizationEnum pol_out = Polarization::DiagonalInstrumental;

  ExtractWeights(actual.data(), 0, nchan,
                 std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extractweight_circular_pol_to_diag_instrumental) {
  size_t nchan = 4;
  size_t pol_per_file = 2;
  std::vector<float> actual(nchan * pol_per_file);
  std::vector<float> expected{0.4f, 0.4f, 0.8f, 0.8f, 0.0f, 0.0f, 1.6f, 1.6f};
  PolarizationEnum pol_out = Polarization::DiagonalInstrumental;

  ExtractWeights(actual.data(), 0, nchan,
                 std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                            Polarization::LR, Polarization::LL},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_SUITE_END()
