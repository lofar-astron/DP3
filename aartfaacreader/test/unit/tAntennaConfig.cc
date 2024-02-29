#include <aartfaacreader/AntennaConfig.h>
#include <base/RcuMode.h>
#include <casacore/measures/Measures/MPosition.h>
#include <filesystem>
#include <fstream>
#include <array>
#include <boost/filesystem.hpp>  // for the unique_path generation
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

namespace {
bool ComparePositions(const casacore::MPosition &left,
                      const casacore::MPosition &right,
                      const casacore::MPosition &reference) {
  return left.getValue() == (right.getValue() + reference.getValue());
}

casacore::MPosition PositionFromITRFCoordinates(double x, double y, double z) {
  return casacore::MPosition(casacore::MVPosition{x, y, z},
                             casacore::MPosition::ITRF);
}

void TestPosition(const casacore::MPosition &left,
                  const casacore::MPosition &right,
                  const casacore::MPosition &reference) {
  BOOST_CHECK(ComparePositions(left, right, reference));
}

void TestArrayPositions(
    const std::vector<casacore::MPosition> &positions,
    const std::vector<casacore::MPosition> &expected_position,
    const casacore::MPosition &reference_antenna) {
  BOOST_CHECK_EQUAL(positions.size(), expected_position.size());
  for (size_t i = 0; i < expected_position.size(); i++) {
    TestPosition(positions[i], expected_position[i], reference_antenna);
  }
}

const int kAxesDimension = 3 * 3;
void TestAxes(const std::array<double, kAxesDimension> &axes,
              const std::array<double, kAxesDimension> &expected_axes) {
  for (int i = 0; i < kAxesDimension; i++) {
    BOOST_CHECK_EQUAL(axes[i], expected_axes[i]);
  }
}
const std::string configuration_example =
    "\
#\n\
# AntennaPositions for AARTFAAC-12 LBA_OUTER antennas\n\
# ITRF2005 target_date = 2015.5\n\
# Created: 2023-10-16 14:35:38\n\
#\n\
\n\
NORMAL_VECTOR LBA\n\
3 [   0.598753   0.072099   0.797682 ]\n\
\n\
ROTATION_MATRIX LBA\n\
3 x 3 [\n\
 -0.1195950000  -0.7919540000   0.5987530000 \n\
  0.9928230000  -0.0954190000   0.0720990000 \n\
  0.0000330000   0.6030780000   0.7976820000 \n\
]\n\
\n\
LBA\n\
3 [ 3826577.462000000 461022.624000000 5064892.526 ]\n\
4 x 2 x 3 [\n\
    2.38300  -17.79500   -0.18000       2.38300  -17.79500   -0.18000\n\
    0.95600  -20.19400    1.10800       0.95600  -20.19400    1.10800\n\
  -10.83100  -14.47100    9.43800     -10.83100  -14.47100    9.43800\n\
  -15.87100    0.64500   11.85500     -15.87100    0.64500   11.85500\n\
]\n\
\n\
HBA\n\
3 [ 3826577.462000000 461022.624000000 5064892.526 ]\n\
4 x 2 x 3 [\n\
   -2.38300  -17.79500   -0.18000      -2.38300  -17.79500   -0.18000\n\
   -0.95600  -20.19400    1.10800      -0.95600  -20.19400    1.10800\n\
  +10.83100  -14.47100    9.43800     +10.83100  -14.47100    9.43800\n\
  +15.87100    0.64500   11.85500     +15.87100    0.64500   11.85500\n\
]\n\
\n\
NORMAL_VECTOR HBA0\n\
3 [   0.598753   0.072099   0.797682 ]\n\
\n\
ROTATION_MATRIX HBA0\n\
3 x 3 [\n\
 -0.1195950000  -0.7919540000   0.5987530000 \n\
  0.9928230000  -0.0954190000   0.0720990000 \n\
  0.0000330000   0.6030780000   0.7976820000 \n\
]\n\
\n\
HBA0\n\
3 [ 0.000000000 0.000000000 0.000 ]\n\
\n\
NORMAL_VECTOR HBA1\n\
3 [   0.598753   0.072099   0.797682 ]\n\
\n\
ROTATION_MATRIX HBA1\n\
3 x 3 [\n\
 -0.1195950000  -0.7919540000   0.5987530000 \n\
  0.9928230000  -0.0954190000   0.0720990000 \n\
  0.0000330000   0.6030780000   0.7976820000 \n\
]\n\
\n\
HBA1\n\
3 [ 0.000000000 0.000000000 0.000 ]\n\
\n";

struct AntennaConfigFixture {
  AntennaConfigFixture() {
    const std::filesystem::path tmp_path =
        std::filesystem::temp_directory_path() /
        boost::filesystem::unique_path("tmp%%%%%%%.config").string();
    path = tmp_path.string();
    std::ofstream config;
    config.open(path);
    config << configuration_example;
    config.close();
  };

  ~AntennaConfigFixture() { std::filesystem::remove(path); };

  std::string path;
};

const casacore::MPosition kLbaReference{PositionFromITRFCoordinates(
    3826577.462000000, 461022.62400000, 5064892.526)};

const std::vector<casacore::MPosition> kLbaShifts = {
    PositionFromITRFCoordinates(2.38300, -17.79500, -0.18000),
    PositionFromITRFCoordinates(0.95600, -20.19400, 1.10800),
    PositionFromITRFCoordinates(-10.83100, -14.47100, 9.43800),
    PositionFromITRFCoordinates(-15.87100, 0.64500, 11.85500)};

const casacore::MPosition kHbaReference{PositionFromITRFCoordinates(
    3826577.462000000, 461022.62400000, 5064892.526)};

const std::vector<casacore::MPosition> kHbaShifts = {
    PositionFromITRFCoordinates(-2.38300, -17.79500, -0.18000),
    PositionFromITRFCoordinates(-0.95600, -20.19400, 1.10800),
    PositionFromITRFCoordinates(+10.83100, -14.47100, 9.43800),
    PositionFromITRFCoordinates(+15.87100, 0.64500, 11.85500)};

const std::array<double, kAxesDimension> kLbaAxes{
    -0.1195950000, -0.7919540000, 0.5987530000, 0.9928230000, -0.0954190000,
    0.0720990000,  0.0000330000,  0.6030780000, 0.7976820000};

const std::array<double, kAxesDimension> kHba0Axes{
    -0.1195950000, -0.7919540000, 0.5987530000, 0.9928230000, -0.0954190000,
    0.0720990000,  0.0000330000,  0.6030780000, 0.7976820000};

const std::array<double, kAxesDimension> kHba1Axes{
    -0.1195950000, -0.7919540000, 0.5987530000, 0.9928230000, -0.0954190000,
    0.0720990000,  0.0000330000,  0.6030780000, 0.7976820000};

const int kFirstLbaMode = 1;
const int kLastLbaMode = 4;
const int kFirstHbaMode = 5;
const int kLastHbaMode = 7;
}  // namespace
BOOST_AUTO_TEST_SUITE(aartfaacantennaconfig)

BOOST_FIXTURE_TEST_CASE(constructor_and_parsing, AntennaConfigFixture) {
  BOOST_REQUIRE_NO_THROW(
      dp3::aartfaacreader::AntennaConfig antennaConfig(path.c_str()));
};

BOOST_FIXTURE_TEST_CASE(get_lba_positions, AntennaConfigFixture) {
  dp3::aartfaacreader::AntennaConfig antennaConfig(path.c_str());

  std::vector<casacore::MPosition> positions = antennaConfig.GetLBAPositions();
  TestArrayPositions(positions, kLbaShifts, kLbaReference);
}

BOOST_FIXTURE_TEST_CASE(get_hba_positions, AntennaConfigFixture) {
  dp3::aartfaacreader::AntennaConfig antennaConfig(path.c_str());

  std::vector<casacore::MPosition> positions = antennaConfig.GetHBAPositions();
  TestArrayPositions(positions, kHbaShifts, kHbaReference);
}
BOOST_FIXTURE_TEST_CASE(get_lba_axes, AntennaConfigFixture) {
  dp3::aartfaacreader::AntennaConfig antennaConfig(path.c_str());

  std::array<double, kAxesDimension> lba_axes = antennaConfig.GetLBAAxes();
  TestAxes(lba_axes, kLbaAxes);
}
BOOST_FIXTURE_TEST_CASE(get_hba0_axes, AntennaConfigFixture) {
  dp3::aartfaacreader::AntennaConfig antennaConfig(path.c_str());

  std::array<double, kAxesDimension> hba0_axes = antennaConfig.GetHBA0Axes();
  TestAxes(hba0_axes, kHba0Axes);
}
BOOST_FIXTURE_TEST_CASE(get_hba1_axes, AntennaConfigFixture) {
  dp3::aartfaacreader::AntennaConfig antennaConfig(path.c_str());

  std::array<double, kAxesDimension> hba1_axes = antennaConfig.GetHBA1Axes();
  TestAxes(hba1_axes, kHba1Axes);
}

BOOST_FIXTURE_TEST_CASE(get_axes_from_mode, AntennaConfigFixture) {
  dp3::aartfaacreader::AntennaConfig antennaConfig(path.c_str());

  for (int i = kFirstLbaMode; i <= kLastLbaMode; i++) {
    TestAxes(antennaConfig.GetAxesFromMode(dp3::base::RcuMode::FromNumber(i)),
             kLbaAxes);
  }
  for (int i = kFirstHbaMode; i <= kLastHbaMode; i++) {
    TestAxes(antennaConfig.GetAxesFromMode(dp3::base::RcuMode::FromNumber(i)),
             kHba0Axes);
  }

  BOOST_CHECK_THROW(
      antennaConfig.GetAxesFromMode(dp3::base::RcuMode::FromNumber(0)),
      std::runtime_error);
}
BOOST_FIXTURE_TEST_CASE(get_array_from_mode, AntennaConfigFixture) {
  dp3::aartfaacreader::AntennaConfig antennaConfig(path.c_str());

  for (int i = kFirstLbaMode; i <= kLastLbaMode; i++) {
    TestArrayPositions(
        antennaConfig.GetArrayFromMode(dp3::base::RcuMode::FromNumber(i)),
        kLbaShifts, kLbaReference);
  }
  for (int i = kFirstHbaMode; i <= kLastHbaMode; i++) {
    TestArrayPositions(
        antennaConfig.GetArrayFromMode(dp3::base::RcuMode::FromNumber(i)),
        kHbaShifts, kHbaReference);
  }
  BOOST_CHECK_THROW(
      antennaConfig.GetArrayFromMode(dp3::base::RcuMode::FromNumber(0)),
      std::runtime_error);
}
BOOST_AUTO_TEST_SUITE_END()