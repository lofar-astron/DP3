#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <base/RcuMode.h>

BOOST_AUTO_TEST_SUITE(rcumode)

using dp3::base::RcuMode;

BOOST_AUTO_TEST_CASE(create_from_mode_number) {
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(0).mode, RcuMode::Unused);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(1).mode, RcuMode::LBAOuter10_90);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(2).mode, RcuMode::LBAOuter30_90);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(3).mode, RcuMode::LBAInner10_90);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(4).mode, RcuMode::LBAInner30_90);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(5).mode, RcuMode::HBA110_190);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(6).mode, RcuMode::HBA170_230);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(7).mode, RcuMode::HBA210_270);

  BOOST_CHECK_THROW(RcuMode::FromNumber(-1).Bandwidth(), std::runtime_error);

  BOOST_CHECK_THROW(RcuMode::FromNumber(50).Bandwidth(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(frequency_and_bandwidth_from_mode_number) {
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(1).CentralFrequency(), 50.);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(2).CentralFrequency(), 60.);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(3).CentralFrequency(), 50.);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(4).CentralFrequency(), 60.);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(5).CentralFrequency(), 150.);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(6).CentralFrequency(), 200.);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(7).CentralFrequency(), 240.);

  BOOST_CHECK_EQUAL(RcuMode::FromNumber(1).Bandwidth(), 195312.5);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(2).Bandwidth(), 195312.5);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(3).Bandwidth(), 195312.5);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(4).Bandwidth(), 195312.5);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(5).Bandwidth(), 195312.5);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(6).Bandwidth(), 156250.);
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(7).Bandwidth(), 195312.5);

  BOOST_CHECK_THROW(RcuMode::FromNumber(0).Bandwidth(), std::runtime_error);
  BOOST_CHECK_THROW(RcuMode::FromNumber(-1).Bandwidth(), std::runtime_error);
  BOOST_CHECK_THROW(RcuMode::FromNumber(50).Bandwidth(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(antenna_type_from_number) {
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(1).AntennaType(), "LBA");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(7).AntennaType(), "HBA");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(10).AntennaType(), "?");
}

BOOST_AUTO_TEST_CASE(string_from_code) {
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(0).ToString(), "unused");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(1).ToString(), "LBA_OUTER 10-90 MHz");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(2).ToString(), "LBA_OUTER 30-90 MHz");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(3).ToString(), "LBA_INNER 10-90 MHz");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(4).ToString(), "LBA_INNER 30-90 MHz");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(5).ToString(), "HBA 110-190 MHz");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(6).ToString(), "HBA 170-230 MHz");
  BOOST_CHECK_EQUAL(RcuMode::FromNumber(7).ToString(), "HBA 210-270 MHz");
}

BOOST_AUTO_TEST_SUITE_END()