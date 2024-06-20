#include <boost/test/unit_test.hpp>
#include <steps/IDGImager.h>
#include <common/ParameterSet.h>
#include "../../../base/test/LoggerFixture.h"

BOOST_AUTO_TEST_SUITE(
    IDGImager, *boost::unit_test::fixture<dp3::base::test::LoggerFixture>())

#ifdef HAVE_IDG
BOOST_AUTO_TEST_CASE(Init) {
  dp3::common::ParameterSet parset;
  std::string step_name = "test_step_name";
  parset.add("test_step_name.image_size", "1234");
  parset.add("test_step_name.type", "imager");
  dp3::steps::IDGImager imager(parset, step_name + ".");

  BOOST_CHECK_EQUAL(imager.GetImageSize(), 1234);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
