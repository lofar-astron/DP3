#include "../../DS9FacetFile.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ds9facetfile)

BOOST_AUTO_TEST_CASE(direction_comment) {
  std::string dir = DS9FacetFile::parseDirection(DS9FacetFile::Comment,
                                                 "text=CygA, color=green");

  BOOST_TEST(dir == "CygA");
}

BOOST_AUTO_TEST_CASE(direction_comment_reversed) {
  std::string dir = DS9FacetFile::parseDirection(DS9FacetFile::Comment,
                                                 "color=green, text=CygA");

  BOOST_TEST(dir == "CygA");
}

BOOST_AUTO_TEST_CASE(direction_empty) {
  std::string dir =
      DS9FacetFile::parseDirection(DS9FacetFile::Comment, "some other comment");

  BOOST_TEST(dir.empty());
}

BOOST_AUTO_TEST_CASE(direction_wrong_type) {
  std::string dir =
      DS9FacetFile::parseDirection(DS9FacetFile::Word, "text=CygA");

  BOOST_TEST(dir.empty());
}

BOOST_AUTO_TEST_SUITE_END()
