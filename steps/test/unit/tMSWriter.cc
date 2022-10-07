// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MSWriter.h"

#include <boost/test/unit_test.hpp>

#include "mock/MockInput.h"

using dp3::steps::MSWriter;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(mswriter)

BOOST_AUTO_TEST_CASE(insert_number_in_filename) {
  BOOST_CHECK_EQUAL(MSWriter::InsertNumberInFilename("", 1982), "-1982");
  BOOST_CHECK_EQUAL(MSWriter::InsertNumberInFilename("test", 0), "test-000");
  BOOST_CHECK_EQUAL(MSWriter::InsertNumberInFilename("afile", 12), "afile-012");
  BOOST_CHECK_EQUAL(MSWriter::InsertNumberInFilename("another-file", 123),
                    "another-file-123");
  BOOST_CHECK_EQUAL(
      MSWriter::InsertNumberInFilename("file-with-extension.ms", 42),
      "file-with-extension-042.ms");
  BOOST_CHECK_EQUAL(MSWriter::InsertNumberInFilename("/and.also/path.ms", 42),
                    "/and.also/path-042.ms");
}

BOOST_AUTO_TEST_CASE(fields) {
  const dp3::common::Fields kAlwaysWritten =
      Step::kDataField | Step::kFlagsField | Step::kWeightsField |
      Step::kUvwField;

  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;

  {
    const MSWriter writer(input, "test_mswriter.ms", parset, "");
    BOOST_TEST(writer.getRequiredFields() ==
               (kAlwaysWritten | Step::kFullResFlagsField));
    BOOST_TEST(writer.getProvidedFields() == dp3::common::Fields());
  }

  {
    parset.add("writefullresflag", "false");
    const MSWriter writer(input, "test_mswriter.ms", parset, "");
    BOOST_TEST(writer.getRequiredFields() == kAlwaysWritten);
    BOOST_TEST(writer.getProvidedFields() == dp3::common::Fields());
  }
}

BOOST_AUTO_TEST_SUITE_END()
