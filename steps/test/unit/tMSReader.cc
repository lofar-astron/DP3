// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MSReader.h"

#include <vector>

#include <xtensor/xstrided_view.hpp>
#include <xtensor/xview.hpp>

#include <boost/test/unit_test.hpp>

#include <casacore/tables/Tables/ArrayColumn.h>

#include <dp3/base/DPBuffer.h>

#include "../../../common/ParameterSet.h"

#include "../../ResultStep.h"

using dp3::steps::MSReader;

BOOST_AUTO_TEST_SUITE(msreader)

BOOST_AUTO_TEST_CASE(provided_fields) {
  const dp3::common::Fields kFieldsToRead(
      dp3::common::Fields::Single::kFullResFlags);
  const dp3::common::Fields kWeightsField(
      dp3::common::Fields::Single::kWeights);
  // TODO(AST-1061): Use an MS that can be autoweighted, and re-add
  // 'autoweight' to kAutoweightSettings.
  const casacore::MeasurementSet ms("tNDPPP_tmp.MS");
  const dp3::common::ParameterSet parset;
  const std::vector<std::string> kAutoweightSettings = {/*"autoweight",*/
                                                        "forceautoweight"};

  MSReader msreader(ms, parset, "");
  BOOST_TEST(msreader.getProvidedFields() == dp3::common::Fields());
  msreader.setFieldsToRead(kFieldsToRead);
  BOOST_TEST(msreader.getProvidedFields() == kFieldsToRead);

  for (const std::string& setting : kAutoweightSettings) {
    dp3::common::ParameterSet parset_autoweight;
    parset_autoweight.add(setting, "true");
    MSReader autoweight(ms, parset_autoweight, "");
    BOOST_TEST(autoweight.getProvidedFields() == kWeightsField);
    autoweight.setFieldsToRead(kFieldsToRead);
    BOOST_TEST(autoweight.getProvidedFields() ==
               (kFieldsToRead | kWeightsField));
  }
}

BOOST_AUTO_TEST_CASE(fill_full_res_flags_available_in_ms) {
  const std::string ms_name = "tNDPPP-generic.MS";
  const casacore::MeasurementSet ms(ms_name);
  BOOST_REQUIRE(ms.tableDesc().isColumn("LOFAR_FULL_RES_FLAG"));
  const dp3::common::ParameterSet parset;

  // Read the LOFAR_FULL_RES_FLAG column using the process() function.
  std::unique_ptr<dp3::base::DPBuffer> expected_buffer;
  {
    MSReader reader(ms, parset, "");
    reader.setFieldsToRead(
        dp3::common::Fields(dp3::common::Fields::Single::kFullResFlags));
    auto result = std::make_shared<dp3::steps::ResultStep>();
    reader.setNextStep(result);
    reader.process(std::make_unique<dp3::base::DPBuffer>());
    expected_buffer = result->extract();
    BOOST_REQUIRE(expected_buffer);
  }

  // Read the LOFAR_FULL_RES_FLAG column using FillFullResFlags.
  dp3::base::DPBuffer buffer;
  buffer.setRowNrs(expected_buffer->getRowNrs());
  {
    MSReader reader(ms, parset, "");
    reader.FillFullResFlags(buffer);
  }

  // Compare the flags obtained with the FillFullResFlags function with the
  // flags read from the MS.
  BOOST_CHECK(buffer.GetFullResFlags() == expected_buffer->GetFullResFlags());
}

BOOST_AUTO_TEST_CASE(fill_full_res_flags_not_available_in_ms) {
  const int kNCorr = 4;
  const int kNChan = 16;
  const int kNBaseline = 6;
  const std::string ms_name = "tNDPPP_tmp.MS";
  const casacore::MeasurementSet ms(ms_name);
  BOOST_REQUIRE(!ms.tableDesc().isColumn("LOFAR_FULL_RES_FLAG"));

  dp3::common::ParameterSet parset;
  MSReader reader(ms, parset, "");

  dp3::base::DPBuffer buffer;

  // Set some random flags
  casacore::Cube<bool> flags(kNCorr, kNChan, kNBaseline);
  flags = false;
  flags(0, 0, 0) = true;
  flags(0, 0, 1) = true;
  buffer.setFlags(flags);

  reader.FillFullResFlags(buffer);

  // Create a view of the buffer's flags into expected_full_resolution_flags
  // with the expected shape (there is no polarization axis in the full
  // resolution flags, the XX polarization flag is taken)
  const auto xx_flags = xt::view(buffer.GetFlags(), xt::all(), xt::all(), 0);
  const auto expected_full_resolution_flags =
      xt::reshape_view(xx_flags, {kNBaseline, 1, kNChan});

  // Compare the flags obtained with the FillFullResFlags function with the
  // flags copied from the buffer's flags
  BOOST_CHECK(buffer.GetFullResFlags() == expected_full_resolution_flags);
}

BOOST_AUTO_TEST_SUITE_END()
