// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MSUpdater.h"

#include <boost/test/unit_test.hpp>

#include "../../../common/ParameterSet.h"

#include "mock/MockInput.h"

using dp3::steps::MSUpdater;
using dp3::steps::Step;

namespace {

const std::string kInputDataColumnName = "input_data";
const std::string kInputFlagColumnName = "input_flags";
const std::string kInputWeightColumnName = "input_weights";

/**
 * Test input class that returns fixed column names.
 */
class TestInput : public dp3::steps::MockInput {
 public:
  const std::string& dataColumnName() const override {
    return kInputDataColumnName;
  }

  const std::string& flagColumnName() const override {
    return kInputFlagColumnName;
  }

  const std::string& weightColumnName() const override {
    return kInputWeightColumnName;
  }
};

void TestFields(const std::string& data_column_name,
                const std::string& flag_column_name,
                const std::string& weight_column_name,
                const dp3::common::Fields expected_required_fields) {
  TestInput input;
  dp3::common::ParameterSet parset;
  parset.add("datacolumn", data_column_name);
  parset.add("flagcolumn", flag_column_name);
  parset.add("weightcolumn", weight_column_name);
  const MSUpdater updater(&input, "test_msupdater.ms", parset, "");
  BOOST_TEST(updater.getRequiredFields() == expected_required_fields);
  BOOST_TEST(updater.getProvidedFields() == dp3::common::Fields());
}

}  // namespace

BOOST_AUTO_TEST_SUITE(msupdater)

BOOST_AUTO_TEST_CASE(fields) {
  TestFields(kInputDataColumnName, kInputFlagColumnName, kInputWeightColumnName,
             dp3::common::Fields());
  TestFields("custom_data", kInputFlagColumnName, kInputWeightColumnName,
             Step::kDataField);
  TestFields(kInputDataColumnName, "custom_flag", kInputWeightColumnName,
             Step::kFlagsField);
  TestFields(kInputDataColumnName, kInputFlagColumnName, "custom_weight",
             Step::kWeightsField);
  TestFields("all_columns", "have", "different_names",
             Step::kDataField | Step::kFlagsField | Step::kWeightsField);
}

BOOST_AUTO_TEST_SUITE_END()
