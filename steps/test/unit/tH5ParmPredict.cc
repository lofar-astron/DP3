// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../H5ParmPredict.h"

#include <boost/test/unit_test.hpp>

#include "../../Predict.h"
#include "../../../common/ParameterSet.h"

#include "mock/MockInput.h"
#include "H5ParmFixture.h"
#include "tPredict.h"

using dp3::steps::H5ParmPredict;
using dp3::steps::Step;
using schaapcommon::h5parm::H5Parm;

BOOST_AUTO_TEST_SUITE(h5parm_predict)

BOOST_FIXTURE_TEST_CASE(fields, dp3::steps::test::H5ParmFixture) {
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSourceDB);
  parset.add("applycal.parmdb", kParmDb);
  parset.add("applycal.correction", kSoltabName);
  const H5ParmPredict h5predict(&input, parset, "");

  // H5ParmPredict creates Predict steps which have internal OnePredict
  // sub-steps as next step.
  const dp3::steps::Predict predict(input, parset, "",
                                    {dp3::steps::test::kPredictDirection});
  // TODO(AST-1032) Determine Predict fields using generic DP3 functions.
  const dp3::common::Fields predict_required =
      predict.getNextStep()->getRequiredFields();
  const dp3::common::Fields predict_provided =
      predict.getNextStep()->getProvidedFields();
  BOOST_TEST(h5predict.getRequiredFields() == predict_required);
  BOOST_TEST(h5predict.getProvidedFields() == predict_provided);
}

BOOST_AUTO_TEST_SUITE_END()
