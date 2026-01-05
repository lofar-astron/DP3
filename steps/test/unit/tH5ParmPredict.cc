// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "steps/H5ParmPredict.h"

#include <boost/test/unit_test.hpp>

#include "base/DP3.h"

#include "steps/Predict.h"
#include "common/ParameterSet.h"

#include "H5ParmFixture.h"
#include "tPredict.h"

using dp3::steps::H5ParmPredict;
using dp3::steps::Step;
using schaapcommon::h5parm::H5Parm;

BOOST_AUTO_TEST_SUITE(h5parm_predict)

BOOST_FIXTURE_TEST_CASE(fields, dp3::steps::test::H5ParmFixture) {
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkymodel);
  parset.add("applycal.parmdb", kParmDb);
  parset.add("applycal.correction", kSoltabName);
  const H5ParmPredict h5predict(parset, "");

  // H5ParmPredict creates Predict steps which have internal OnePredict
  // sub-steps as next step.
  const dp3::steps::Predict predict(parset, "",
                                    {dp3::steps::test::kPredictDirection});

  const dp3::common::Fields predict_required =
      dp3::base::GetChainRequiredFields(
          std::make_shared<dp3::steps::Predict>(predict));

  // TODO(AST-1033) Determine Predict provided fields using generic DP3
  // functions.
  const dp3::common::Fields predict_provided =
      predict.getNextStep()->getProvidedFields();
  BOOST_TEST(h5predict.getRequiredFields() == predict_required);
  BOOST_TEST(h5predict.getProvidedFields() == predict_provided);
}

BOOST_AUTO_TEST_SUITE_END()
