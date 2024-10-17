// GainCal.cc: DP3 step class to H5ParmPredict visibilities
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "H5ParmPredict.h"

#include <iostream>
#include <utility>
#include <vector>

#include <schaapcommon/h5parm/h5parm.h>

#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

using dp3::common::operator<<;

namespace dp3 {
namespace steps {

H5ParmPredict::H5ParmPredict(const common::ParameterSet& parset,
                             const std::string& prefix)
    : itsName(),
      itsPredictSteps(),
      itsResultStep(),
      itsH5ParmName(parset.getString(prefix + "applycal.parmdb")),
      itsDirections(parset.getStringVector(prefix + "directions",
                                           std::vector<std::string>())),
      itsTimer() {
  H5Parm h5parm = H5Parm(itsH5ParmName, false);
  std::string soltabName = parset.getString(prefix + "applycal.correction");
  if (soltabName == "fulljones") soltabName = "amplitude000";
  SolTab soltab = h5parm.GetSolTab(soltabName);

  std::vector<std::string> h5directions = soltab.GetStringAxis("dir");
  if (h5directions.empty())
    throw std::runtime_error("H5Parm has empty dir axis");

  std::string operation = parset.getString(prefix + "operation", "replace");

  if (itsDirections.empty()) {
    itsDirections = h5directions;
  } else {
    // Check that all specified directions are in the h5parm
    for (const std::string& dirStr : itsDirections) {
      if (find(h5directions.begin(), h5directions.end(), dirStr) ==
          h5directions.end()) {
        throw std::runtime_error("Direction " + dirStr + " not found in " +
                                 itsH5ParmName);
      }
    }
  }

  for (size_t i = 0; i < itsDirections.size(); ++i) {
    std::string directionStr = itsDirections[i];
    std::vector<std::string> directionVec;
    // each direction should be like '[patch1,patch2]'
    if (directionStr.size() <= 2 || directionStr[0] != '[' ||
        directionStr[directionStr.size() - 1] != ']')
      throw std::runtime_error(
          "Invalid direction string: expecting array between square brackets, "
          "e.g. [a, b]");
    directionVec = common::stringtools::tokenize(
        directionStr.substr(1, directionStr.size() - 2), ",");
    auto predictStep = std::make_shared<Predict>(parset, prefix, directionVec);

    if (operation == "replace" && i > 0) {
      predictStep->SetOperation("add");
    } else {
      predictStep->SetOperation(operation);
    }

    if (!itsPredictSteps.empty()) {
      itsPredictSteps.back()->setNextStep(predictStep);
    }
    itsPredictSteps.push_back(std::move(predictStep));
  }

  itsResultStep = std::make_shared<ResultStep>();
  itsPredictSteps.back()->setNextStep(itsResultStep);
}

void H5ParmPredict::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);

  itsPredictSteps.front()->setInfo(infoIn);
}

void H5ParmPredict::show(std::ostream& os) const {
  os << "H5ParmPredict " << itsName << '\n';
  os << "  H5Parm:     " << itsH5ParmName << '\n';
  os << "  directions: " << itsDirections << '\n';

  for (std::shared_ptr<Step> step = itsPredictSteps.front(); step;
       step = step->getNextStep()) {
    step->show(os);
  }
}

void H5ParmPredict::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " H5ParmPredict " << itsName << '\n';
}

bool H5ParmPredict::process(std::unique_ptr<DPBuffer> buffer) {
  itsTimer.start();

  itsPredictSteps.front()->process(std::move(buffer));
  buffer = itsResultStep->take();

  itsTimer.stop();
  getNextStep()->process(std::move(buffer));
  return false;
}

void H5ParmPredict::finish() {
  // Let the next steps finish.
  itsPredictSteps.front()->finish();
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
