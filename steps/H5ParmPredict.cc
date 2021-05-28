// GainCal.cc: DPPP step class to H5ParmPredict visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "H5ParmPredict.h"

#include "../base/Exceptions.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"
#include "../common/StringTools.h"
#include "../common/Timer.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <iostream>

#include <schaapcommon/h5parm/h5parm.h>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;
using dp3::base::PredictBuffer;

using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

using dp3::common::operator<<;

namespace dp3 {
namespace steps {

H5ParmPredict::H5ParmPredict(InputStep* input,
                             const common::ParameterSet& parset,
                             const string& prefix)
    : itsInput(input),
      itsName(),
      itsPredictSteps(),
      itsPredictBuffer(std::make_shared<PredictBuffer>()),
      itsResultStep(),
      itsH5ParmName(parset.getString(prefix + "applycal.parmdb")),
      itsDirections(
          parset.getStringVector(prefix + "directions", std::vector<string>())),
      itsTimer(),
      itsThreadPool() {
  H5Parm h5parm = H5Parm(itsH5ParmName, false);
  std::string soltabName = parset.getString(prefix + "applycal.correction");
  if (soltabName == "fulljones") soltabName = "amplitude000";
  SolTab soltab = h5parm.GetSolTab(soltabName);

  std::vector<string> h5directions = soltab.GetStringAxis("dir");
  if (h5directions.empty())
    throw std::runtime_error("H5Parm has empty dir axis");

  string operation = parset.getString(prefix + "operation", "replace");

  if (itsDirections.empty()) {
    itsDirections = h5directions;
  } else {
    // Check that all specified directions are in the h5parm
    for (const string& dirStr : itsDirections) {
      if (find(h5directions.begin(), h5directions.end(), dirStr) ==
          h5directions.end()) {
        throw Exception("Direction " + dirStr + " not found in " +
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
    auto predictStep =
        std::make_shared<Predict>(*input, parset, prefix, directionVec);

    if (operation == "replace" && i > 0) {
      predictStep->SetOperation("add");
    } else {
      predictStep->SetOperation(operation);
    }
    predictStep->SetThreadData(itsThreadPool, nullptr);
    predictStep->SetPredictBuffer(itsPredictBuffer);

    if (!itsPredictSteps.empty()) {
      itsPredictSteps.back()->setNextStep(predictStep);
    }
    itsPredictSteps.push_back(predictStep);
  }

  itsResultStep = std::make_shared<ResultStep>();
  itsPredictSteps.back()->setNextStep(itsResultStep);
}

H5ParmPredict::~H5ParmPredict() {}

void H5ParmPredict::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  info().setNeedVisData();
  info().setWriteData();

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

bool H5ParmPredict::process(const DPBuffer& bufin) {
  itsThreadPool.SetNThreads(getInfo().nThreads());

  itsTimer.start();
  itsBuffer.copy(bufin);
  itsInput->fetchUVW(bufin, itsBuffer, itsTimer);
  itsInput->fetchWeights(bufin, itsBuffer, itsTimer);

  itsPredictSteps.front()->process(itsBuffer);
  itsBuffer = itsResultStep->get();

  itsTimer.stop();
  getNextStep()->process(itsBuffer);
  return false;
}

void H5ParmPredict::finish() {
  // Let the next steps finish.
  itsPredictSteps.front()->finish();
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
