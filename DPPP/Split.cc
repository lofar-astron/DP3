// Split.cc: DPPP step class to Split visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "Exceptions.h"
#include "Split.h"
#include "DPRun.h"

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"
#include "../Common/StreamUtil.h"

#include <stddef.h>

#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using namespace casacore;

namespace DP3 {
namespace DPPP {

Split::Split(DPInput* input, const ParameterSet& parset, const string& prefix)
    : itsAddedToMS(false) {
  itsReplaceParms = parset.getStringVector(prefix + "replaceparms");
  vector<vector<string> > replaceParmValues;  // For each of the parameters, the
                                              // values for each substep
  replaceParmValues.resize(itsReplaceParms.size());

  vector<vector<string> >::iterator replaceParmValueIt =
      replaceParmValues.begin();
  unsigned int numSteps = 0;
  for (const string& replaceParm : itsReplaceParms) {
    vector<string> parmValues = parset.getStringVector(replaceParm);
    *(replaceParmValueIt++) = parmValues;
    if (numSteps > 0) {
      if (parmValues.size() != numSteps)
        throw Exception(
            "Each parameter in replaceparms should have the same number of "
            "items (expected " +
            std::to_string(numSteps) + ", got " +
            std::to_string(parmValues.size()) + " for step " + replaceParm);
    } else {
      numSteps = parmValues.size();
    }
  }

  // Make a shallow copy to work around constness of parset
  ParameterSet parsetCopy(parset);

  // Create the substeps
  unsigned int numParameters = itsReplaceParms.size();
  for (unsigned int i = 0; i < numSteps; ++i) {
    for (unsigned int j = 0; j < numParameters; ++j) {
      parsetCopy.replace(itsReplaceParms[j], replaceParmValues[j][i]);
    }
    DPStep::ShPtr firstStep = DPRun::makeSteps(parsetCopy, prefix, input);
    firstStep->setPrevStep(this);
    itsSubsteps.push_back(firstStep);
  }
}

Split::~Split() {}

void Split::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;

  vector<DPStep::ShPtr>::iterator it;
  for (it = itsSubsteps.begin(); it != itsSubsteps.end(); ++it) {
    (*it)->setInfo(infoIn);
  }
}

void Split::show(std::ostream& os) const {
  os << "Split " << itsName << '\n'
     << "  replace parameters:" << itsReplaceParms << '\n';
  // Show the steps.
  for (unsigned int i = 0; i < itsSubsteps.size(); ++i) {
    os << "Split substep " << (i + 1) << " of " << itsSubsteps.size() << endl;
    DPStep::ShPtr step = itsSubsteps[i];
    DPStep::ShPtr lastStep;
    while (step) {
      step->show(os);
      lastStep = step;
      step = step->getNextStep();
    }
  }
}

void Split::showTimings(std::ostream& os, double duration) const {
  for (unsigned int i = 0; i < itsSubsteps.size(); ++i) {
    DPStep::ShPtr step = itsSubsteps[i];
    while (step) {
      step->showTimings(os, duration);
      step = step->getNextStep();
    }
  }
}

bool Split::process(const DPBuffer& bufin) {
  for (unsigned int i = 0; i < itsSubsteps.size(); ++i) {
    itsSubsteps[i]->process(bufin);
  }
  return false;
}

void Split::finish() {
  // Let the next steps finish.
  for (unsigned int i = 0; i < itsSubsteps.size(); ++i) {
    itsSubsteps[i]->finish();
  }
}

void Split::addToMS(const string& msName) {
  if (itsAddedToMS) {
    getPrevStep()->addToMS(msName);
  } else {
    itsAddedToMS = true;
    for (auto& subStep : itsSubsteps) {
      DPStep::ShPtr step, lastStep;
      step = subStep;
      while (step) {
        lastStep = step;
        step = step->getNextStep();
      }
      lastStep->addToMS(msName);
    }
  }
}
}  // namespace DPPP
}  // namespace DP3
