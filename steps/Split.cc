// Split.cc: DPPP step class to Split visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "Split.h"

#include "../base/Exceptions.h"
#include "../base/DP3.h"

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/StreamUtil.h"

#include <stddef.h>

#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

Split::Split(InputStep* input, const common::ParameterSet& parset,
             const string& prefix)
    : itsAddedToMS(false) {
  itsReplaceParms = parset.getStringVector(prefix + "replaceparms");
  // For each of the parameters, the values for each substep
  std::vector<std::vector<string>> replaceParmValues(itsReplaceParms.size());

  // numSplits, the number of 'new streams' that the data are split into is
  // determined from the replaced parameters.
  size_t numSplits = 0;

  for (size_t parmIndex = 0; parmIndex != itsReplaceParms.size(); ++parmIndex) {
    const std::vector<string> parmValues =
        parset.getStringVector(itsReplaceParms[parmIndex]);
    replaceParmValues[parmIndex] = parmValues;
    if (numSplits == 0) {
      numSplits = parmValues.size();
    } else {
      if (parmValues.size() != numSplits)
        throw Exception(
            "Each parameter in replaceparms should have the same number of "
            "items (expected " +
            std::to_string(numSplits) + ", got " +
            std::to_string(parmValues.size()) + " for step " +
            itsReplaceParms[parmIndex]);
    }
  }

  // Make a shallow copy to work around constness of parset
  common::ParameterSet parsetCopy(parset);

  // Create the substeps
  const size_t numParameters = itsReplaceParms.size();
  for (size_t i = 0; i != numSplits; ++i) {
    for (size_t j = 0; j != numParameters; ++j) {
      parsetCopy.replace(itsReplaceParms[j], replaceParmValues[j][i]);
    }
    Step::ShPtr firstStep =
        base::DP3::makeStepsFromParset(parsetCopy, prefix, "steps", *input,
                                       true, steps::Step::MsType::kRegular);
    firstStep->setPrevStep(this);

    itsSubsteps.push_back(std::move(firstStep));
  }
}

Split::~Split() {}

void Split::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;

  for (Step::ShPtr& step : itsSubsteps) {
    step->setInfo(infoIn);
  }
}

void Split::show(std::ostream& os) const {
  os << "Split " << itsName << '\n'
     << "  replace parameters:" << itsReplaceParms << '\n';
  // Show the steps.
  for (unsigned int i = 0; i < itsSubsteps.size(); ++i) {
    os << "Split substep " << (i + 1) << " of " << itsSubsteps.size() << '\n';
    Step::ShPtr step = itsSubsteps[i];
    Step::ShPtr lastStep;
    while (step) {
      step->show(os);
      lastStep = step;
      step = step->getNextStep();
    }
  }
}

void Split::showTimings(std::ostream& os, double duration) const {
  for (unsigned int i = 0; i < itsSubsteps.size(); ++i) {
    Step::ShPtr step = itsSubsteps[i];
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
      Step::ShPtr step, lastStep;
      step = subStep;
      while (step) {
        lastStep = step;
        step = step->getNextStep();
      }
      lastStep->addToMS(msName);
    }
  }
}
}  // namespace steps
}  // namespace dp3
