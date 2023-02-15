// GainCal.cc: DPPP step class to ApplyCal visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "ApplyCal.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "../common/ParameterSet.h"
#include "../common/ParameterValue.h"
#include "../common/Timer.h"

inline bool isfinite(casacore::DComplex val) {
  return casacore::isFinite(val.real()) && casacore::isFinite(val.imag());
}

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

ApplyCal::ApplyCal(const common::ParameterSet& parset, const string& prefix,
                   bool substep, string predictDirection)
    : is_sub_step_(substep) {
  std::vector<std::string> subStepNames;
  common::ParameterValue namesPar(parset.getString(prefix + "steps", ""));

  if (namesPar.isVector()) {
    subStepNames = namesPar.getStringVector();
  } else {
    subStepNames.push_back(namesPar.getString());
  }

  std::vector<std::string>::const_iterator subStepNameIter;
  for (subStepNameIter = subStepNames.begin();
       subStepNameIter != subStepNames.end(); ++subStepNameIter) {
    string subStepName = (*subStepNameIter);
    string subStepPrefix;
    if (subStepName.empty()) {
      // No substeps given, use parameters of this step
      subStepPrefix = prefix;
    } else {
      // Substeps given, use named parameters like applycal.applySol.parmdb
      subStepPrefix = prefix + subStepName + ".";
    }
    apply_cals_.push_back(std::make_shared<OneApplyCal>(
        parset, subStepPrefix, prefix, substep, predictDirection));
  }

  Step::setNextStep(apply_cals_.front());
  for (unsigned int step = 0; step < apply_cals_.size() - 1; ++step) {
    apply_cals_[step]->setNextStep(apply_cals_[step + 1]);
  }
}

void ApplyCal::setNextStep(std::shared_ptr<Step> nextStep) {
  apply_cals_.back()->setNextStep(nextStep);
}

void ApplyCal::show(std::ostream& os) const {
  // If not a substep, show will be called by DPRun,
  // through the nextStep() mechanism
  if (is_sub_step_) {
    for (const std::shared_ptr<OneApplyCal>& apply_cal : apply_cals_) {
      apply_cal->show(os);
    }
  }
}

void ApplyCal::showTimings(std::ostream& os, double duration) const {
  if (is_sub_step_) {
    for (const std::shared_ptr<OneApplyCal>& apply_cal : apply_cals_) {
      apply_cal->showTimings(os, duration);
    }
  }
}

bool ApplyCal::process(std::unique_ptr<DPBuffer> buffer) {
  getNextStep()->process(std::move(buffer));
  return true;
}

void ApplyCal::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

void ApplyCal::applyDiag(const casacore::Complex* gainA,
                         const casacore::Complex* gainB, casacore::Complex* vis,
                         float* weight, bool* flag, unsigned int bl,
                         unsigned int chan, bool updateWeights,
                         base::FlagCounter& flagCounter) {
  // If parameter is NaN or inf, do not apply anything and flag the data
  if (!(isfinite(gainA[0]) && isfinite(gainB[0]) && isfinite(gainA[1]) &&
        isfinite(gainB[1]))) {
    // Only update flagcounter for first correlation
    if (!flag[0]) {
      flagCounter.incrChannel(chan);
      flagCounter.incrBaseline(bl);
    }
    for (unsigned int corr = 0; corr < 4; ++corr) {
      flag[corr] = true;
    }
    return;
  }

  vis[0] *= gainA[0] * conj(gainB[0]);
  vis[1] *= gainA[0] * conj(gainB[1]);
  vis[2] *= gainA[1] * conj(gainB[0]);
  vis[3] *= gainA[1] * conj(gainB[1]);

  if (updateWeights) {
    weight[0] /= norm(gainA[0]) * norm(gainB[0]);
    weight[1] /= norm(gainA[0]) * norm(gainB[1]);
    weight[2] /= norm(gainA[1]) * norm(gainB[0]);
    weight[3] /= norm(gainA[1]) * norm(gainB[1]);
  }
}

void ApplyCal::applyScalar(const casacore::Complex* gainA,
                           const casacore::Complex* gainB,
                           casacore::Complex* vis, float* weight, bool* flag,
                           unsigned int bl, unsigned int chan,
                           bool updateWeights, base::FlagCounter& flagCounter) {
  // If parameter is NaN or inf, do not apply anything and flag the data
  if (!(isfinite(gainA[0]) && isfinite(gainB[0]))) {
    // Only update flagcounter for first correlation
    if (!flag[0]) {
      flagCounter.incrChannel(chan);
      flagCounter.incrBaseline(bl);
    }
    for (unsigned int corr = 0; corr < 4; ++corr) {
      flag[corr] = true;
    }
    return;
  }

  vis[0] *= gainA[0] * conj(gainB[0]);
  vis[1] *= gainA[0] * conj(gainB[0]);
  vis[2] *= gainA[0] * conj(gainB[0]);
  vis[3] *= gainA[0] * conj(gainB[0]);

  if (updateWeights) {
    weight[0] /= norm(gainA[0]) * norm(gainB[0]);
    weight[1] /= norm(gainA[0]) * norm(gainB[0]);
    weight[2] /= norm(gainA[0]) * norm(gainB[0]);
    weight[3] /= norm(gainA[0]) * norm(gainB[0]);
  }
}

// Inverts complex 2x2 input matrix
template <typename NumType>
void ApplyCal::invert(std::complex<NumType>* v, NumType sigmaMMSE) {
  // Add the variance of the nuisance term to the elements on the diagonal.
  const NumType variance = sigmaMMSE * sigmaMMSE;
  std::complex<NumType> v0 = v[0] + variance;
  std::complex<NumType> v3 = v[3] + variance;
  // Compute inverse in the usual way.
  std::complex<NumType> invDet(NumType(1.0) / (v0 * v3 - v[1] * v[2]));
  v[0] = v3 * invDet;
  v[2] = v[2] * -invDet;
  v[1] = v[1] * -invDet;
  v[3] = v0 * invDet;
}
template void ApplyCal::invert(std::complex<double>* v, double sigmaMMSE);
template void ApplyCal::invert(std::complex<float>* v, float sigmaMMSE);

void ApplyCal::applyFull(const casacore::Complex* gainA,
                         const casacore::Complex* gainB, casacore::Complex* vis,
                         float* weight, bool* flag, unsigned int bl,
                         unsigned int chan, bool doUpdateWeights,
                         base::FlagCounter& flagCounter) {
  casacore::Complex gainAxvis[4];

  // If parameter is NaN or inf, do not apply anything and flag the data
  bool anyinfnan = false;
  for (unsigned int corr = 0; corr < 4; ++corr) {
    if (!(isfinite(gainA[corr]) && isfinite(gainB[corr]))) {
      anyinfnan = true;
      break;
    }
  }
  if (anyinfnan) {
    // Only update flag counter for first correlation
    if (!flag[0]) {
      flagCounter.incrChannel(chan);
      flagCounter.incrBaseline(bl);
    }
    for (unsigned int corr = 0; corr < 4; ++corr) {
      flag[corr] = true;
    }
    return;
  }

  // gainAxvis = gainA * vis
  for (unsigned int row = 0; row < 2; ++row) {
    for (unsigned int col = 0; col < 2; ++col) {
      gainAxvis[2 * row + col] =
          gainA[2 * row + 0] * casacore::Complex(vis[2 * 0 + col]) +
          gainA[2 * row + 1] * casacore::Complex(vis[2 * 1 + col]);
    }
  }

  // vis = gainAxvis * gainB^H
  for (unsigned int row = 0; row < 2; ++row) {
    for (unsigned int col = 0; col < 2; ++col) {
      vis[2 * row + col] =
          gainAxvis[2 * row + 0] * std::conj(gainB[2 * col + 0]) +
          gainAxvis[2 * row + 1] * std::conj(gainB[2 * col + 1]);
    }
  }

  if (doUpdateWeights) {
    applyWeights(gainA, gainB, weight);
  }
}

void ApplyCal::applyWeights(const casacore::Complex* gainA,
                            const casacore::Complex* gainB, float* weight) {
  float cov[4], normGainA[4], normGainB[4];
  for (unsigned int i = 0; i < 4; ++i) {
    cov[i] = 1. / weight[i];
    normGainA[i] = std::norm(gainA[i]);
    normGainB[i] = std::norm(gainB[i]);
  }

  weight[0] = cov[0] * (normGainA[0] * normGainB[0]) +
              cov[1] * (normGainA[0] * normGainB[1]) +
              cov[2] * (normGainA[1] * normGainB[0]) +
              cov[3] * (normGainA[1] * normGainB[1]);
  weight[0] = 1. / weight[0];

  weight[1] = cov[0] * (normGainA[0] * normGainB[2]) +
              cov[1] * (normGainA[0] * normGainB[3]) +
              cov[2] * (normGainA[1] * normGainB[2]) +
              cov[3] * (normGainA[1] * normGainB[3]);
  weight[1] = 1. / weight[1];

  weight[2] = cov[0] * (normGainA[2] * normGainB[0]) +
              cov[1] * (normGainA[2] * normGainB[1]) +
              cov[2] * (normGainA[3] * normGainB[0]) +
              cov[3] * (normGainA[3] * normGainB[1]);
  weight[2] = 1. / weight[2];

  weight[3] = cov[0] * (normGainA[2] * normGainB[2]) +
              cov[1] * (normGainA[2] * normGainB[3]) +
              cov[2] * (normGainA[3] * normGainB[2]) +
              cov[3] * (normGainA[3] * normGainB[3]);
  weight[3] = 1. / weight[3];
}

}  // namespace steps
}  // namespace dp3
