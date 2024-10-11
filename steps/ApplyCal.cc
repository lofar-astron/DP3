// ApplyCal.cc: DP3 step class that applies gain solutions to visibilities.
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "ApplyCal.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <xtensor/xview.hpp>

#include "../common/ParameterSet.h"
#include "../common/ParameterValue.h"
#include "../common/Timer.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace {

inline bool isfinite(casacore::DComplex val) {
  return casacore::isFinite(val.real()) && casacore::isFinite(val.imag());
}

inline void CheckBuffer(const DPBuffer& buffer, std::size_t baseline,
                        std::size_t channel) {
  assert(baseline < buffer.GetData().shape(0));
  assert(channel < buffer.GetData().shape(1));
  assert(4 == buffer.GetData().shape(2));

  assert(buffer.GetFlags().shape(0) == buffer.GetData().shape(0));
  assert(buffer.GetFlags().shape(1) == buffer.GetData().shape(1));
}

}  // namespace

namespace dp3 {
namespace steps {

ApplyCal::ApplyCal(const common::ParameterSet& parset,
                   const std::string& prefix, bool substep,
                   std::string predictDirection)
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

void ApplyCal::ApplyDiag(const aocommon::MC2x2FDiag& gain_a,
                         const aocommon::MC2x2FDiag& gain_b, DPBuffer& buffer,
                         unsigned int baseline, unsigned int channel,
                         bool update_weights, base::FlagCounter& flag_counter,
                         const std::string& direction) {
  CheckBuffer(buffer, baseline, channel);

  // If parameter is NaN or inf, do not apply anything and flag the data
  if (!(isfinite(gain_a[0]) && isfinite(gain_b[0]) && isfinite(gain_a[1]) &&
        isfinite(gain_b[1]))) {
    // Only update flagcounter for first correlation
    if (!buffer.GetFlags()(baseline, channel, 0)) {
      flag_counter.incrChannel(channel);
      flag_counter.incrBaseline(baseline);
    }
    xt::view(buffer.GetFlags(), baseline, channel, xt::all()).fill(true);
    return;
  }

  aocommon::MC2x2F visibilities(
      &buffer.GetData(direction)(baseline, channel, 0));
  visibilities = gain_a * visibilities * gain_b.HermTranspose();
  visibilities.AssignTo(&buffer.GetData(direction)(baseline, channel, 0));

  if (update_weights) {
    DPBuffer::WeightsType& weights = buffer.GetWeights();
    weights(baseline, channel, 0) /=
        std::norm(gain_a[0]) * std::norm(gain_b[0]);
    weights(baseline, channel, 1) /=
        std::norm(gain_a[0]) * std::norm(gain_b[1]);
    weights(baseline, channel, 2) /=
        std::norm(gain_a[1]) * std::norm(gain_b[0]);
    weights(baseline, channel, 3) /=
        std::norm(gain_a[1]) * std::norm(gain_b[1]);
  }
}

void ApplyCal::ApplyScalar(const std::complex<float>& gain_a,
                           const std::complex<float>& gain_b, DPBuffer& buffer,
                           unsigned int baseline, unsigned int channel,
                           bool update_weights,
                           base::FlagCounter& flag_counter) {
  CheckBuffer(buffer, baseline, channel);

  // If parameter is NaN or inf, do not apply anything and flag the data
  if (!(isfinite(gain_a) && isfinite(gain_b))) {
    // Only update flagcounter for first correlation
    if (!buffer.GetFlags()(baseline, channel, 0)) {
      flag_counter.incrChannel(channel);
      flag_counter.incrBaseline(baseline);
    }
    xt::view(buffer.GetFlags(), baseline, channel, xt::all()).fill(true);
    return;
  }

  aocommon::MC2x2F visibilities(&buffer.GetData()(baseline, channel, 0));
  visibilities = visibilities * (gain_a * std::conj(gain_b));
  visibilities.AssignTo(&buffer.GetData()(baseline, channel, 0));

  if (update_weights) {
    DPBuffer::WeightsType& weights = buffer.GetWeights();
    const float norm_gain_a_b = std::norm(gain_a) * std::norm(gain_b);
    weights(baseline, channel, 0) /= norm_gain_a_b;
    weights(baseline, channel, 1) /= norm_gain_a_b;
    weights(baseline, channel, 2) /= norm_gain_a_b;
    weights(baseline, channel, 3) /= norm_gain_a_b;
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

void ApplyCal::ApplyFull(const aocommon::MC2x2F& gain_a,
                         const aocommon::MC2x2F& gain_b, DPBuffer& buffer,
                         unsigned int baseline, unsigned int channel,
                         bool update_weights, base::FlagCounter& flag_counter,
                         const std::string& direction) {
  CheckBuffer(buffer, baseline, channel);

  // If parameter is NaN or inf, do not apply anything and flag the data
  bool anyinfnan = false;
  for (unsigned int corr = 0; corr < 4; ++corr) {
    if (!(isfinite(gain_a[corr]) && isfinite(gain_b[corr]))) {
      anyinfnan = true;
      break;
    }
  }
  if (anyinfnan) {
    // Only update flag counter for first correlation
    if (!buffer.GetFlags()(baseline, channel, 0)) {
      flag_counter.incrChannel(channel);
      flag_counter.incrBaseline(baseline);
    }
    xt::view(buffer.GetFlags(), baseline, channel, xt::all()).fill(true);
    return;
  }

  aocommon::MC2x2F visibilities(
      &buffer.GetData(direction)(baseline, channel, 0));
  visibilities = gain_a * visibilities.MultiplyHerm(gain_b);
  visibilities.AssignTo(&buffer.GetData(direction)(baseline, channel, 0));

  if (update_weights) {
    ApplyWeights(gain_a, gain_b, &buffer.GetWeights()(baseline, channel, 0));
  }
}

void ApplyCal::ApplyWeights(const aocommon::MC2x2F& gain_a,
                            const aocommon::MC2x2F& gain_b, float* weight) {
  float cov[4], normGainA[4], normGainB[4];
  for (unsigned int i = 0; i < 4; ++i) {
    cov[i] = 1. / weight[i];
    normGainA[i] = std::norm(gain_a[i]);
    normGainB[i] = std::norm(gain_b[i]);
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
