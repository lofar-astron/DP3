// H5ParmPredict.h: DPPP step class to H5ParmPredict visibilities from a source
// model
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to H5ParmPredict visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_H5PARM_PREDICT_H_
#define DP3_STEPS_H5PARM_PREDICT_H_

#include <utility>

#include <aocommon/threadpool.h>

#include "InputStep.h"
#include "Predict.h"
#include "ResultStep.h"

#include <dp3/base/DP3.h>
#include <dp3/base/DPBuffer.h>

#include "../base/ModelComponent.h"
#include "../base/Patch.h"
#include "../base/PredictBuffer.h"

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// This class is a Step class to H5ParmPredict visibilities with optionally
/// beam

typedef std::pair<size_t, size_t> Baseline;
typedef std::pair<std::shared_ptr<const base::ModelComponent>,
                  std::shared_ptr<const base::Patch>>
    Source;

/// @brief DPPP step class to H5ParmPredict visibilities from a source model
class H5ParmPredict : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  H5ParmPredict(InputStep*, const common::ParameterSet&, const string& prefix);

  ~H5ParmPredict() override;

  common::Fields getRequiredFields() const override {
    // Combine the result of all sub steps.
    return base::DP3::GetChainRequiredFields(itsPredictSteps.front());
  }

  common::Fields getProvidedFields() const override {
    // Combine the result of all sub steps.
    common::Fields fields;
    for (std::shared_ptr<Step> step = itsPredictSteps.front(); step;
         step = step->getNextStep()) {
      fields |= step->getProvidedFields();
    }
    return fields;
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

 private:
  InputStep* itsInput;
  std::string itsName;
  base::DPBuffer itsBuffer;

  std::vector<std::shared_ptr<Predict>> itsPredictSteps;
  /// Buffers used by the predict steps. Normally, each Predict step
  /// allocates its own buffers. However, since the Predict steps run
  /// sequentially, they are told to share this single buffer to save memory.
  std::shared_ptr<base::PredictBuffer> itsPredictBuffer;
  std::shared_ptr<ResultStep> itsResultStep;

  std::string itsH5ParmName;
  std::vector<std::string> itsDirections;

  common::NSTimer itsTimer;
  aocommon::ThreadPool itsThreadPool;
};

}  // namespace steps
}  // namespace dp3

#endif
