// H5ParmPredict.h: DPPP step class to H5ParmPredict visibilities from a source
// model
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to H5ParmPredict visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DPPP_H5ParmPredict_H
#define DPPP_H5ParmPredict_H

#include "InputStep.h"
#include "Predict.h"

#include "../base/DPBuffer.h"
#include "../base/ModelComponent.h"
#include "../base/Patch.h"
#include "../base/PredictBuffer.h"

#include <aocommon/threadpool.h>

#include <utility>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// This class is a Step class to H5ParmPredict visibilities with optionally
/// beam

typedef std::pair<size_t, size_t> Baseline;
typedef std::pair<base::ModelComponent::ConstPtr, base::Patch::ConstPtr> Source;

/// @brief DPPP step class to H5ParmPredict visibilities from a source model
class H5ParmPredict : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  H5ParmPredict(InputStep*, const common::ParameterSet&, const string& prefix);

  virtual ~H5ParmPredict();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

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
