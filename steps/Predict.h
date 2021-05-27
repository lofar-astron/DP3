// Predict.h: DPPP step class to predict visibilities from a source model
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to predict visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DPPP_PREDICT_H
#define DPPP_PREDICT_H

#include "ApplyBeam.h"
#include "ApplyCal.h"
#include "InputStep.h"

#include "../base/DPBuffer.h"
#include "../base/ModelComponent.h"
#include "../base/Patch.h"
#include "../base/PredictBuffer.h"
#include "../base/SourceDBUtil.h"

#include <EveryBeam/station.h>
#include <EveryBeam/common/types.h>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Quanta/MVEpoch.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <mutex>
#include <utility>

namespace aocommon {
class ThreadPool;
}  // namespace aocommon

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// This class is a Step class to predict visibilities with optionally beam

typedef std::pair<size_t, size_t> Baseline;
typedef std::pair<base::ModelComponent::ConstPtr, base::Patch::ConstPtr> Source;
typedef std::complex<double> dcomplex;

/// @brief DPPP step class to predict visibilities from a source model
class Predict : public ModelDataStep {
 public:
  typedef std::shared_ptr<Predict> ShPtr;

  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Predict(InputStep*, const common::ParameterSet&, const std::string& prefix);

  /// Constructor with explicit sourcelist
  Predict(InputStep*, const common::ParameterSet&, const std::string& prefix,
          const std::vector<string>& sourcePatterns);

  /// The actual constructor
  void init(InputStep*, const common::ParameterSet&, const std::string& prefix,
            const std::vector<string>& sourcePatterns);

  /// Set the applycal substep
  void setApplyCal(InputStep*, const common::ParameterSet&,
                   const std::string& prefix);

  /// Set the operation type
  void setOperation(const std::string& type);

  /// When multiple Predict steps are running in parallel from multiple threads,
  /// they require synchronisation. This is done with these two synchronisation
  /// structures. When multiple Predicts steps run serially (like currently in
  /// H5ParmPredict), this function should not be called, as otherwise they
  /// will synchronize needlessly.
  ///
  /// It is also possible to make the predict steps share the same threadpool
  /// without further synchronisation, by setting measuresMutex to nullptr.
  void setThreadData(aocommon::ThreadPool& pool, std::mutex* measuresMutex) {
    itsThreadPool = &pool;
    itsMeasuresMutex = measuresMutex;
  }

  void setPredictBuffer(std::shared_ptr<base::PredictBuffer> predictBuffer) {
    itsPredictBuffer = std::move(predictBuffer);
  }

  virtual ~Predict();

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

  /// Prepare the sources
  void setSources(const std::vector<string>& sourcePatterns);

  /// Return the direction of the first patch
  std::pair<double, double> GetFirstDirection() const override;

 private:
  void initializeThreadData();
  everybeam::vector3r_t dir2Itrf(const casacore::MDirection& dir,
                                 casacore::MDirection::Convert& measConverter);
  void addBeamToData(base::Patch::ConstPtr patch, double time,
                     const everybeam::vector3r_t& refdir,
                     const everybeam::vector3r_t& tiledir, size_t thread,
                     size_t nBeamValues, dcomplex* data0, bool stokesIOnly);
  InputStep* itsInput;
  std::string itsName;
  base::DPBuffer itsBuffer;
  std::string itsSourceDBName;
  bool itsCorrectFreqSmearing;
  std::string itsOperation;
  bool itsApplyBeam;
  bool itsUseChannelFreq;
  bool itsOneBeamPerPatch;
  /// If two sources are closer together than given by this setting, they
  /// will be grouped into one patch. Value is in arcsec; zero means don't
  /// group.
  double itsBeamProximityLimit;
  bool itsStokesIOnly;
  base::Position itsPhaseRef;
  bool itsMovingPhaseRef;

  bool itsDoApplyCal;
  ApplyCal itsApplyCalStep;
  std::shared_ptr<ResultStep>
      itsResultStep;  ///< For catching results from ApplyCal

  unsigned int itsDebugLevel;

  std::vector<Baseline> itsBaselines;

  /// Vector containing info on converting baseline uvw to station uvw
  std::vector<int> itsUVWSplitIndex;

  /// UVW coordinates per station (3 coordinates per station)
  casacore::Matrix<double> itsStationUVW;

  /// The info needed to calculate the station beams.
  std::shared_ptr<base::PredictBuffer> itsPredictBuffer;
  base::BeamCorrectionMode itsBeamMode;
  everybeam::ElementResponseModel itsElementResponseModel;
  std::vector<casacore::MeasFrame> itsMeasFrames;
  std::vector<casacore::MDirection::Convert> itsMeasConverters;

  std::string itsDirectionsStr;  ///< Definition of patches, to pass to applycal
  std::vector<base::Patch::ConstPtr> itsPatchList;
  std::vector<Source> itsSourceList;

  common::NSTimer itsTimer;
  common::NSTimer itsTimerPredict;

  aocommon::ThreadPool* itsThreadPool;
  std::mutex* itsMeasuresMutex;
};

}  // namespace steps
}  // namespace dp3

#endif