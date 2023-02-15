// OneApplyCal.h: DP3 step class to apply a calibration correction to the data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to apply a calibration correction to the data
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_ONEAPPLYCAL_H_
#define DP3_STEPS_ONEAPPLYCAL_H_

#include <mutex>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>

#include <dp3/steps/Step.h>
#include <dp3/base/DPBuffer.h>
#include "../base/FlagCounter.h"
#include "../common/Timer.h"
#include "../parmdb/ParmFacade.h"
#include "../parmdb/ParmSet.h"
#include "../parmdb/Parm.h"

using schaapcommon::h5parm::JonesParameters;

namespace dp3 {
namespace steps {
/// @brief DPPP step class to apply a calibration correction to the data

/// This class is a Step class applying calibration parameters to the data.
/// It only applies one correction.

class OneApplyCal : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  OneApplyCal(const common::ParameterSet&, const std::string& prefix,
              const std::string& defaultPrefix, bool substep = false,
              std::string predictDirection = "");

  ~OneApplyCal() override;

  common::Fields getRequiredFields() const override {
    return kDataField | kWeightsField | kFlagsField;
  }

  common::Fields getProvidedFields() const override {
    common::Fields fields = kDataField | kFlagsField;
    if (itsUpdateWeights) fields |= kWeightsField;
    return fields;
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  bool invert() { return itsInvert; }

 private:
  /// Read parameters from the associated parmdb and store them in
  /// itsJonesParameters
  void updateParmsParmDB(const double bufStartTime);

  /// Read parameters from the associated h5 and store them in
  /// itsJonesParameters
  void updateParmsH5(const double bufStartTime);

  /// If needed, show the flag counts.
  void showCounts(std::ostream&) const override;

  void initDataArrays();

  /// Check the number of polarizations in the parmdb or h5parm
  unsigned int nPol(const std::string& parmName);

  /// Replace values by NaN on places where weight is zero
  static void applyFlags(std::vector<double>& values,
                         const std::vector<double>& weights);

  std::string itsName;
  std::string itsParmDBName;
  // itsParmDBOnDisk specifies the existence of a parmdb on disk. If this is
  // false, it means the solutions are instead stored in the DPBuffer, which can
  // be fetched with DPBuffer::GetSolutions.
  bool itsParmDBOnDisk;
  bool itsUseH5Parm;
  std::string itsSolSetName;
  std::shared_ptr<parmdb::ParmFacade> itsParmDB;
  schaapcommon::h5parm::H5Parm itsH5Parm;
  std::string itsSolTabName;
  JonesParameters::MissingAntennaBehavior itsMissingAntennaBehavior;
  schaapcommon::h5parm::SolTab itsSolTab;
  schaapcommon::h5parm::SolTab itsSolTab2;  ///< in the case of full Jones, amp
                                            ///< and phase table need to be open
  JonesParameters::CorrectType itsCorrectType;
  bool itsInvert;
  JonesParameters::InterpolationType itsInterpolationType;
  unsigned int itsTimeSlotsPerParmUpdate;
  double itsSigmaMMSE;
  bool itsUpdateWeights;

  unsigned int itsCount;  ///< number of steps

  /// Expressions to search for in itsParmDB
  std::vector<casacore::String> itsParmExprs;

  /// itsJonesParameters contains the gridded parameters, first for all
  /// parameters (e.g. Gain:0:0 and Gain:1:1), next all antennas, next over freq
  /// * time as returned by ParmDB numparms, antennas, time x frequency
  std::unique_ptr<JonesParameters> itsJonesParameters;
  unsigned int itsTimeStep;  ///< time step within current chunk
  unsigned int itsNCorr;
  double itsTimeInterval;
  double itsLastTime;  ///< last time of current chunk
  base::FlagCounter itsFlagCounter;
  bool itsUseAP;  ///< use ampl/phase or real/imag
  hsize_t itsDirection;
  common::NSTimer itsTimer;

  static std::mutex theirHDF5Mutex;  ///< Prevent parallel access to HDF5
  void setCorrectType(std::vector<std::string>& solTabs);
};

}  // namespace steps
}  // namespace dp3

#endif
