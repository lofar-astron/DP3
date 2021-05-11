// OneApplyCal.h: DPPP step class to apply a calibration correction to the data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to apply a calibration correction to the data
/// @author Tammo Jan Dijkema

#ifndef DPPP_ONEAPPLYCAL_H
#define DPPP_ONEAPPLYCAL_H

#include "InputStep.h"

#include "../base/DPBuffer.h"
#include "../base/FlagCounter.h"

#include "../parmdb/ParmFacade.h"
#include "../parmdb/ParmSet.h"
#include "../parmdb/Parm.h"

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>

#include <mutex>

using schaapcommon::h5parm::JonesParameters;

namespace dp3 {
namespace steps {
/// @brief DPPP step class to apply a calibration correction to the data

/// This class is a Step class applying calibration parameters to the data.
/// It only applies one correction.

class OneApplyCal : public Step {
 public:
  /// Define the shared pointer for this type.
  typedef std::shared_ptr<OneApplyCal> ShPtr;

  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  OneApplyCal(InputStep*, const common::ParameterSet&,
              const std::string& prefix, const std::string& defaultPrefix,
              bool substep = false, std::string predictDirection = "");

  virtual ~OneApplyCal();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer& buffer);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  bool invert() { return itsInvert; }

 private:
  /// Read parameters from the associated parmdb and store them in
  /// itsJonesParameters
  void updateParmsParmDB(const double bufStartTime);

  /// Read parameters from the associated h5 and store them in
  /// itsJonesParameters
  void updateParmsH5(const double bufStartTime);

  /// If needed, show the flag counts.
  virtual void showCounts(std::ostream&) const;

  void initDataArrays();

  /// Check the number of polarizations in the parmdb or h5parm
  unsigned int nPol(const std::string& parmName);

  /// Replace values by NaN on places where weight is zero
  static void applyFlags(std::vector<double>& values,
                         const std::vector<double>& weights);

  InputStep* itsInput;
  base::DPBuffer itsBuffer;
  std::string itsName;
  std::string itsParmDBName;
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
};

}  // namespace steps
}  // namespace dp3

#endif
