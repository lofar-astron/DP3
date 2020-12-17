// OneApplyCal.h: DPPP step class to apply a calibration correction to the data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to apply a calibration correction to the data
/// @author Tammo Jan Dijkema

#ifndef DPPP_ONEAPPLYCAL_H
#define DPPP_ONEAPPLYCAL_H

#include "DPInput.h"
#include "DPBuffer.h"
// #include "H5Parm.h"
#include "FlagCounter.h"

#include "../ParmDB/ParmFacade.h"
#include "../ParmDB/ParmSet.h"
#include "../ParmDB/Parm.h"

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <schaapcommon/h5parm/h5parm.h>

#include <mutex>

namespace DP3 {
namespace DPPP {
/// @brief DPPP step class to apply a calibration correction to the data

/// This class is a DPStep class applying calibration parameters to the data.
/// It only applies one correction.

class OneApplyCal : public DPStep {
 public:
  /// Define the shared pointer for this type.
  typedef std::shared_ptr<OneApplyCal> ShPtr;

  enum class InterpolationType { NEAREST, LINEAR };

  enum CorrectType {
    GAIN,
    FULLJONES,
    TEC,
    CLOCK,
    ROTATIONANGLE,
    SCALARPHASE,
    PHASE,
    ROTATIONMEASURE,
    SCALARAMPLITUDE,
    AMPLITUDE
  };

  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  OneApplyCal(DPInput*, const ParameterSet&, const std::string& prefix,
              const std::string& defaultPrefix, bool substep = false,
              std::string predictDirection = "");

  virtual ~OneApplyCal();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const DPBuffer& buffer);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  bool invert() { return itsInvert; }

 private:
  /// Read parameters from the associated parmdb and store them in itsParms
  void updateParms(const double bufStartTime);

  /// If needed, show the flag counts.
  virtual void showCounts(std::ostream&) const;

  void initDataArrays();

  /// Check the number of polarizations in the parmdb or h5parm
  unsigned int nPol(const std::string& parmName);

  /// Replace values by NaN on places where weight is zero
  static void applyFlags(std::vector<double>& values,
                         const std::vector<double>& weights);

  static std::string correctTypeToString(CorrectType);
  static CorrectType stringToCorrectType(const string&);

  DPInput* itsInput;
  DPBuffer itsBuffer;
  string itsName;
  string itsParmDBName;
  bool itsUseH5Parm;
  string itsSolSetName;
  std::shared_ptr<BBS::ParmFacade> itsParmDB;
  schaapcommon::h5parm::H5Parm itsH5Parm;
  string itsSolTabName;
  schaapcommon::h5parm::SolTab itsSolTab;
  schaapcommon::h5parm::SolTab itsSolTab2;  ///< in the case of full Jones, amp
                                            ///< and phase table need to be open
  CorrectType itsCorrectType;
  bool itsInvert;
  InterpolationType itsInterpolationType;
  unsigned int itsTimeSlotsPerParmUpdate;
  double itsSigmaMMSE;
  bool itsUpdateWeights;

  unsigned int itsCount;  ///< number of steps

  /// Expressions to search for in itsParmDB
  std::vector<casacore::String> itsParmExprs;

  /// parameters, numparms, antennas, time x frequency
  casacore::Cube<casacore::DComplex> itsParms;
  unsigned int itsTimeStep;  ///< time step within current chunk
  unsigned int itsNCorr;
  double itsTimeInterval;
  double itsLastTime;  ///< last time of current chunk
  FlagCounter itsFlagCounter;
  bool itsUseAP;  ///< use ampl/phase or real/imag
  hsize_t itsDirection;
  NSTimer itsTimer;

  static std::mutex theirHDF5Mutex;  ///< Prevent parallel access to HDF5
};

}  // namespace DPPP
}  // namespace DP3

#endif
