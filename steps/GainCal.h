// GainCal.h: DPPP step class to calibrate (direction independent) gains
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to apply a calibration correction to the data
/// @author Tammo Jan Dijkema

#ifndef DPPP_GAINCAL_H
#define DPPP_GAINCAL_H

#include "ApplyBeam.h"
#include "InputStep.h"
#include "Predict.h"
#include "UVWFlagger.h"

#include "../base/DPBuffer.h"
#include "../base/PhaseFitter.h"
#include "../base/BaselineSelection.h"
#include "../base/GainCalAlgorithm.h"
#include "../base/Patch.h"
#include "../base/SourceDBUtil.h"

#include "../parmdb/Parm.h"
#include "../parmdb/ParmFacade.h"
#include "../parmdb/ParmSet.h"

#include <aocommon/parallelfor.h>
#include <aocommon/threadpool.h>

#include <EveryBeam/station.h>
#include <EveryBeam/common/types.h>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Quanta/MVEpoch.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <schaapcommon/h5parm/h5parm.h>

/// Convince HDF5 to use new API, even when system is configured to use 1.6 API
#define H5Acreate_vers 2
#define H5Tarray_create_vers 2
#define H5Dcreate_vers 2
#define H5Gcreate_vers 2
#include <H5Cpp.h>

namespace dp3 {
namespace base {
enum class CalType;
}
namespace common {
class ParameterSet;
}

namespace steps {

/// @brief This class is a Step class to calibrate (direction independent)
/// gains.
class GainCal final : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  GainCal(InputStep& input, const common::ParameterSet& parset,
          const std::string& prefix);

  virtual ~GainCal();

  /// Process the data. It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer&) override;

  virtual void finish() override;

  virtual void updateInfo(const base::DPInfo&) override;

  virtual void show(std::ostream&) const override;

  virtual void showTimings(std::ostream&, double duration) const override;

  virtual bool modifiesData() const override { return itsApplySolution; }

  /// Make a soltab with the given type
  static std::vector<schaapcommon::h5parm::SolTab> makeSolTab(
      schaapcommon::h5parm::H5Parm& h5parm, base::CalType caltype,
      std::vector<schaapcommon::h5parm::AxisInfo>& axes);

 private:
  /// Perform gaincal (polarized or unpolarized)
  void calibrate();

  /// Check for scalar mode
  static bool scalarMode(base::CalType caltype);

  /// Check for diagonal mode
  static bool diagonalMode(base::CalType caltype);

  /// Apply the solution
  void applySolution(base::DPBuffer& buf,
                     const casacore::Cube<casacore::Complex>& invsol);

  /// Invert solution (for applying it)
  casacore::Cube<casacore::Complex> invertSol(
      const casacore::Cube<casacore::Complex>& sol);

  /// Fills the matrices itsVis and itsMVis
  void fillMatrices(casacore::Complex* model, casacore::Complex* data,
                    float* weight, const casacore::Bool* flag);

  /// Initialize the parmdb
  void initParmDB();

  /// Get parmdbname from itsMode
  std::string parmName();

  /// Determine which stations are used
  void setAntennaUsed();

  /// Write out the solutions of the current parameter chunk
  /// (timeslotsperparmupdate) Variant for writing ParmDB
  void writeSolutionsParmDB(double startTime);

  /// Write out the solutions of the current parameter chunk
  /// (timeslotsperparmupdate) Variant for writing H5Parm
  void writeSolutionsH5Parm(double startTime);

  InputStep& itsInput;
  std::string itsName;
  std::vector<base::DPBuffer> itsBuf;
  bool itsUseModelColumn;
  std::string itsModelColumnName;
  std::string itsParmDBName;
  bool itsUseH5Parm;
  std::shared_ptr<parmdb::ParmDB> itsParmDB;
  std::string itsParsetString;  ///< Parset, for logging in H5Parm

  base::CalType itsMode;

  unsigned int itsDebugLevel;
  bool itsDetectStalling;

  bool itsApplySolution;

  std::vector<casacore::Cube<casacore::Complex>>
      itsSols;  ///< for every timeslot, nCr x nSt x nFreqCells
  std::vector<casacore::Matrix<double>>
      itsTECSols;  ///< for every timeslot, 2 x nSt (alpha and beta)
  std::vector<double> itsFreqData;  ///< Mean frequency for every freqcell

  std::vector<std::unique_ptr<PhaseFitter>> itsPhaseFitters;  ///< Length nSt

  std::vector<base::GainCalAlgorithm> iS;

  UVWFlagger itsUVWFlagStep;
  std::shared_ptr<ResultStep>
      itsDataResultStep;  ///< Result step for data after UV-flagging

  aocommon::ParallelFor<size_t> itsParallelFor;
  aocommon::ThreadPool itsThreadPool;

  /// The series of sub-steps ends with itsResultStep.
  std::unique_ptr<Step> itsFirstSubStep;
  std::shared_ptr<ResultStep> itsResultStep;
  bool itsApplyBeamToModelColumn;

  base::BaselineSelection itsBaselineSelection;  ///< Filter
  casacore::Vector<bool> itsSelectedBL;   ///< Vector (length nBl) telling which
                                          ///< baselines are selected
  casacore::Vector<bool> itsAntennaUsed;  ///< Vector (length nSt) telling which
                                          ///< stations are solved for

  std::map<std::string, int> itsParmIdMap;  ///< -1 = new parm name

  unsigned int itsMaxIter;
  double itsTolerance;
  bool itsPropagateSolutions;
  unsigned int itsSolInt;  ///< Time cell size
  unsigned int itsNChan;   ///< Frequency cell size
  unsigned int itsNFreqCells;

  unsigned int itsTimeSlotsPerParmUpdate;
  unsigned int itsConverged;
  unsigned int itsNonconverged;
  unsigned int itsFailed;
  unsigned int itsStalled;
  std::vector<unsigned int>
      itsNIter;  ///< Total iterations made (for converged, stalled,
                 ///< nonconverged, failed)
  unsigned int itsStepInParmUpdate;  ///< Timestep within parameter update
  double itsChunkStartTime;          ///< First time value of chunk to be stored
  unsigned int itsStepInSolInt;      ///< Timestep within solint

  casacore::Array<casacore::DComplex>
      itsAllSolutions;  ///< Array that holds all solutions for all iterations

  base::FlagCounter itsFlagCounter;

  common::NSTimer itsTimer;
  common::NSTimer itsTimerPredict;
  common::NSTimer itsTimerSolve;
  common::NSTimer itsTimerPhaseFit;
  common::NSTimer itsTimerWrite;
  common::NSTimer itsTimerFill;
};

}  // namespace steps
}  // namespace dp3

#endif
