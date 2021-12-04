// Demixer.h: DPPP step class to subtract A-team sources
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to average in time and/or freq
/// @author Ger van Diepen

#ifndef DPPP_DEMIXER_H
#define DPPP_DEMIXER_H

#include "InputStep.h"
#include "PhaseShift.h"
#include "Filter.h"

#include "../base/Baseline.h"
#include "../base/DPBuffer.h"
#include "../base/Patch.h"

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

typedef std::vector<base::Patch::ConstPtr> PatchList;

/// @brief DPPP step class to subtract A-team sources
/// This class is a Step class to subtract the strong A-team sources.
/// It is based on the demixing.py script made by Bas vd Tol and operates
/// per time chunk as follows:
/// <ul>
///  <li> The data are phase-shifted and averaged for each source.
///  <li> Demixing is done using the combined results.
///  <li> For each source a BBS solve, smooth, and predict is done.
///  <li> The predicted results are subtracted from the averaged data.
/// </ul>

class Demixer : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Demixer(InputStep*, const common::ParameterSet&, const string& prefix);

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

  /// Show the counts.
  virtual void showCounts(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

 private:
  /// Add the decorrelation factor contribution for each time slot.
  void addFactors(const base::DPBuffer& newBuf,
                  casacore::Array<casacore::DComplex>& factorBuf);

  /// Calculate the decorrelation factors by averaging them.
  /// Apply the P matrix to deproject the sources without a model.
  void makeFactors(const casacore::Array<casacore::DComplex>& bufIn,
                   casacore::Array<casacore::DComplex>& bufOut,
                   const casacore::Cube<float>& weightSums,
                   unsigned int nChanOut, unsigned int nChanAvg);

  /// Do the demixing.
  void handleDemix();

  /// Deproject the sources without a model.
  void deproject(casacore::Array<casacore::DComplex>& factors,
                 unsigned int resultIndex);

  /// Solve gains and subtract sources.
  void demix();

  /// Export the solutions to a ParmDB.
  void dumpSolutions();

  /// Merge the data of the selected baselines from the subtract buffer
  /// into the full buffer.
  void mergeSubtractResult();

  InputStep* itsInput;
  string itsName;
  base::DPBuffer itsBufTmp;
  string itsSkyName;
  string itsInstrumentName;
  double itsDefaultGain;
  size_t itsMaxIter;
  base::BaselineSelection itsSelBL;
  Filter itsFilter;
  std::vector<std::shared_ptr<PhaseShift>> itsPhaseShifts;
  bool itsMovingPhaseRef;
  casacore::MeasFrame itsMeasFrame;
  /// Phase shift and average steps for demix.
  std::vector<Step::ShPtr> itsFirstSteps;
  /// Result of phase shifting and averaging the directions of interest
  /// at the demix resolution.
  std::vector<std::shared_ptr<MultiResultStep>> itsAvgResults;
  std::shared_ptr<Step> itsAvgStepSubtr;
  std::shared_ptr<Filter> itsFilterSubtr;
  /// Result of averaging the target at the subtract resolution.
  std::shared_ptr<MultiResultStep> itsAvgResultFull;
  std::shared_ptr<MultiResultStep> itsAvgResultSubtr;
  /// Ignore target in demixing?
  bool itsIgnoreTarget;
  /// Name of the target. Empty if no model is available for the target.
  string itsTargetSource;
  std::vector<string> itsSubtrSources;
  std::vector<string> itsModelSources;
  std::vector<string> itsExtraSources;
  std::vector<string> itsAllSources;
  //      std::vector<double>                        itsCutOffs;
  bool itsPropagateSolutions;
  unsigned int itsNDir;
  unsigned int itsNModel;
  unsigned int itsNStation;
  unsigned int itsNBl;
  unsigned int itsNCorr;
  unsigned int itsNChanIn;
  unsigned int itsNTimeIn;
  unsigned int itsNTimeDemix;
  unsigned int itsNChanAvgSubtr;
  unsigned int itsNTimeAvgSubtr;
  unsigned int itsNChanOutSubtr;
  unsigned int itsNTimeOutSubtr;
  unsigned int itsNTimeChunk;
  unsigned int itsNTimeChunkSubtr;
  unsigned int itsNChanAvg;
  unsigned int itsNTimeAvg;
  unsigned int itsNChanOut;
  unsigned int itsNTimeOut;
  double itsTimeIntervalAvg;

  bool itsUseLBFGS;  ///< if this is not false, use LBFGS solver instead of
                     ///< LSQfit.
  unsigned int
      itsLBFGShistory;       ///< the size of LBFGS memory(history), specified
                             ///< as a multiple of the size of parameter vector.
  double itsLBFGSrobustdof;  ///< the degrees of freedom used in robust noise
                             ///< model.

  /// Accumulator used for computing the demixing weights at the demix
  /// resolution. The shape of this buffer is #correlations x #channels
  /// x #baselines x #directions x #directions (fastest axis first).
  casacore::Array<casacore::DComplex> itsFactorBuf;
  /// Buffer of demixing weights at the demix resolution. Each Array is a
  /// cube of shape #correlations x #channels x #baselines of matrices of
  /// shape #directions x #directions.
  std::vector<casacore::Array<casacore::DComplex>> itsFactors;

  /// Accumulator used for computing the demixing weights. The shape of this
  /// buffer is #correlations x #channels x #baselines x #directions
  /// x #directions (fastest axis first).
  casacore::Array<casacore::DComplex> itsFactorBufSubtr;
  /// Buffer of demixing weights at the subtract resolution. Each Array is a
  /// cube of shape #correlations x #channels x #baselines of matrices of
  /// shape #directions x #directions.
  std::vector<casacore::Array<casacore::DComplex>> itsFactorsSubtr;

  PatchList itsPatchList;
  base::Direction itsPhaseRef;
  std::vector<base::Baseline> itsBaselines;
  std::vector<int> itsUVWSplitIndex;
  casacore::Vector<double> itsFreqDemix;
  casacore::Vector<double> itsFreqSubtr;
  std::vector<double> itsUnknowns;
  std::vector<double> itsPrevSolution;
  unsigned int itsTimeIndex;
  unsigned int itsNConverged;
  base::FlagCounter itsFlagCounter;

  common::NSTimer itsTimer;
  common::NSTimer itsTimerPhaseShift;
  common::NSTimer itsTimerDemix;
  common::NSTimer itsTimerSolve;
  common::NSTimer itsTimerDump;
};

}  // namespace steps
}  // namespace dp3

#endif
