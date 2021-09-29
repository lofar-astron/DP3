// DemixWorker.h: Demixer helper class processing a time chunk
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to average in time and/or freq
/// @author Ger van Diepen

#ifndef DPPP_DEMIXWORKER_H
#define DPPP_DEMIXWORKER_H

#include "DemixInfo.h"
#include "DPBuffer.h"
#include "Patch.h"
#include "EstimateNew.h"

#include "../steps/InputStep.h"
#include "../steps/PhaseShift.h"
#include "../steps/Filter.h"

#include "../parmdb/ParmDB.h"

#include <EveryBeam/station.h>
#include <aocommon/matrix2x2.h>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MeasureHolder.h>
#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>

namespace dp3 {

namespace base {

/// @brief Demixer helper class processing a time chunk

/// DemixWorker::process processes a single time window (say, 2 minutes).
/// It predicts the A-team and target sources to determine which sources
/// have to be taken into account and which antennae have to be solved for.
/// Multiple DemixWorker::process can be executed in parallel by the parent
/// class DemixerNew.
//
/// Each DemixWorker object references a DemixInfo object containing the
/// general info and parameters.

class DemixWorker {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DemixWorker(steps::InputStep*, const string& prefix, const DemixInfo& info,
              const DPInfo& dpinfo, int workernr);

  /// Process the data in the input buffers and store the result in the
  /// output buffers.
  void process(const DPBuffer* bufin, unsigned int nbufin, DPBuffer* bufout,
               std::vector<double>* solutions, unsigned int chunkNr);

  /// Get the number of solves.
  unsigned int nSolves() const { return itsNrSolves; }
  /// Get the number of converged solves.
  unsigned int nConverged() const { return itsNrConverged; }
  /// Get the total nr of iterations used.
  unsigned int nIterations() const { return itsNrIter; }
  /// Get the number of times no demix was needed.
  unsigned int nNoDemix() const { return itsNrNoDemix; }
  unsigned int nIncludeStrongTarget() const { return itsNrIncludeStrongTarget; }
  unsigned int nIncludeCloseTarget() const { return itsNrIncludeCloseTarget; }
  unsigned int nIgnoreTarget() const { return itsNrIgnoreTarget; }
  unsigned int nDeprojectTarget() const { return itsNrDeprojectTarget; }
  /// Get nr of times a source was demixed.
  const casacore::Vector<unsigned int>& nsourcesDemixed() const {
    return itsNrSourcesDemixed;
  }
  /// Get nr of times a station was demixed.
  const casacore::Vector<unsigned int>& nstationsDemixed() const {
    return itsNrStationsDemixed;
  }
  /// Get nr of times a station/source was demixed.
  const casacore::Matrix<unsigned int>& statSourceDemixed() const {
    return itsStatSourceDemixed;
  }
  const casacore::Matrix<double>& amplSubtrMean() const {
    return itsAmplSubtrMean;
  }
  const casacore::Matrix<double>& amplSubtrM2() const { return itsAmplSubtrM2; }
  const casacore::Matrix<size_t>& amplSubtrNr() const { return itsAmplSubtrNr; }

  /// Get the timings of the various processing steps.
  ///@{
  double getTotalTime() const { return itsTimer.getElapsed(); }
  double getCoarseTime() const { return itsTimerCoarse.getElapsed(); }
  double getPhaseShiftTime() const { return itsTimerPhaseShift.getElapsed(); }
  double getDemixTime() const { return itsTimerDemix.getElapsed(); }
  double getPredictTime() const { return itsTimerPredict.getElapsed(); }
  double getSolveTime() const { return itsTimerSolve.getElapsed(); }
  double getSubtractTime() const { return itsTimerSubtract.getElapsed(); }
  ///@}

 private:
  /// Setup the demix processing steps for this piece of data.
  /// It fills itsFirstSteps, etc. for the sources to be demixed.
  /// It also determines how to handle the target (include,deproject,ignore).
  void setupDemix(unsigned int chunkNr);

  /// Find the median ampltitude for the selected baselines.
  /// It uses itsTmpAmpl as temporary buffer.
  float findMedian(const casacore::Cube<float>& ampl, const bool* selbl);

  /// Average the baseline UVWs in bufin and split them into UVW per station.
  /// It returns the number of time averages.
  unsigned int avgSplitUVW(const DPBuffer* bufin, unsigned int nbufin,
                           unsigned int ntimeAvg,
                           const std::vector<unsigned int>& selbl);

  /// Predict the target StokesI amplitude.
  /// It applies the beam at each target patch.
  void predictTarget(const std::vector<Patch::ConstPtr>& patchList,
                     unsigned int ntime, double time, double timeStep);

  /// Predict the StokesI amplitude of the Ateam patches and determine
  /// which antennae and sources to use when demixing.
  /// It applies the beam at each patch center.
  void predictAteam(const std::vector<Patch::ConstPtr>& patchList,
                    unsigned int ntime, double time, double timeStep);

  /// Add the StokesI of itsPredictVis to ampl.
  void addStokesI(casacore::Matrix<float>& ampl);

  /// Calculate the beam for demix resolution and apply to itsPredictVis.
  /// If apply==False, nothing is done.
  void applyBeam(double time, const Direction& direction, bool apply);

  /// Calculate the beam for the given sky direction and frequencies.
  /// Apply it to the data.
  /// If apply==False, nothing is done.
  void applyBeam(double time, const Direction& direction, bool apply,
                 const casacore::Vector<double>& chanFreqs, dcomplex* data);

  /// Convert a direction to ITRF.
  everybeam::vector3r_t dir2Itrf(const casacore::MDirection&);

  /// Calculate the StokesI amplitude from the predicted visibilities.
  /// (0.5 * (XX+YY))
  void calcStokesI(casacore::Matrix<float>& ampl);

  /// Simply average the data if no demixing needs to bedone.
  void average(const DPBuffer* bufin, unsigned int nbufin, DPBuffer* bufout);

  /// Add the decorrelation factor contribution for each time slot.
  void addFactors(const DPBuffer& newBuf,
                  casacore::Array<casacore::DComplex>& factorBuf);

  /// Calculate the decorrelation factors by averaging them.
  /// Apply the P matrix to deproject the sources without a model.
  void makeFactors(const casacore::Array<casacore::DComplex>& bufIn,
                   casacore::Array<casacore::DComplex>& bufOut,
                   const casacore::Cube<float>& weightSums,
                   unsigned int nChanOut, unsigned int nChanAvg);

  /// Deproject the sources without a model.
  void deproject(casacore::Array<casacore::DComplex>& factors,
                 std::vector<steps::MultiResultStep*> avgResults,
                 unsigned int resultIndex);

  /// Do the demixing.
  void handleDemix(DPBuffer* bufout, std::vector<double>* solutions,
                   double time, double timeStep);

  /// Solve gains and subtract sources.
  void demix(std::vector<double>* solutions, double time, double timeStep);

  /// Add amplitude subtracted to the arrays for mean and stddev.
  void addMeanM2(const std::vector<float>& sourceAmpl, unsigned int src);

  /// Merge the data of the selected baselines from the subtract buffer
  /// into the full buffer.
  void mergeSubtractResult();

  int itsWorkerNr;
  const DemixInfo* itsMix;
  std::vector<steps::PhaseShift*> itsOrigPhaseShifts;
  /// Phase shift and average steps for demix.
  std::vector<steps::Step::ShPtr> itsOrigFirstSteps;
  /// Result of phase shifting and averaging the directions of interest
  /// at the demix resolution.
  std::vector<steps::MultiResultStep*> itsAvgResults;
  std::vector<steps::PhaseShift*> itsPhaseShifts;
  std::vector<steps::Step::ShPtr> itsFirstSteps;
  steps::Step::ShPtr itsAvgStepSubtr;
  steps::Filter itsFilter;
  std::shared_ptr<steps::Filter> itsFilterSubtr;
  /// Result of averaging the target at the subtract resolution.
  std::shared_ptr<steps::MultiResultStep> itsAvgResultFull;
  std::shared_ptr<steps::MultiResultStep> itsAvgResultSubtr;
  /// The sources to demix (excluding target).
  std::vector<Patch::ConstPtr> itsDemixList;
  // TODO: unique_ptr?
  std::shared_ptr<everybeam::telescope::Telescope> telescope_;
  /// Measure objects unique to this worker (thread).
  /// This is needed because they are not thread-safe.
  casacore::MPosition itsArrayPos;
  casacore::MDirection itsDelayCenter;
  casacore::MDirection itsTileBeamDir;

  /// Variables set by setupDemix and used by handleDemix.
  unsigned int itsNDir;
  unsigned int itsNModel;
  unsigned int itsNSubtr;
  bool itsIgnoreTarget;
  bool itsIncludeTarget;
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

  /// Variables for conversion of directions to ITRF.
  casacore::MeasFrame itsMeasFrame;
  casacore::MDirection::Convert itsMeasConverter;
  std::vector<aocommon::MC2x2> itsBeamValues;  ///< [nst,nch]

  /// Indices telling which Ateam sources to use.
  std::vector<unsigned int> itsSrcSet;
  casacore::Cube<double> itsStationUVW;    ///< UVW per station
  casacore::Matrix<double> itsAvgUVW;      ///< temp buffer
  casacore::Cube<dcomplex> itsPredictVis;  ///< temp buffer
  /// #nfreq x #bl x #time StokesI amplitude per A-source.
  std::vector<casacore::Cube<float>> itsAteamAmpl;
  /// #bl x #src telling if baseline has sufficient Ateam flux.
  casacore::Matrix<bool> itsAteamAmplSel;
  /// #nfreq x #bl x #time StokesI amplitude of target.
  casacore::Cube<float> itsTargetAmpl;
  /// Temporary buffer to determine medians.
  std::vector<float> itsTmpAmpl;
  /// Per A-source and for target the min and max amplitude.
  ///@{
  std::vector<double> itsAteamMinAmpl;
  std::vector<double> itsAteamMaxAmpl;
  ///@}
  double itsTargetMinAmpl;
  double itsTargetMaxAmpl;
  /// Per A-source the stations to use (matching the minimum amplitude).
  std::vector<std::vector<unsigned int>> itsStationsToUse;
  casacore::Block<bool> itsSolveStation;  ///< solve station i?
  /// Per station and source the index in the unknowns vector.
  /// Note there are 8 unknowns (4 pol, ampl/phase) per source/station.
  std::vector<std::vector<int>> itsUnknownsIndex;
  /// The estimater (solver).
  EstimateNew itsEstimate;
  /// Variables for the predict.
  ///@{
  casacore::Matrix<double> itsUVW;
  std::vector<casacore::Cube<dcomplex>> itsModelVisDemix;
  std::vector<casacore::Cube<dcomplex>> itsModelVisSubtr;
  unsigned int itsNTimeOut;
  unsigned int itsNTimeOutSubtr;
  unsigned int itsTimeIndex;
  std::vector<float> itsObservedAmpl;
  std::vector<float> itsSourceAmpl;
  std::vector<float> itsSumSourceAmpl;
  ///@}
  /// Statistics
  ///@{
  unsigned int itsNrSolves;
  unsigned int itsNrConverged;
  unsigned int itsNrIter;
  unsigned int itsNrNoDemix;
  unsigned int itsNrIncludeStrongTarget;
  unsigned int itsNrIncludeCloseTarget;
  unsigned int itsNrIgnoreTarget;
  unsigned int itsNrDeprojectTarget;
  ///@}
  /// Nr of times a source is demixed.
  casacore::Vector<unsigned int> itsNrSourcesDemixed;
  /// Nr of times a station is demixed.
  casacore::Vector<unsigned int> itsNrStationsDemixed;
  /// Nr of times a source/station is demixed.
  casacore::Matrix<unsigned int> itsStatSourceDemixed;
  /// Average amplitude subtracted for middle channel [nbl,nsrc]
  casacore::Matrix<double> itsAmplSubtrMean;
  /// M2n to calculate stddev online in stable way (see Wikipedia)
  casacore::Matrix<double> itsAmplSubtrM2;
  /// N for mean/stddev amplitude calculations.
  casacore::Matrix<size_t> itsAmplSubtrNr;

  common::NSTimer itsTimer;
  common::NSTimer itsTimerCoarse;
  common::NSTimer itsTimerPhaseShift;
  common::NSTimer itsTimerDemix;
  common::NSTimer itsTimerPredict;
  common::NSTimer itsTimerSolve;
  common::NSTimer itsTimerSubtract;
};

}  // namespace base
}  // namespace dp3

#endif
