//# Demixer.cc: DPPP step class to subtract A-team sources
//# Copyright (C) 2011
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/Demixer.h>
#include <DPPP/Averager.h>
#include <DPPP/PhaseShift.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <ParmDB/Axis.h>
#include <ParmDB/SourceDB.h>
#include <Common/LofarLogger.h>
#include <Common/StreamUtil.h>
#include <Common/OpenMP.h>
#include <BBSKernel/MeasurementAIPS.h>

#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <scimath/Mathematics/MatrixMathLA.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MCPosition.h>
#include <casa/Quanta/MVAngle.h>

#include <iostream>
#include <iomanip>

using namespace casa;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    Demixer::Demixer (DPInput* input,
                      const ParSet& parset, const string& prefix)
      : itsInput         (input),
        itsName          (prefix),
        itsSkyName       (parset.getString(prefix+"skymodel", "sky")),
        itsInstrumentName(parset.getString(prefix+"instrumentmodel",
                                           "instrument")),
        itsElevCutoff    (getAngle(parset.getString(prefix+"elevationcutoff",
                                                    "0."))),
        itsBBSExpr       (*input, itsSkyName, itsInstrumentName, itsElevCutoff),
        itsTargetSource  (parset.getString(prefix+"targetsource", string())),
        itsSubtrSources  (parset.getStringVector (prefix+"subtractsources")),
        itsModelSources  (parset.getStringVector (prefix+"modelsources",
                                                  vector<string>())),
        itsExtraSources  (parset.getStringVector (prefix+"othersources",
                                                  vector<string>())),
///        itsJointSolve    (parset.getBool  (prefix+"jointsolve", true)),
        itsNTimeIn       (0),
        itsNChanAvgSubtr (parset.getUint  (prefix+"freqstep", 1)),
        itsNTimeAvgSubtr (parset.getUint  (prefix+"timestep", 1)),
        itsNTimeOutSubtr (0),
        itsNChanAvg      (parset.getUint  (prefix+"demixfreqstep",
                                           itsNChanAvgSubtr)),
        itsNTimeAvg      (parset.getUint  (prefix+"demixtimestep",
                                           itsNTimeAvgSubtr)),
        itsNTimeChunk    (parset.getUint  (prefix+"ntimechunk", 0)),
        itsNTimeOut      (0)
    {
      // Get and set solver options.
      itsSolveOpt.maxIter =
        parset.getUint  (prefix+"Solve.Options.MaxIter", 300);
      itsSolveOpt.epsValue =
        parset.getDouble(prefix+"Solve.Options.EpsValue", 1e-9);
      itsSolveOpt.epsDerivative =
        parset.getDouble(prefix+"Solve.Options.EpsDerivative", 1e-9);
      itsSolveOpt.colFactor =
        parset.getDouble(prefix+"Solve.Options.ColFactor", 1e-9);
      itsSolveOpt.lmFactor  =
        parset.getDouble(prefix+"Solve.Options.LMFactor", 1.0);
      itsSolveOpt.balancedEq =
        parset.getBool  (prefix+"Solve.Options.BalancedEqs", false);
      itsSolveOpt.useSVD  =
        parset.getBool  (prefix+"Solve.Options.UseSVD", true);
      itsBBSExpr.setOptions (itsSolveOpt);
      /// Maybe optionally a parset parameter directions to give the
      /// directions of unknown sources.
      /// Or make sources a vector of vectors like [name, ra, dec] where
      /// ra and dec are optional.

      // Default nr of time chunks is maximum number of threads.
      if (itsNTimeChunk == 0) {
        itsNTimeChunk = OpenMP::maxThreads();
      }
      // Check that time windows fit nicely.
      ASSERTSTR ((itsNTimeChunk * itsNTimeAvg) % itsNTimeAvgSubtr == 0,
                 "time window should fit final averaging integrally");
      itsNTimeChunkSubtr = (itsNTimeChunk * itsNTimeAvg) / itsNTimeAvgSubtr;
      // Collect all source names.
      itsNrModel = itsSubtrSources.size() + itsModelSources.size();
      itsNrDir   = itsNrModel + itsExtraSources.size() + 1;
      itsAllSources.reserve (itsNrDir);
      itsAllSources.insert (itsAllSources.end(),
                            itsSubtrSources.begin(), itsSubtrSources.end());
      itsAllSources.insert (itsAllSources.end(),
                            itsModelSources.begin(), itsModelSources.end());
      itsAllSources.insert (itsAllSources.end(),
                            itsExtraSources.begin(), itsExtraSources.end());
      itsAllSources.push_back (itsTargetSource);
      // If the target source is given, add it to the model.
      // Because the target source has to be the last direction, it means
      // that (for the time being) no extra sources can be given.
      if (! itsTargetSource.empty()) {
        itsNrModel++;
        ASSERTSTR (itsExtraSources.empty(), "Currently no extrasources can "
                   "be given if the targetsource is given");
      }
      // Size buffers.
      itsFactors.resize      (itsNTimeChunk);
      itsFactorsSubtr.resize (itsNTimeChunkSubtr);
      itsPhaseShifts.reserve (itsNrDir-1);
      itsFirstSteps.reserve  (itsNrDir);
      itsAvgResults.reserve  (itsNrDir);

      // Get the patch names and positions from the SourceDB table.
      SourceDB sdb(ParmDBMeta("casa", itsSkyName));
      // Create the steps for the sources to be removed.
      // Demixing consists of the following steps:
      // - phaseshift data to each demix source
      // - average data in each direction, also for original phasecenter.
      // - determine demix factors for all directions
      // - use BBS to predict and solve in each direction. It is possible to
      //   predict more directions than to solve (for strong sources in field).
      // - use BBS to subtract the solved sources using the demix factors.
      //   The averaging used here can be smaller than used when solving.
      for (uint i=0; i<itsNrDir-1; ++i) {
        // First make the phaseshift and average steps for each demix source.
        // The resultstep gets the result.
        // The phasecenter can be given in a parameter. Its name is the default.
        // Look up the source direction in the patch table.
        // If found, turn it into a vector of strings.
        vector<PatchInfo> patchInfo (sdb.getPatchInfo (-1, itsAllSources[i]));
        vector<string> sourceVec (1, itsAllSources[i]);
        if (! patchInfo.empty()) {
          sourceVec[0] = toString(patchInfo[0].getRa());
          sourceVec.push_back (toString(patchInfo[0].getDec()));
        }
        PhaseShift* step1 = new PhaseShift (input, parset,
                                            prefix + itsAllSources[i] + '.',
                                            sourceVec);
        itsFirstSteps.push_back (DPStep::ShPtr(step1));
        itsPhaseShifts.push_back (step1);
        DPStep::ShPtr step2 (new Averager(input, parset, prefix));
        step1->setNextStep (step2);
        MultiResultStep* step3 = new MultiResultStep(itsNTimeChunk);
        step2->setNextStep (DPStep::ShPtr(step3));
        // There is a single demix factor step which needs to get all results.
        itsAvgResults.push_back (step3);
      }

      // Now create the step to average the data themselves.
      DPStep::ShPtr targetAvg(new Averager(input, prefix,
                                           itsNChanAvg, itsNTimeAvg));
      itsFirstSteps.push_back (targetAvg);
      MultiResultStep* targetAvgRes = new MultiResultStep(itsNTimeChunk);
      targetAvg->setNextStep (DPStep::ShPtr(targetAvgRes));
      itsAvgResults.push_back (targetAvgRes);

      // Do the same for the subtract.
      itsAvgSubtr = DPStep::ShPtr (new Averager(input, prefix,
                                                itsNChanAvgSubtr,
                                                itsNTimeAvgSubtr));
      itsAvgResultSubtr = new MultiResultStep(itsNTimeChunk);
      itsAvgSubtr->setNextStep (DPStep::ShPtr(itsAvgResultSubtr));
    }

    Demixer::~Demixer()
    {
    }

    void Demixer::updateInfo (DPInfo& info)
    {
      // Get size info.
      itsNChanIn  = info.nchan();
      itsNrBl     = info.nbaselines();
      itsNrCorr   = info.ncorr();
      itsFactorBuf.resize (IPosition(4, itsNrBl, itsNChanIn, itsNrCorr,
                                     itsNrDir*(itsNrDir-1)/2));
      itsFactorBufSubtr.resize (IPosition(4, itsNrBl, itsNChanIn, itsNrCorr,
                                          itsNrDir*(itsNrDir-1)/2));
      itsTimeInterval = info.timeInterval();

      // Adapt averaging to available nr of channels and times.
      // Use a copy of the DPInfo, otherwise it is updated multiple times.
      DPInfo infocp(info);
      itsNTimeAvg = std::min (itsNTimeAvg, infocp.ntime());
      itsNChanAvg = infocp.update (itsNChanAvg, itsNTimeAvg);
      // Let the internal steps update their data.
      infocp = info;
      itsAvgSubtr->updateInfo (infocp);
      for (uint i=0; i<itsFirstSteps.size(); ++i) {
        infocp = info;
        DPStep::ShPtr step = itsFirstSteps[i];
        while (step) {
          step->updateInfo (infocp);
          step = step->getNextStep();
        }
        // Create the BBS model expression for sources with a model.
        if (i < itsNrModel) {
          itsBBSExpr.addModel (itsAllSources[i], infocp.phaseCenter());
        }
      }
      // Keep the averaged time interval.
      itsNChanOut = infocp.nchan();
      itsTimeIntervalAvg = infocp.timeInterval();
      // Update the info of this object.
      info.setNeedVisData();
      info.setNeedWrite();
      itsNTimeAvgSubtr = std::min (itsNTimeAvgSubtr, info.ntime());
      itsNChanAvgSubtr = info.update (itsNChanAvgSubtr, itsNTimeAvgSubtr);
      itsNChanOutSubtr = info.nchan();
      itsTimeIntervalSubtr = info.timeInterval();
      // Construct frequency axis for the demix and subtract averaging.
      itsFreqAxisDemix = makeFreqAxis (itsNChanAvg);
      itsFreqAxisSubtr = makeFreqAxis (itsNChanAvgSubtr);
    }

    void Demixer::show (std::ostream& os) const
    {
      os << "Demixer " << itsName << std::endl;
      os << "  skymodel:       " << itsSkyName << std::endl;
      os << "  instrumentmodel:" << itsInstrumentName << std::endl;
      os << "  targetsource:   " << itsTargetSource << std::endl;
      os << "  subtractsources:" << itsSubtrSources << std::endl;
      os << "  modelsources:   " << itsModelSources << std::endl;
      os << "  extrasources:   " << itsExtraSources << std::endl;
///      os << "  jointsolve:     " << itsJointSolve << std::endl;
      os << "  freqstep:       " << itsNChanAvgSubtr << std::endl;
      os << "  timestep:       " << itsNTimeAvgSubtr << std::endl;
      os << "  demixfreqstep:  " << itsNChanAvg << std::endl;
      os << "  demixtimestep:  " << itsNTimeAvg << std::endl;
      os << "  ntimechunk:     " << itsNTimeChunk << std::endl;
      os << "  Solve.Options.MaxIter:       " << itsSolveOpt.maxIter << endl;
      os << "  Solve.Options.EpsValue:      " << itsSolveOpt.epsValue << endl;
      os << "  Solve.Options.EpsDerivative: " << itsSolveOpt.epsDerivative << endl;
      os << "  Solve.Options.ColFactor:     " << itsSolveOpt.colFactor << endl;
      os << "  Solve.Options.LMFactor:      " << itsSolveOpt.lmFactor << endl;
      os << "  Solve.Options.BalancedEqs:   " << itsSolveOpt.balancedEq << endl;
      os << "  Solve.Options.UseSVD:        " << itsSolveOpt.useSVD <<endl;
    }

    void Demixer::showTimings (std::ostream& os, double duration) const
    {
      const double self = itsTimer.getElapsed();

      os << "  ";
      FlagCounter::showPerc1 (os, self, duration);
      os << " Demixer " << itsName << endl;

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerPhaseShift.getElapsed(), self);
      os << " of it spent in phase shifting/averaging data" << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerDemix.getElapsed(), self);
      os << " of it spent in calculating decorrelation factors" << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerSolve.getElapsed(), self);
      os << " of it spent in solving source gains" << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerSubtract.getElapsed(), self);
      os << " of it spent in subtracting modeled sources" << endl;
    }

    bool Demixer::process (const DPBuffer& buf)
    {
      itsTimer.start();
      // Set start time of the chunk.
      if (itsNTimeIn == 0) {
        itsStartTimeChunk = buf.getTime() - itsInput->timeInterval() * 0.5;
      }
      // Update the count.
      itsNTimeIn++;
      // Make sure all required data arrays are filled in.
      DPBuffer newBuf(buf);
      RefRows refRows(newBuf.getRowNrs());
      if (newBuf.getUVW().empty()) {
        newBuf.setUVW(itsInput->fetchUVW(newBuf, refRows, itsTimer));
      }
      if (newBuf.getWeights().empty()) {
        newBuf.setWeights(itsInput->fetchWeights(newBuf, refRows, itsTimer));
      }
      if (newBuf.getFullResFlags().empty()) {
        newBuf.setFullResFlags(itsInput->fetchFullResFlags(newBuf, refRows,
                                                           itsTimer));
      }

      // Do the initial steps (phaseshift and average).
      itsTimerPhaseShift.start();
///#pragma omp parallel for
      for (int i=0; i<int(itsFirstSteps.size()); ++i) {
        itsFirstSteps[i]->process(newBuf);
      }
      itsAvgSubtr->process(newBuf);
      itsTimerPhaseShift.stop();

      // For each itsNTimeAvg times, calculate the
      // phase rotation per direction.
      itsTimerDemix.start();
      addFactors (newBuf, itsFactorBuf);
      if (itsNTimeIn % itsNTimeAvg == 0) {
        makeFactors (itsFactorBuf, itsFactors[itsNTimeOut],
                     itsAvgResults[0]->get()[itsNTimeOut].getWeights(),
                     itsNChanOut);
        // Deproject sources without a model.
        deproject (itsFactors[itsNTimeOut], itsAvgResults, itsNTimeOut);
        itsFactorBuf = Complex();   // clear summation buffer
        itsNTimeOut++;
      }
      // Subtract is done with different averaging parameters, so
      // calculate the factors for it.
      addFactors (newBuf, itsFactorBufSubtr);
      if (itsNTimeIn % itsNTimeAvgSubtr == 0) {
        makeFactors (itsFactorBufSubtr, itsFactorsSubtr[itsNTimeOutSubtr],
                     itsAvgResultSubtr->get()[itsNTimeOutSubtr].getWeights(),
                     itsNChanOutSubtr);
        itsFactorBufSubtr = Complex();   // clear summation buffer
        itsNTimeOutSubtr++;
      }
      itsTimerDemix.stop();

      // Do BBS solve, etc. when sufficient time slots have been collected.
      if (itsNTimeOut == itsNTimeChunk) {
        demix();
        itsNTimeIn       = 0;
        itsNTimeOut      = 0;
        itsNTimeOutSubtr = 0;
      }

      itsTimer.stop();
      return true;
    }

    void Demixer::finish()
    {
      // Process remaining entries.
      if (itsNTimeIn > 0) {
        itsTimer.start();

        // Finish the initial steps (phaseshift and average).
        itsTimerPhaseShift.start();
        ///#pragma omp parallel for
        for (int i=0; i<int(itsFirstSteps.size()); ++i) {
          itsFirstSteps[i]->finish();
        }
        itsAvgSubtr->finish();
        itsTimerPhaseShift.stop();
        // Only average if there is some data.
        itsTimerDemix.start();
        if (itsNTimeIn % itsNTimeAvg != 0) {
          makeFactors (itsFactorBuf, itsFactors[itsNTimeOut],
                       itsAvgResults[0]->get()[itsNTimeOut].getWeights(),
                       itsNChanOut);
          // Deproject sources without a model.
          deproject (itsFactors[itsNTimeOut], itsAvgResults, itsNTimeOut);
          itsNTimeOut++;
        }
        if (itsNTimeIn % itsNTimeAvgSubtr != 0) {
          makeFactors (itsFactorBufSubtr, itsFactorsSubtr[itsNTimeOutSubtr],
                       itsAvgResultSubtr->get()[itsNTimeOutSubtr].getWeights(),
                       itsNChanOutSubtr);
          itsNTimeOutSubtr++;
        }
        itsTimerDemix.stop();
        // Resize lists of mixing factors to the number of valid entries.
        itsFactors.resize(itsNTimeOut);
        itsFactorsSubtr.resize(itsNTimeOutSubtr);
        // Demix the source directions.
        demix();
        itsTimer.stop();
      }
      // Let the next steps finish.
      getNextStep()->finish();
    }

    void Demixer::addFactors (const DPBuffer& newBuf,
                              Array<DComplex>& factorBuf)
    {
      // Nothing to do if only target direction.
      if (itsNrDir <= 1) return;
///#pragma omp parallel
      {
///#pragma omp for
        uint ncorr  = newBuf.getData().shape()[0];
        uint nchan  = newBuf.getData().shape()[1];
        uint nbl    = newBuf.getData().shape()[2];
        DComplex* factorPtr = factorBuf.data();
        //# If ever in the future a time dependent phase center is used,
        //# the machine must be reset for each new time, thus each new call
        //# to process.
        for (uint i1=0; i1<itsNrDir-1; ++i1) {
          for (uint i0=i1+1; i0<itsNrDir; ++i0) {
            const double* uvw       = newBuf.getUVW().data();
            const bool*   flagPtr   = newBuf.getFlags().data();
            const float*  weightPtr = newBuf.getWeights().data();
            const DComplex* phasor1 = itsPhaseShifts[i1]->getPhasors().data();
            if (i0 == itsNrDir-1) {
              for (uint i=0; i<nbl; ++i) {
                for (uint j=0; j<nchan; ++j) {
                  DComplex factor = conj(*phasor1++);
                  for (uint k=0; k<ncorr; ++k) {
                    if (! *flagPtr) {
                      *factorPtr += factor * double(*weightPtr);
                    }
                    flagPtr++;
                    weightPtr++;
                    factorPtr++;
                  }
                }
                uvw += 3;
              }
            } else {
              const DComplex* phasor0 = itsPhaseShifts[i0]->getPhasors().data();
              for (uint i=0; i<nbl; ++i) {
                for (uint j=0; j<nchan; ++j) {
                  // Probably multiply with conj
                  DComplex factor = *phasor0++ / *phasor1++;
                  for (uint k=0; k<ncorr; ++k) {
                    if (! *flagPtr) {
                      *factorPtr += factor * double(*weightPtr);
                    }
                    flagPtr++;
                    weightPtr++;
                    factorPtr++;
                  }
                }
                uvw += 3;
              }
            }
          }
        }
      }
    }

    void Demixer::makeFactors (const Array<DComplex>& bufIn,
                               Array<DComplex>& bufOut,
                               const Cube<float>& weightSums,
                               uint nchanOut)
    {
      // Nothing to do if only target direction.
      if (itsNrDir <= 1) return;
      ASSERT (! weightSums.empty());
      bufOut.resize (IPosition(5, itsNrDir, itsNrDir,
                               itsNrCorr, nchanOut, itsNrBl));
      bufOut = DComplex(1,0);
      const DComplex* phin = bufIn.data();
      for (uint d0=0; d0<itsNrDir; ++d0) {
        for (uint d1=d0+1; d1<itsNrDir; ++d1) {
          DComplex* ph1 = bufOut.data() + d0*itsNrDir + d1;
          DComplex* ph2 = bufOut.data() + d1*itsNrDir + d0;
          // Average for all channels and divide by the summed weights.
          const float* weightPtr = weightSums.data();
          for (uint k=0; k<itsNrBl; ++k) {
            for (uint c0=0; c0<nchanOut; ++c0) {
              DComplex sum[4];
              uint nch = std::min(itsNChanAvg, itsNChanIn-c0*itsNChanAvg);
              for (uint c1=0; c1<nch; ++c1) {
                for (uint j=0; j<itsNrCorr; ++j) {
                  sum[j] += *phin++;
                }
              }
              for (uint j=0; j<itsNrCorr; ++j) {
                *ph1 = sum[j] / double(*weightPtr++);
                *ph2 = conj(*ph1);
                ph1 += itsNrDir*itsNrDir;
                ph2 += itsNrDir*itsNrDir;
              }
            }
          }
        }
      }
    }

    void Demixer::deproject (Array<DComplex>& factors,
                             vector<MultiResultStep*> avgResults,
                             uint resultIndex)
    {
      // Nothing to do if only target direction or if all sources are modeled.
      if (itsNrDir <= 1  ||  itsNrDir == itsNrModel) return;
      // Get pointers to the data for the various directions.
      vector<Complex*> resultPtr(itsNrDir);
      for (uint j=0; j<itsNrDir; ++j) {
        resultPtr[j] = avgResults[j]->get()[resultIndex].getData().data();
      }
      // Sources without a model have to be deprojected.
      uint nrDeproject = itsNrDir - itsNrModel;
      // The projection matrix is given by
      //     P = I - A * inv(A.T.conj * A) * A.T.conj
      // where A is the last column of the demixing matrix M.
      // The BBS equations get:
      //     P * M' * v_predict = P * v_averaged
      // where M' is obtained by removing the last column of demixing matrix M.
      // The dimensions of the matrices/vectors are:
      //     P : NxN
      //     M' : Nx(N-1)
      //     v_predict : (N-1) x 1
      //     v_averaged: N x 1
      // where N is the number of modeled sources to use in demixing.
      // In the general case S sources might not have a source model.
      // In that case A is the NxS matrix containing all these columns
      // from M and M' is the Nx(N-S) matrix without all these columns.

      // Calculate P for all baselines,channels,correlations.
      IPosition shape = factors.shape();
      int nvis = shape[2] * shape[3] * shape[4];
      shape[1] = itsNrModel;
      Array<DComplex> newFactors (shape);
      IPosition inShape (2, itsNrDir, itsNrDir);
      IPosition outShape(2, itsNrDir, itsNrModel);
      // omp parallel
      casa::Matrix<DComplex> a(itsNrDir, nrDeproject);
      casa::Matrix<DComplex> ma(itsNrDir, itsNrModel);
      vector<DComplex> vec(itsNrDir);
      // omp for
      for (int i=0; i<nvis; ++i) {
        // Split the matrix into the modeled and deprojected sources.
        // Copy the columns to the individual matrices.
        const DComplex* inptr  = factors.data() + i*itsNrDir*itsNrDir;
        DComplex* outptr = newFactors.data() + i*itsNrDir*itsNrModel;
        casa::Matrix<DComplex> out (outShape, outptr, SHARE);
        // Copying a bit of data is probably faster than taking a matrix subset.
        objcopy (ma.data(), inptr, itsNrDir*itsNrModel);
        objcopy (a.data(), inptr + itsNrDir*itsNrModel, itsNrDir*nrDeproject);
        // Calculate conjugated transpose of A, multiply with A, and invert.
        casa::Matrix<DComplex> at(adjoint(a));
        casa::Matrix<DComplex> ata(invert(product(at, a)));
        if (ata.empty()) {
          ata.resize (nrDeproject, nrDeproject);
        }
        DBGASSERT(ata.ncolumn()==nrDeproject && ata.nrow()==nrDeproject);
        // Calculate P = I - A * ata * A.T.conj
        casa::Matrix<DComplex> aata(product(a,ata));
        casa::Matrix<DComplex> p (-product(product(a, ata), at));
        casa::Vector<DComplex> diag(p.diagonal());
        diag += DComplex(1,0);
        // Multiply the demixing factors with P (get stored in newFactors).
        out = product(p, ma);
        ///        cout << "p matrix: " << p;
        // Multiply the averaged data point with P.
        std::fill (vec.begin(), vec.end(), DComplex());
        for (uint j=0; j<itsNrDir; ++j) {
          for (uint k=0; k<itsNrDir; ++k) {
            vec[k] += DComplex(resultPtr[j][i]) * p(k,j);
          }
        }
        // Put result back in averaged data for those sources.
        for (uint j=0; j<itsNrDir; ++j) {
          resultPtr[j][i] = vec[j];
        }
        ///        cout << vec << endl;
      }
      // Set the new demixing factors.
      factors.reference (newFactors);
    }

    void Demixer::demix()
    {
      size_t targetIndex = itsAvgResults.size() - 1;

      // Only solve and subtract if multiple directions.
      if (itsNrDir > 1) {
        // Collect buffers for each direction.
        vector<vector<DPBuffer> > buffers;
        for (size_t i=0; i<itsAvgResults.size(); ++i) {
          buffers.push_back (itsAvgResults[i]->get());
          // Do not clear target buffer, because it is shared with
          // itsAvgResultSubtr if averaging of demix and subtract is the same.
          if (i != targetIndex) {
            itsAvgResults[i]->clear();
          }
        }

        itsTimerSolve.start();
        // Construct grids for parameter estimation.
        Axis::ShPtr timeAxis (new RegularAxis (itsStartTimeChunk,
                                               itsTimeIntervalAvg,
                                               itsNTimeOut));
        Grid visGrid(itsFreqAxisDemix, timeAxis);
        // Solve for each time slot over all channels.
        Grid solGrid(itsFreqAxisDemix->compress(itsFreqAxisDemix->size()),
                     timeAxis);

        // Estimate model parameters.
        itsBBSExpr.estimate(buffers, visGrid, solGrid, itsFactors);
        itsTimerSolve.stop();

        itsTimerSubtract.start();
        // Construct subtract grid if it differs from the grid used for
        // parameter estimation.
        Axis::ShPtr timeAxisSubtr (new RegularAxis (itsStartTimeChunk,
                                                    itsTimeIntervalSubtr,
                                                    itsNTimeOutSubtr));
        visGrid = Grid(itsFreqAxisSubtr, timeAxisSubtr);
        // Subtract the sources.
        itsBBSExpr.subtract (itsAvgResultSubtr->get(), visGrid, itsFactorsSubtr,
                             targetIndex, itsSubtrSources.size());
        itsTimerSubtract.stop();
      }

      // Let the next step process the data.
      itsTimer.stop();
      for (uint i=0; i<itsNTimeOutSubtr; ++i) {
        getNextStep()->process (itsAvgResultSubtr->get()[i]);
        itsAvgResultSubtr->get()[i].clear();
        itsAvgResults[targetIndex]->get()[i].clear();
      }
      itsAvgResultSubtr->get().clear();
      itsAvgResults[targetIndex]->get().clear();
      itsTimer.start();
    }

    Axis::ShPtr Demixer::makeFreqAxis (uint nchanAvg)
    {
      casa::Vector<double> freq = itsInput->chanFreqs(nchanAvg);
      casa::Vector<double> width = itsInput->chanWidths(nchanAvg);
      ASSERT (allEQ (width, width[0]));

      return Axis::ShPtr(new RegularAxis (freq[0] - width[0] * 0.5, width[0],
        freq.size()));
    }

    string Demixer::toString (double value) const
    {
      ostringstream os;
      os << setprecision(16) << value;
      return os.str();
    }

    double Demixer::getAngle (const String& value) const
    {
      double angle;
      Quantity q;
      ASSERTSTR (Quantity::read (q, value),
                 "Demixer: " + value + " is not a proper angle");
      if (q.getUnit().empty()) {
        angle = q.getValue() / 180. * C::pi;
      } else {
        ASSERTSTR (q.getFullUnit().getValue() == UnitVal::ANGLE,
                   "Demixer: " + value + " is not a proper angle");
        angle = q.getValue("rad");
      }
      return angle;
    }

  } //# end namespace
}
