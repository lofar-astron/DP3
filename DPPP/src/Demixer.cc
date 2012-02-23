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
        itsBBSExpr       (*input, itsSkyName, itsInstrumentName),
        itsTarget        (parset.getString(prefix+"target", "target")),
        itsSubtrSources  (parset.getStringVector (prefix+"subtractsources")),
        itsModelSources  (parset.getStringVector (prefix+"modelsources",
                                                  vector<string>())),
        itsExtraSources  (parset.getStringVector (prefix+"othersources",
                                                  vector<string>())),
        itsJointSolve    (parset.getBool  (prefix+"jointsolve", true)),
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
      // See if averaging in demix and subtract differs.
      itsCalcSubtr = (itsNChanAvg != itsNChanAvgSubtr  ||
                      itsNTimeAvg != itsNTimeAvgSubtr);
      // Collect all source names.
      itsNrModel = itsSubtrSources.size() + itsModelSources.size();
      itsNrDir   = itsNrModel + itsExtraSources.size() + 1;
      itsAllSources.reserve (itsNrDir);
      itsAllSources.insert (itsAllSources.end(),
                            itsModelSources.begin(), itsModelSources.end());
      itsAllSources.insert (itsAllSources.end(),
                            itsSubtrSources.begin(), itsSubtrSources.end());
      itsAllSources.insert (itsAllSources.end(),
                            itsExtraSources.begin(), itsExtraSources.end());
//      itsAllSources.push_back("target"); //// probably not needed
      itsAllSources.push_back(itsTarget);

      // Size buffers.
      itsFactors.resize      (itsNTimeChunk);
      itsFactorsSubtr.resize (itsNTimeChunkSubtr);
      itsPhaseShifts.reserve (itsNrDir-1);
      itsFirstSteps.reserve  (itsNrDir);
      itsAvgResults.reserve  (itsNrDir);

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
        // Note the PhaseShift knows about source names CygA, etc.
        PhaseShift* step1 = new PhaseShift (input, parset,
                                            prefix + '.' + itsAllSources[i],
                                            itsAllSources[i]);
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

      // Do the same for the subtract if it has different averaging.
      if (itsCalcSubtr) {
        itsAvgSubtr = DPStep::ShPtr (new Averager(input, prefix,
                                                  itsNChanAvgSubtr,
                                                  itsNTimeAvgSubtr));
        itsAvgResultSubtr = new MultiResultStep(itsNTimeChunk);
        itsAvgSubtr->setNextStep (DPStep::ShPtr(itsAvgResultSubtr));
      } else {
        // Same averaging, so use demix averaging for the subtract.
        itsAvgSubtr = DPStep::ShPtr (new NullStep());
        itsAvgResultSubtr = targetAvgRes;
      }
      // Construct frequency axis for the demix and subtract averaging.
      itsFreqAxisDemix = makeFreqAxis (itsNChanAvg);
      itsFreqAxisSubtr = makeFreqAxis (itsNChanAvgSubtr);
    }

    Demixer::~Demixer()
    {
    }

    Axis::ShPtr Demixer::makeFreqAxis (uint nchanAvg)
    {
      casa::Vector<double> chanw = itsInput->chanWidths(nchanAvg);
      double chanWidth = chanw[0];
      ASSERT (allEQ (chanw, chanWidth));
      return Axis::ShPtr
        (new RegularAxis (itsInput->chanFreqs(nchanAvg)[0] - chanWidth*0.5,
                          chanWidth, chanw.size()));
    }

/*
    MDirection Demixer::handleCenter(const vector<string> &center) const
    {
      // A case-insensitive name can be given for a moving source (e.g. SUN)
      // or a known source (e.g. CygA).
      if (center.size() == 1) {
        return MDirection::makeMDirection (center[0]);
      }
      // The phase center must be given in J2000 as two values (ra,dec).
      // In the future time dependent frames can be done as in UVWFlagger.
      ASSERTSTR (center.size() == 2,
                 "2 values must be given in PhaseShift phasecenter");
      ///ASSERTSTR (center.size() < 4,
      ///"Up to 3 values can be given in UVWFlagger phasecenter");
      MDirection phaseCenter;
      if (center.size() == 1) {
        string str = toUpper(center[0]);
        MDirection::Types tp;
        ASSERTSTR (MDirection::getType(tp, str),
                   str << " is an invalid source type"
                   " in UVWFlagger phasecenter");
        return MDirection(tp);
      }
      Quantity q0, q1;
      ASSERTSTR (MVAngle::read (q0, center[0]),
                 center[0] << " is an invalid RA or longitude"
                 " in UVWFlagger phasecenter");
      ASSERTSTR (MVAngle::read (q1, center[1]),
                 center[1] << " is an invalid DEC or latitude"
                 " in UVWFlagger phasecenter");
      MDirection::Types type = MDirection::J2000;
      if (center.size() > 2) {
        string str = toUpper(center[2]);
        MDirection::Types tp;
        ASSERTSTR (MDirection::getType(tp, str),
                   str << " is an invalid direction type in UVWFlagger"
                   " in UVWFlagger phasecenter");
      }
      return MDirection(q0, q1, type);
    }
*/

    void Demixer::updateInfo (DPInfo& info)
    {
      info.setNeedVisData();
      info.setNeedWrite();
      itsNChanIn = info.nchan();
      itsNrBl    = info.nbaselines();
      itsNrCorr  = info.ncorr();
      itsFactorBuf.resize (IPosition(4, itsNrBl, itsNChanIn, itsNrCorr,
                                     itsNrDir*(itsNrDir-1)/2));
      itsTimeInterval = info.timeInterval();
      double refFreq = itsFreqAxisDemix->center(itsFreqAxisDemix->size() / 2);
      // Let the internal steps update their data.
      // Use a copy of the DPInfo, otherwise it is updated multiple times.
      DPInfo infocp;
      for (uint i=0; i<itsFirstSteps.size(); ++i) {
        infocp = info;
        DPStep::ShPtr step = itsFirstSteps[i];
        while (step) {
          step->updateInfo (infocp);
          step = step->getNextStep();
        }
        // Create the BBS model expression for sources with a model.
        if (i < itsNrModel) {
          itsBBSExpr.addModel (*itsInput, infocp, itsAllSources[i], refFreq);
        }
      }
      // Keep the averaged time interval.
      itsNChanOut = infocp.nchan();
      itsTimeIntervalAvg = infocp.timeInterval();
      // Update the info of this object.
      info.update (itsNChanAvgSubtr, itsNTimeAvgSubtr);
      itsNChanOutSubtr = info.nchan();
      itsTimeIntervalSubtr = info.timeInterval();
    }

    void Demixer::show (std::ostream& os) const
    {
      os << "Demixer " << itsName << std::endl;
      os << "  skymodel:       " << itsSkyName << std::endl;
      os << "  instrumentmodel:" << itsInstrumentName << std::endl;
      os << "  target:         " << itsTarget << std::endl;
      os << "  subtractsources:" << itsSubtrSources << std::endl;
      os << "  modelsources:   " << itsModelSources << std::endl;
      os << "  extrasources:   " << itsExtraSources << std::endl;
      os << "  jointsolve:     " << itsJointSolve << std::endl;
      os << "  freqstep:       " << itsNChanAvgSubtr << std::endl;
      os << "  timestep:       " << itsNTimeAvgSubtr << std::endl;
      os << "  demixfreqstep:  " << itsNChanAvg << std::endl;
      os << "  demixtimestep:  " << itsNTimeAvg << std::endl;
      os << "  timechunk:      " << itsNTimeChunk << std::endl;
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
      itsTimerPhaseShift.stop();

      // For each itsNTimeAvg times, calculate the
      // phase rotation per direction.
      itsTimerDemix.start();
      addFactors (newBuf, itsFactorBuf);
      if (itsNTimeIn % itsNTimeAvg == 0) {
        makeFactors (itsFactorBuf, itsFactors[itsNTimeOut],
                     itsAvgResults[0]->get()[itsNTimeOut].getWeights(),
                     itsNChanOut);
        // If needed, keep the original factors for subtraction.
        if (!itsCalcSubtr) {
          itsFactorsSubtr[itsNTimeOut].reference (itsFactors[itsNTimeOut].copy());
        }
        // Deproject sources without a model.
        deproject (itsFactors[itsNTimeOut], itsAvgResults, itsNTimeOut);
        itsFactorBuf = Complex();   // clear summation buffer
        itsNTimeOut++;
      }
      if (itsCalcSubtr) {
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
        itsTimerPhaseShift.stop();

        itsTimerDemix.start();
        makeFactors (itsFactorBuf, itsFactors[itsNTimeOut],
                     itsAvgResults[0]->get()[itsNTimeOut].getWeights(),
                     itsNChanOut);
        // Deproject sources without a model.
        // If needed, keep the original factors for subtraction.
        if (!itsCalcSubtr) {
          itsFactorsSubtr[itsNTimeOut].reference (itsFactors[itsNTimeOut].copy());
        }
        deproject (itsFactors[itsNTimeOut], itsAvgResults, itsNTimeOut);
        itsNTimeOut++;
        if (itsCalcSubtr) {
          makeFactors (itsFactorBufSubtr, itsFactorsSubtr[itsNTimeOutSubtr],
                       itsAvgResultSubtr->get()[itsNTimeOutSubtr].getWeights(),
                       itsNChanOutSubtr);
          itsNTimeOutSubtr++;
        }
        itsTimerDemix.stop();

        demix();
        itsTimer.stop();
      }
      // Let the next steps finish.
      getNextStep()->finish();
    }

    void Demixer::addFactors (const DPBuffer& newBuf,
                              Array<DComplex>& factorBuf)
    {
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
      itsTimerSolve.start();
      // Collect buffers for each direction.
      vector<vector<DPBuffer> > buffers;
      size_t targetIndex = itsAvgResults.size() - 1;
      for (size_t i=0; i<itsAvgResults.size(); ++i) {
        buffers.push_back (itsAvgResults[i]->get());
        // Do not clear target buffer, because it is shared with
        // itsAvgResultSubtr if averaging of demix and subtract is the same.
        if (i != targetIndex) {
          itsAvgResults[i]->clear();
        }
      }
      // Solve for the gains in the various directions.
      itsBBSExpr.setSolvables();
      // Make time axis and grid.
      Axis::ShPtr timeAxis (new RegularAxis (itsStartTimeChunk,
                                             itsTimeIntervalAvg,
                                             itsNTimeOut));
      Grid visGrid(itsFreqAxisDemix, timeAxis);
      // Solve for each time slot over all channels.
      Grid solGrid(itsFreqAxisDemix->compress(itsFreqAxisDemix->size()),
                   timeAxis);

      LOG_DEBUG_STR("SHAPES: " << itsFactors[0].shape() << " " << itsFreqAxisDemix->size() << " " << buffers[0][0].getData().shape());

      // Estimate model parameters.
      itsBBSExpr.estimate(buffers, visGrid, solGrid, itsFactors);
      itsTimerSolve.stop();
      // Subtract the modeled sources.
      itsTimerSubtract.start();
      LOG_DEBUG_STR("subtracting....");
      itsBBSExpr.subtract (itsAvgResultSubtr->get(), visGrid, itsFactorsSubtr,
                           itsNrDir-1, itsSubtrSources.size());
      itsTimerSubtract.stop();
      // Let the next step process the data.
      itsTimer.stop();
      for (uint i=0; i<itsNTimeOutSubtr; ++i) {
        getNextStep()->process (itsAvgResultSubtr->get()[i]);
        itsAvgResultSubtr->get()[i].clear();
        itsAvgResults[targetIndex]->get()[i].clear();
      }
      itsTimer.start();
    }

  } //# end namespace
}
