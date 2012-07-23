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
#include <DPPP/Apply.h>
#include <DPPP/Averager.h>
#include <DPPP/CursorUtilCasa.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/EstimateMixed.h>
#include <DPPP/PhaseShift.h>
#include <DPPP/ParSet.h>
#include <DPPP/Simulate.h>
#include <DPPP/SourceDBUtil.h>
#include <DPPP/SubtractMixed.h>

#include <DPPP/PointSource.h>

#include <ParmDB/Axis.h>
#include <ParmDB/SourceDB.h>
#include <ParmDB/ParmDB.h>
#include <ParmDB/ParmSet.h>
#include <ParmDB/ParmCache.h>
#include <ParmDB/Parm.h>

#include <Common/LofarLogger.h>
#include <Common/OpenMP.h>
#include <Common/StreamUtil.h>
#include <Common/lofar_iomanip.h>
#include <Common/lofar_iostream.h>
#include <Common/lofar_fstream.h>

#include <casa/Quanta/MVAngle.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <scimath/Mathematics/MatrixMathLA.h>

#include <boost/multi_array.hpp>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    using LOFAR::operator<<;

//    namespace {
//      void pack(const vector<size_t> &directions, size_t nBl, size_t nCh,
//        size_t nCr, const_cursor<dcomplex> in, cursor<dcomplex> out);
//      void pack(const vector<size_t> &directions, size_t nStation,
//        const_cursor<double> in, double *out);
//      void unpack(const vector<size_t> &directions, size_t nStation,
//        const double *in, cursor<double> out);
//      Patch::Ptr makePatch(const string &name);
//    } //# end unnamed namespace

    Demixer::Demixer (DPInput* input,
                      const ParSet& parset, const string& prefix)
      : itsInput         (input),
        itsName          (prefix),
        itsSkyName       (parset.getString(prefix+"skymodel", "sky")),
        itsInstrumentName(parset.getString(prefix+"instrumentmodel",
                                           "instrument")),
        itsTargetSource  (parset.getString(prefix+"targetsource", string())),
        itsSubtrSources  (parset.getStringVector (prefix+"subtractsources")),
        itsModelSources  (parset.getStringVector (prefix+"modelsources",
                                                  vector<string>())),
        itsExtraSources  (parset.getStringVector (prefix+"othersources",
                                                  vector<string>())),
//        itsCutOffs       (parset.getUintVector (prefix+"elevationcutoffs",
//                                                  vector<uint>())),
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
        itsNTimeOut      (0),
        itsNConverged    (0),
        itsTimeCount     (0)
    {
      // Get and set solver options.
//      itsSolveOpt.maxIter =
//        parset.getUint  (prefix+"Solve.Options.MaxIter", 300);
//      itsSolveOpt.epsValue =
//        parset.getDouble(prefix+"Solve.Options.EpsValue", 1e-9);
//      itsSolveOpt.epsDerivative =
//        parset.getDouble(prefix+"Solve.Options.EpsDerivative", 1e-9);
//      itsSolveOpt.colFactor =
//        parset.getDouble(prefix+"Solve.Options.ColFactor", 1e-9);
//      itsSolveOpt.lmFactor  =
//        parset.getDouble(prefix+"Solve.Options.LMFactor", 1.0);
//      itsSolveOpt.balancedEq =
//        parset.getBool  (prefix+"Solve.Options.BalancedEqs", false);
//      itsSolveOpt.useSVD  =
//        parset.getBool  (prefix+"Solve.Options.UseSVD", true);

      // Maybe optionally a parset parameter directions to give the
      // directions of unknown sources.
      // Or make sources a vector of vectors like [name, ra, dec] where
      // ra and dec are optional.

      ASSERTSTR (!(itsSkyName.empty() || itsInstrumentName.empty()),
                 "An empty name is given for the sky and/or instrument model");
      // Default nr of time chunks is maximum number of threads.
      if (itsNTimeChunk == 0) {
        itsNTimeChunk = OpenMP::maxThreads();
      }
      // Check that time windows fit nicely.
      ASSERTSTR ((itsNTimeChunk * itsNTimeAvg) % itsNTimeAvgSubtr == 0,
                 "time window should fit final averaging integrally");
      itsNTimeChunkSubtr = (itsNTimeChunk * itsNTimeAvg) / itsNTimeAvgSubtr;
      // Collect all source names.
      itsNModel = itsSubtrSources.size() + itsModelSources.size();
      itsNDir   = itsNModel + itsExtraSources.size() + 1;
      itsAllSources.reserve (itsNDir);
      itsAllSources.insert (itsAllSources.end(),
                            itsSubtrSources.begin(), itsSubtrSources.end());
      itsAllSources.insert (itsAllSources.end(),
                            itsModelSources.begin(), itsModelSources.end());
      itsAllSources.insert (itsAllSources.end(),
                            itsExtraSources.begin(), itsExtraSources.end());
      itsAllSources.push_back (itsTargetSource);

      // Get the source info of all patches from the SourceDB table.
      BBS::SourceDB sourceDB(BBS::ParmDBMeta("", itsSkyName), false);
      vector<string> patchNames(itsAllSources);
      // If the target source is given, add it to the model.
      // Because the target source has to be the last direction, it means
      // that (for the time being) no extra sources can be given.
      if (! itsTargetSource.empty()) {
        patchNames[itsNModel++] = itsTargetSource;
        // The target has to be the last demix direction.
        // If it has a source model, there cannot be any extra source
        // because the sources to be predicted have to be a consecutive vector.
        ASSERTSTR (itsExtraSources.empty(), "Currently no extrasources can "
                   "be given if the targetsource is given");
      }
      itsPatchList = makePatches (sourceDB, patchNames, itsNModel);
      ASSERT(itsPatchList.size() == itsNModel);

      // Size buffers.
      itsFactors.resize      (itsNTimeChunk);
      itsFactorsSubtr.resize (itsNTimeChunkSubtr);
      itsPhaseShifts.reserve (itsNDir-1);
      itsFirstSteps.reserve  (itsNDir+1);   // one extra for itsAvgSubtr
      itsAvgResults.reserve  (itsNDir);

      // Create the steps for the sources to be removed.
      // Demixing consists of the following steps:
      // - phaseshift data to each demix source
      // - average data in each direction, also for original phasecenter.
      // - determine demix factors for all directions
      // - use BBS to predict and solve in each direction. It is possible to
      //   predict more directions than to solve (for strong sources in field).
      // - use BBS to subtract the solved sources using the demix factors.
      //   The averaging used here can be smaller than used when solving.
      for (uint i=0; i<itsNDir-1; ++i) {
        // First make the phaseshift and average steps for each demix source.
        // The resultstep gets the result.
        // The phasecenter can be given in a parameter. Its name is the default.
        // Look up the source direction in the patch table.
        // If found, turn it into a vector of strings.
        vector<string> sourceVec (1, itsAllSources[i]);
        if(i < itsNModel) {
          sourceVec[0] = toString(itsPatchList[i]->position()[0]);
          sourceVec.push_back(toString(itsPatchList[i]->position()[1]));
        }
        PhaseShift* step1 = new PhaseShift (input, parset,
                                            prefix + itsAllSources[i] + '.',
                                            sourceVec);
        itsFirstSteps.push_back (DPStep::ShPtr(step1));
        itsPhaseShifts.push_back (step1);
        DPStep::ShPtr step2 (new Averager(input, prefix, itsNChanAvg,
                                          itsNTimeAvg));
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

      // Create the data average step for the subtract.
      DPStep::ShPtr targetAvgSubtr(new Averager(input, prefix,
                                                itsNChanAvgSubtr,
                                                itsNTimeAvgSubtr));
      itsAvgResultSubtr = new MultiResultStep(itsNTimeChunkSubtr);
      targetAvgSubtr->setNextStep (DPStep::ShPtr(itsAvgResultSubtr));
      itsFirstSteps.push_back (targetAvgSubtr);

//      while(itsCutOffs.size() < itsNModel) {
//        itsCutOffs.push_back(0);
//      }
//      itsCutOffs.resize(itsNModel);

//      itsFrames.reserve(OpenMP::maxThreads());
//      itsConverters.reserve(OpenMP::maxThreads());
//      for(size_t i = 0; i < OpenMP::maxThreads(); ++i) {
//        MeasFrame frame;
//        frame.set(input->arrayPos());
//        frame.set(MEpoch());
//        itsFrames.push_back(frame);

//        // Using AZEL because AZELGEO produces weird results.
//        itsConverters.push_back(MDirection::Convert(MDirection::Ref(MDirection::J2000),
//          MDirection::Ref(MDirection::AZEL, itsFrames.back())));
//      }
    }

    void Demixer::initUnknowns()
    {
      itsUnknowns.resize(IPosition(4, 8, itsNStation, itsNModel,
        itsNTimeDemix));
      itsUnknowns = 0.0;

//      itsErrors.resize(IPosition(4, 8, itsNStation, itsNModel,
//        itsNTimeDemix));
//      itsErrors = 0.0;

      itsLastKnowns.resize(IPosition(3, 8, itsNStation, itsNModel));
      for(uint dr = 0; dr < itsNModel; ++dr) {
        for(uint st = 0; st < itsNStation; ++st) {
          itsLastKnowns(IPosition(3, 0, st, dr)) = 1.0;
          itsLastKnowns(IPosition(3, 1, st, dr)) = 0.0;
          itsLastKnowns(IPosition(3, 2, st, dr)) = 0.0;
          itsLastKnowns(IPosition(3, 3, st, dr)) = 0.0;
          itsLastKnowns(IPosition(3, 4, st, dr)) = 0.0;
          itsLastKnowns(IPosition(3, 5, st, dr)) = 0.0;
          itsLastKnowns(IPosition(3, 6, st, dr)) = 1.0;
          itsLastKnowns(IPosition(3, 7, st, dr)) = 0.0;
        }
      }
    }

    Demixer::~Demixer()
    {
    }

    void Demixer::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      // Get size info.
      itsNStation = infoIn.antennaNames().size();
      itsNChanIn  = infoIn.nchan();
      itsNBl      = infoIn.nbaselines();
      itsNCorr    = infoIn.ncorr();
      ASSERTSTR (itsNCorr==4, "Demixing requires data with 4 polarizations");
      itsFactorBuf.resize (IPosition(4, itsNCorr, itsNChanIn, itsNBl,
                                     itsNDir*(itsNDir-1)/2));
      itsFactorBufSubtr.resize (IPosition(4, itsNCorr, itsNChanIn, itsNBl,
                                     itsNDir*(itsNDir-1)/2));

      // Adapt averaging to available nr of channels and times.
      // Use a copy of the DPInfo, otherwise it is updated multiple times.
      DPInfo infoDemix(infoIn);
      itsNTimeAvg = std::min (itsNTimeAvg, infoIn.ntime());
      itsNChanAvg = infoDemix.update (itsNChanAvg, itsNTimeAvg);
      itsNChanOut = infoDemix.nchan();
      itsTimeIntervalAvg = infoDemix.timeInterval();
      itsNTimeDemix      = infoDemix.ntime();
      for (size_t i=0; i<infoIn.getAnt1().size(); ++i) {
        itsBaselines.push_back (Baseline(infoIn.getAnt1()[i],
                                         infoIn.getAnt2()[i]));
      }

      // Let the internal steps update their data.
      for (uint i=0; i<itsFirstSteps.size(); ++i) {
        itsFirstSteps[i]->setInfo (infoIn);
      }
      // Update the info of this object.
      info().setNeedVisData();
      info().setNeedWrite();
      itsNTimeAvgSubtr = std::min (itsNTimeAvgSubtr, infoIn.ntime());
      itsNChanAvgSubtr = info().update (itsNChanAvgSubtr, itsNTimeAvgSubtr);
      itsNChanOutSubtr = info().nchan();
      ASSERTSTR (itsNChanAvg % itsNChanAvgSubtr == 0,
		 "Demix averaging " << itsNChanAvg
		 << " must be multiple of output averaging "
		 << itsNChanAvgSubtr);
      ASSERTSTR (itsNTimeAvg % itsNTimeAvgSubtr == 0,
		 "Demix averaging " << itsNTimeAvg
		 << " must be multiple of output averaging "
		 << itsNTimeAvgSubtr);
      // Store channel frequencies for the demix and subtract resolutions.
      itsFreqDemix = infoDemix.chanFreqs();
      itsFreqSubtr = getInfo().chanFreqs();

      // Store phase center position in J2000.
      MDirection dirJ2000(MDirection::Convert(infoIn.phaseCenter(),
                                              MDirection::J2000)());
      Quantum<Vector<Double> > angles = dirJ2000.getAngle();
      itsPhaseRef = Position(angles.getBaseValue()[0],
                             angles.getBaseValue()[1]);

      initUnknowns();
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
//      os << "  jointsolve:     " << itsJointSolve << std::endl;
      os << "  freqstep:       " << itsNChanAvgSubtr << std::endl;
      os << "  timestep:       " << itsNTimeAvgSubtr << std::endl;
      os << "  demixfreqstep:  " << itsNChanAvg << std::endl;
      os << "  demixtimestep:  " << itsNTimeAvg << std::endl;
      os << "  ntimechunk:     " << itsNTimeChunk << std::endl;
//      os << "  Solve.Options.MaxIter:       " << itsSolveOpt.maxIter << endl;
//      os << "  Solve.Options.EpsValue:      " << itsSolveOpt.epsValue << endl;
//      os << "  Solve.Options.EpsDerivative: " << itsSolveOpt.epsDerivative << endl;
//      os << "  Solve.Options.ColFactor:     " << itsSolveOpt.colFactor << endl;
//      os << "  Solve.Options.LMFactor:      " << itsSolveOpt.lmFactor << endl;
//      os << "  Solve.Options.BalancedEqs:   " << itsSolveOpt.balancedEq << endl;
//      os << "  Solve.Options.UseSVD:        " << itsSolveOpt.useSVD <<endl;
    }

    void Demixer::showCounts (std::ostream& os) const
    {
      os << endl << "Statistics for Demixer " << itsName;
      os << endl << "======================" << endl;
      os << endl << "Converged: " << itsNConverged << "/" << itsNTimeDemix
        << " cells" << endl;
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
      os << " of it spent in estimating gains and computing residuals" << endl;
    }

    bool Demixer::process (const DPBuffer& buf)
    {
      itsTimer.start();
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
                     itsNChanOut,
                     itsNChanAvg);
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
                     itsNChanOutSubtr,
                     itsNChanAvgSubtr);
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
        itsTimeCount += itsNTimeChunk;
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
        for (int i=0; i<int(itsFirstSteps.size()); ++i) {
          itsFirstSteps[i]->finish();
        }
        itsTimerPhaseShift.stop();
        // Only average if there is some unaveraged data.
        itsTimerDemix.start();
        if (itsNTimeIn % itsNTimeAvg != 0) {
          makeFactors (itsFactorBuf, itsFactors[itsNTimeOut],
                       itsAvgResults[0]->get()[itsNTimeOut].getWeights(),
                       itsNChanOut,
                       itsNChanAvg);
          // Deproject sources without a model.
          deproject (itsFactors[itsNTimeOut], itsAvgResults, itsNTimeOut);
          itsNTimeOut++;
        }
        if (itsNTimeIn % itsNTimeAvgSubtr != 0) {
          makeFactors (itsFactorBufSubtr, itsFactorsSubtr[itsNTimeOutSubtr],
                       itsAvgResultSubtr->get()[itsNTimeOutSubtr].getWeights(),
                       itsNChanOutSubtr,
                       itsNChanAvgSubtr);
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

      // Write solutions to disk in ParmDB format.
      dumpSolutions();

      // Let the next steps finish.
      getNextStep()->finish();
    }

    void Demixer::addFactors (const DPBuffer& newBuf,
                              Array<DComplex>& factorBuf)
    {
      // Nothing to do if only target direction.
      if (itsNDir <= 1) return;
      int ncorr  = newBuf.getData().shape()[0];
      int nchan  = newBuf.getData().shape()[1];
      int nbl    = newBuf.getData().shape()[2];
      int ncc    = ncorr*nchan;
      //# If ever in the future a time dependent phase center is used,
      //# the machine must be reset for each new time, thus each new call
      //# to process.
      // Add the weighted factors for each pair of directions.
      // The input factor is the phaseshift from target direction to
      // source direction. By combining them you get the shift from one
      // source direction to another.
      int dirnr = 0;
      for (uint i1=0; i1<itsNDir-1; ++i1) {
        for (uint i0=i1+1; i0<itsNDir; ++i0) {
          if (i0 == itsNDir-1) {
            // The last direction is the target direction, so no need to
            // combine the factors. Take conj to get shift source to target.
#pragma omp parallel for
            for (int i=0; i<nbl; ++i) {
              const bool*   flagPtr   = newBuf.getFlags().data() + i*ncc;
              const float*  weightPtr = newBuf.getWeights().data() + i*ncc;
              DComplex* factorPtr     = factorBuf.data() + (dirnr*nbl + i)*ncc;
              const DComplex* phasor1 = itsPhaseShifts[i1]->getPhasors().data()
                                        + i*nchan;
              for (int j=0; j<nchan; ++j) {
                DComplex factor = conj(*phasor1++);
                for (int k=0; k<ncorr; ++k) {
                  if (! *flagPtr) {
                    *factorPtr += factor * double(*weightPtr);
                  }
                  flagPtr++;
                  weightPtr++;
                  factorPtr++;
                }
              }
            } // end omp parallel for
          } else {
            // Different source directions; take both phase terms into account.
#pragma omp parallel for
            for (int i=0; i<nbl; ++i) {
              const bool*   flagPtr   = newBuf.getFlags().data() + i*ncc;
              const float*  weightPtr = newBuf.getWeights().data() + i*ncc;
              DComplex* factorPtr     = factorBuf.data() + (dirnr*nbl + i)*ncc;
              const DComplex* phasor0 = itsPhaseShifts[i0]->getPhasors().data()
                                        + i*nchan;
              const DComplex* phasor1 = itsPhaseShifts[i1]->getPhasors().data()
                                        + i*nchan;
              for (int j=0; j<nchan; ++j) {
                DComplex factor = *phasor0++ * conj(*phasor1++);
                for (int k=0; k<ncorr; ++k) {
                  if (! *flagPtr) {
                    *factorPtr += factor * double(*weightPtr);
                  }
                  flagPtr++;
                  weightPtr++;
                  factorPtr++;
                }
              }
            } // end omp parallel for
          }

          // Next direction pair.
          dirnr++;
        }
      }
    }

    void Demixer::makeFactors (const Array<DComplex>& bufIn,
                               Array<DComplex>& bufOut,
                               const Cube<float>& weightSums,
                               uint nChanOut,
                               uint nChanAvg)
    {
      // Nothing to do if only target direction.
      if (itsNDir <= 1) return;
      ASSERT (! weightSums.empty());
      bufOut.resize (IPosition(5, itsNDir, itsNDir,
                               itsNCorr, nChanOut, itsNBl));
      bufOut = DComplex(1,0);
      int ncc = itsNCorr*nChanOut;
      int nccdd = ncc*itsNDir*itsNDir;
      int nccin = itsNCorr*itsNChanIn;
      // Fill the factors for each combination of different directions.
      uint dirnr = 0;
      for (uint d0=0; d0<itsNDir; ++d0) {
        for (uint d1=d0+1; d1<itsNDir; ++d1) {
#pragma omp parallel for
          // Average factors by summing channels.
          // Note that summing in time is done in addFactors.
          // The sum per output channel is divided by the summed weight.
          // Note there is a summed weight per baseline,outchan,corr.
          for (int k=0; k<int(itsNBl); ++k) {
            const DComplex* phin = bufIn.data() + (dirnr*itsNBl + k)*nccin;
            DComplex* ph1 = bufOut.data() + k*nccdd + (d0*itsNDir + d1);
            DComplex* ph2 = bufOut.data() + k*nccdd + (d1*itsNDir + d0);
            const float* weightPtr = weightSums.data() + k*ncc;
            for (uint c0=0; c0<nChanOut; ++c0) {
              // Sum the factors for the input channels to average.
              DComplex sum[4];
              // In theory the last output channel could consist of fewer
              // input channels, so take care of that.
              uint nch = std::min(nChanAvg, itsNChanIn-c0*nChanAvg);
              for (uint c1=0; c1<nch; ++c1) {
                for (uint j=0; j<itsNCorr; ++j) {
                  sum[j] += *phin++;
                }
              }
              for (uint j=0; j<itsNCorr; ++j) {
                *ph1 = sum[j] / double(*weightPtr++);
                *ph2 = conj(*ph1);
                ph1 += itsNDir*itsNDir;
                ph2 += itsNDir*itsNDir;
              }
            }
          } // end omp parallel for
          // Next input direction pair.
          dirnr++;
        }
      }
    }

    void Demixer::deproject (Array<DComplex>& factors,
                             vector<MultiResultStep*> avgResults,
                             uint resultIndex)
    {
      // Nothing to do if only target direction or if all sources are modeled.
      if (itsNDir <= 1 || itsNDir == itsNModel) return;
      // Get pointers to the data for the various directions.
      vector<Complex*> resultPtr(itsNDir);
      for (uint j=0; j<itsNDir; ++j) {
        resultPtr[j] = avgResults[j]->get()[resultIndex].getData().data();
      }
      // Sources without a model have to be deprojected.
      uint nrDeproject = itsNDir - itsNModel;
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
      shape[1] = itsNModel;
      Array<DComplex> newFactors (shape);
      IPosition inShape (2, itsNDir, itsNDir);
      IPosition outShape(2, itsNDir, itsNModel);
///#pragma omp parallel
      {
        Matrix<DComplex> a(itsNDir, nrDeproject);
        Matrix<DComplex> ma(itsNDir, itsNModel);
        vector<DComplex> vec(itsNDir);
        ///#pragma omp for
        for (int i=0; i<nvis; ++i) {
          // Split the matrix into the modeled and deprojected sources.
          // Copy the columns to the individual matrices.
          const DComplex* inptr  = factors.data() + i*itsNDir*itsNDir;
          DComplex* outptr = newFactors.data() + i*itsNDir*itsNModel;
          Matrix<DComplex> out (outShape, outptr, SHARE);
          // Copying a bit of data is probably faster than taking a matrix
          // subset.
          objcopy (ma.data(), inptr, itsNDir*itsNModel);
          objcopy (a.data(), inptr + itsNDir*itsNModel, itsNDir*nrDeproject);
          // Calculate conjugated transpose of A, multiply with A, and invert.
          Matrix<DComplex> at(adjoint(a));
          Matrix<DComplex> ata(invert(product(at, a)));
          if (ata.empty()) {
            ata.resize (nrDeproject, nrDeproject);
          }
          DBGASSERT(ata.ncolumn()==nrDeproject && ata.nrow()==nrDeproject);
          // Calculate P = I - A * ata * A.T.conj
          Matrix<DComplex> aata(product(a,ata));
          Matrix<DComplex> p (-product(product(a, ata), at));
          Vector<DComplex> diag(p.diagonal());
          diag += DComplex(1,0);
          // Multiply the demixing factors with P (get stored in newFactors).
          out = product(p, ma);
          // Multiply the averaged data point with P.
          std::fill (vec.begin(), vec.end(), DComplex());
          for (uint j=0; j<itsNDir; ++j) {
            for (uint k=0; k<itsNDir; ++k) {
              vec[k] += DComplex(resultPtr[j][i]) * p(k,j);
            }
          }
          // Put result back in averaged data for those sources.
          for (uint j=0; j<itsNDir; ++j) {
            resultPtr[j][i] = vec[j];
          }
        }
      }
      // Set the new demixing factors.
      factors.reference (newFactors);
    }

    void Demixer::demix()
    {
      // Only solve and subtract if multiple directions.
      if (itsNDir > 1) {
        // Collect buffers for each direction.
        vector<vector<DPBuffer> > streams;
        for (size_t i=0; i<itsAvgResults.size(); ++i) {
          // Only the phased shifted and averaged data for directions which have
          // an associated model should be used to estimate the directional
          // response.
          if(i < itsNModel) {
            streams.push_back (itsAvgResults[i]->get());
          }
          // Clear the buffers.
          itsAvgResults[i]->clear();
        }

        itsTimerSolve.start();

        const size_t nTime = streams[0].size();
        const size_t nThread = OpenMP::maxThreads();
        const size_t nDr = itsNModel;
        const size_t nSt = itsNStation;
        const size_t nBl = itsBaselines.size();
        const size_t nCh = itsFreqDemix.size();
        const size_t nCr = 4;
        const size_t timeFactor = itsNTimeAvg / itsNTimeAvgSubtr;
        vector<DPBuffer> &target = itsAvgResultSubtr->get();

//        // Thread-private index per direction of the last known "good" solution.
//        boost::multi_array<int, 2> last(boost::extents[nDr][nThread]);
//        fill(&(last[0][0]), &(last[0][0]) + nDr * nThread, -1);

        // Thread-private buffer for the unknowns.
        boost::multi_array<double, 2> unknowns(boost::extents[nThread]
            [nDr * nSt * 8]);
//        boost::multi_array<double, 2> errors(boost::extents[nThread]
//            [nDr * nSt * 8]);

        // Copy solutions from global solution array to thread private solution
        // array (solution propagation between chunks).
        for(size_t i = 0; i < nThread; ++i)
        {
            copy(&(itsLastKnowns(IPosition(3, 0, 0, 0))),
                &(itsLastKnowns(IPosition(3, 0, 0, 0))) + nDr * nSt * 8,
                &(unknowns[i][0]));
        }

        // Thread-private buffer for station UVW coordinates.
        boost::multi_array<double, 3> uvw(boost::extents[nThread][nSt][3]);

        // Thread-private buffer for the simulated visibilities for all
        // directions at the demix resolution.
        boost::multi_array<dcomplex, 5> buffer(boost::extents[nThread][nDr][nBl]
            [nCh][nCr]);

        // Thread-private buffer for the simulated visibilities for a single
        // direction at the residual resolution.
        const size_t nChRes = itsFreqSubtr.size();
        boost::multi_array<dcomplex, 4> residual(boost::extents[nThread][nBl]
            [nChRes][nCr]);

        // Thread-private convergence counter.
        vector<size_t> converged(nThread, 0);

#pragma omp parallel for
        for(size_t i = 0; i < nTime; ++i)
        {
            const size_t thread = OpenMP::threadNum();

            // Compute elevations and determine which sources are visible.
            // TODO: Investigate segfault in MeasRef...
//            vector<size_t> drUp;
//            Quantum<Double> qEpoch(streams[0][i].getTime(), "s");
//            MEpoch mEpoch(qEpoch, MEpoch::UTC);
//            itsFrames[thread].set(mEpoch);

//            for(size_t dr = 0; dr < nDr; ++dr)
//            {
//                MVDirection drJ2000(itsPatchList[dr]->position()[0],
//                    itsPatchList[dr]->position()[1]);
//                MVDirection mvAzel(itsConverters[thread](drJ2000).getValue());
//                Vector<Double> azel = mvAzel.getAngle("deg").getValue();

//                if(azel(1) >= itsCutOffs[dr])
//                {
//                    drUp.push_back(dr);
//                    last[dr][thread] = i;
//                }
//            }

//            const size_t nDrUp = drUp.size();
//            if(nDrUp == 0)
//            {
//                continue;
//            }

            // Simulate.
            size_t strides[3] = {1, nCr, nCh * nCr};
            size_t strides_split[2] = {1, 3};

            const_cursor<Baseline> cr_baseline(&(itsBaselines[0]));
            const_cursor<double> cr_freq(&(itsFreqDemix[0]));
            const_cursor<double> cr_freqRes(&(itsFreqSubtr[0]));

            cursor<double> cr_split(&(uvw[thread][0][0]), 2, strides_split);
//            for(size_t dr = 0; dr < nDrUp; ++dr)
            for(size_t dr = 0; dr < nDr; ++dr)
            {
                fill(&(buffer[thread][dr][0][0][0]),
                    &(buffer[thread][dr][0][0][0]) + nBl * nCh * nCr,
                    dcomplex(0.0, 0.0));

//                const_cursor<double> cr_uvw =
//                    casa_const_cursor(streams[drUp[dr]][i].getUVW());
                const_cursor<double> cr_uvw =
                    casa_const_cursor(streams[dr][i].getUVW());
                splitUVW(nSt, nBl, cr_baseline, cr_uvw, cr_split);

                cursor<dcomplex> cr_model(&(buffer[thread][dr][0][0][0]), 3,
                    strides);
//                simulate(itsPatchList[drUp[dr]]->position(),
//                    itsPatchList[drUp[dr]], nSt, nBl, nCh, cr_baseline, cr_freq,
//                    cr_split, cr_model);
                simulate(itsPatchList[dr]->position(), itsPatchList[dr], nSt,
                    nBl, nCh, cr_baseline, cr_freq, cr_split, cr_model);
            }

            // Estimate.
//            vector<const_cursor<fcomplex> > cr_data(nDrUp);
//            vector<const_cursor<dcomplex> > cr_model(nDrUp);
            vector<const_cursor<fcomplex> > cr_data(nDr);
            vector<const_cursor<dcomplex> > cr_model(nDr);
//            for(size_t dr = 0; dr < nDrUp; ++dr)
            for(size_t dr = 0; dr < nDr; ++dr)
            {
//                cr_data[dr] = casa_const_cursor(streams[drUp[dr]][i].getData());
                cr_data[dr] = casa_const_cursor(streams[dr][i].getData());
                cr_model[dr] =
                    const_cursor<dcomplex>(&(buffer[thread][dr][0][0][0]), 3,
                    strides);
            }

            const_cursor<bool> cr_flag =
                casa_const_cursor(streams[0][i].getFlags());
            const_cursor<float> cr_weight =
                casa_const_cursor(streams[0][i].getWeights());

            bool status = false;
//            if(nDrUp != nDr)
//            {
//                vector<double> unknowns_packed(nDrUp * nSt * 8);
//                vector<double> errors_packed(nDrUp * nSt * 8);
//                vector<dcomplex> mix_packed(nDrUp * nDrUp * nCr * nCh * nBl);

//                // pack unknowns, mixing factors
//                size_t strides_mix[5] = {1, nDrUp, nDrUp * nDrUp, nCr * nDrUp
//                  * nDrUp, nCh * nCr * nDrUp * nDrUp};
//                pack(drUp, nBl, nCh, nCr, casa_const_cursor(itsFactors[i]),
//                  cursor<dcomplex>(&(mix_packed[0]), 5, strides_mix));

//                size_t strides_unknowns[3] = {1, 8, nSt * 8};
//                pack(drUp, nSt, const_cursor<double>(&(unknowns[thread][0]), 3,
//                  strides_unknowns), &(unknowns_packed[0]));

//                const_cursor<dcomplex> cr_mix(&(mix_packed[0]), 5, strides_mix);

//                // estimate
//                status = estimate(nDrUp, nSt, nBl, nCh, cr_baseline, cr_data,
//                    cr_model, cr_flag, cr_weight, cr_mix, &(unknowns_packed[0]),
//                    &(errors_packed[0]));

//                // unpack unknowns, errors
//                unpack(drUp, nSt, &(unknowns_packed[0]),
//                  cursor<double>(&(unknowns[thread][0]), 3, strides_unknowns));
//                unpack(drUp, nSt, &(errors_packed[0]),
//                  cursor<double>(&(errors[thread][0]), 3, strides_unknowns));
//            }
//            else
            {
                const_cursor<dcomplex> cr_mix =
                    casa_const_cursor(itsFactors[i]);
//                status = estimate(nDr, nSt, nBl, nCh, cr_baseline, cr_data,
//                    cr_model, cr_flag, cr_weight, cr_mix,
//                    &(unknowns[thread][0]), &(errors[thread][0]));
                status = estimate(nDr, nSt, nBl, nCh, cr_baseline, cr_data,
                    cr_model, cr_flag, cr_weight, cr_mix,
                    &(unknowns[thread][0]), 0);
            }

            if(status)
            {
                ++converged[thread];
            }

            const size_t nTimeResidual = min(timeFactor, target.size() - i
                * timeFactor);

            const size_t nDrRes = itsSubtrSources.size();

            size_t strides_model[3] = {1, nCr, nChRes * nCr};
            cursor<dcomplex> cr_model_res(&(residual[thread][0][0][0]), 3,
                strides_model);

            for(size_t j = 0; j < nTimeResidual; ++j)
            {
//                for(size_t dr = 0; dr < nDrUp; ++dr)
                for(size_t dr = 0; dr < nDr; ++dr)
                {
//                    const size_t drIdx = drUp[dr];
//                    if(drIdx >= nDrRes)
                    if(dr >= nDrRes)
                    {
                        break;
                    }

                    // Re-simulate for residual if required.
                    if(timeFactor != 1 || nCh != nChRes)
                    {
                        fill(&(residual[thread][0][0][0]),
                            &(residual[thread][0][0][0]) + nBl * nChRes * nCr,
                            dcomplex());

                        const_cursor<double> cr_uvw = casa_const_cursor(target[i
                            * timeFactor + j].getUVW());
                        splitUVW(nSt, nBl, cr_baseline, cr_uvw, cr_split);

//                        rotateUVW(itsPhaseRef, itsPatchList[drIdx]->position(),
//                            nSt, cr_split);
                        rotateUVW(itsPhaseRef, itsPatchList[dr]->position(),
                            nSt, cr_split);

//                        simulate(itsPatchList[drIdx]->position(),
//                            itsPatchList[drIdx], nSt, nBl, nChRes, cr_baseline,
//                            cr_freqRes, cr_split, cr_model_res);
                        simulate(itsPatchList[dr]->position(), itsPatchList[dr],
                            nSt, nBl, nChRes, cr_baseline, cr_freqRes, cr_split,
                            cr_model_res);
                    }
                    else
                    {
                        copy(&(buffer[thread][dr][0][0][0]),
                            &(buffer[thread][dr][0][0][0]) + nBl * nChRes * nCr,
                            &(residual[thread][0][0][0]));
                    }

                    // Apply solutions.
                    size_t strides_coeff[2] = {1, 8};
//                    const_cursor<double> cr_coeff(&(unknowns[thread][drIdx * nSt
//                        * 8]), 2, strides_coeff);
                    const_cursor<double> cr_coeff(&(unknowns[thread][dr * nSt
                        * 8]), 2, strides_coeff);
                    apply(nBl, nChRes, cr_baseline, cr_coeff, cr_model_res);

                    // Subtract.
                    cursor<fcomplex> cr_residual =
                        casa_cursor(target[i * timeFactor + j].getData());

                    const IPosition tmp_strides_mix_res =
                        itsFactorsSubtr[i * timeFactor + j].steps();

                    // The target direction is always last, and therefore it has
                    // index itsNDir - 1. The directions to subtract are always
                    // first, and therefore have indices [0, nDrRes).
//                    const_cursor<dcomplex> cr_mix(&(itsFactorsSubtr[i
//                        * timeFactor + j](IPosition(5, itsNDir - 1, drIdx, 0, 0,
//                        0))), 3, tmp_strides_mix_res.storage() + 2);
                    const_cursor<dcomplex> cr_mix(&(itsFactorsSubtr[i
                        * timeFactor + j](IPosition(5, itsNDir - 1, dr, 0, 0,
                        0))), 3, tmp_strides_mix_res.storage() + 2);
                    subtract(nBl, nChRes, cr_baseline, cr_residual,
                        cr_model_res, cr_mix);
                }
            }

            // Copy solutions to global solution array.
//            for(size_t dr = 0; dr < nDrUp; ++dr)
            for(size_t dr = 0; dr < nDr; ++dr)
            {
//                size_t idx = drUp[dr];
//                copy(&(unknowns[thread][idx * nSt * 8]), &(unknowns[thread][idx
//                    * nSt * 8]) + nSt * 8, &(itsUnknowns(IPosition(4, 0, 0, idx,
//                    itsTimeCount + i))));
//                copy(&(errors[thread][idx * nSt * 8]), &(errors[thread][idx
//                    * nSt * 8]) + nSt * 8, &(itsErrors(IPosition(4, 0, 0, idx,
//                    itsTimeCount + i))));
                copy(&(unknowns[thread][dr * nSt * 8]), &(unknowns[thread][dr
                    * nSt * 8]) + nSt * 8, &(itsUnknowns(IPosition(4, 0, 0, dr,
                    itsTimeCount + i))));
            }
        }

        // Store last known solutions.
        for(size_t dr = 0; dr < nDr; ++dr)
        {
//            int idx = last[dr][0];
//            for(size_t i = 1; i < nThread; ++i)
//            {
//                idx = max(idx, last[dr][i]);
//            }
            int idx = (nTime > 0 ? (nTime - 1) : -1);

            if(idx >= 0)
            {
                copy(&(itsUnknowns(IPosition(4, 0, 0, dr, itsTimeCount + idx))),
                    &(itsUnknowns(IPosition(4, 0, 0, dr, itsTimeCount + idx)))
                    + nSt * 8, &(itsLastKnowns(IPosition(3, 0, 0, dr))));
            }
        }

        // Update convergence count.
        for (size_t i=0; i<nThread; ++i) {
          itsNConverged += converged[i];
        }

        itsTimerSolve.stop();
      }

      // Let the next step process the data.
      itsTimer.stop();
      for (uint i=0; i<itsNTimeOutSubtr; ++i) {
        getNextStep()->process (itsAvgResultSubtr->get()[i]);
        itsAvgResultSubtr->get()[i].clear();
      }
      // Clear the vector in the MultiStep.
      itsAvgResultSubtr->clear();
      itsTimer.start();
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

    void Demixer::dumpSolutions()
    {
      // Construct solution grid.
      const Vector<double>& freq      = getInfo().chanFreqs();
      const Vector<double>& freqWidth = getInfo().chanWidths();
      BBS::Axis::ShPtr freqAxis(new BBS::RegularAxis(freq[0] - freqWidth[0]
        * 0.5, freqWidth[0], 1));
      BBS::Axis::ShPtr timeAxis(new BBS::RegularAxis(getInfo().startTime()
        - getInfo().timeInterval() * 0.5, itsTimeIntervalAvg, itsNTimeDemix));
      BBS::Grid solGrid(freqAxis, timeAxis);

      // Create and initialize ParmDB.
      BBS::ParmDB parmDB(BBS::ParmDBMeta("casa", itsInstrumentName), true);
      BBS::ParmSet parmSet;
      BBS::ParmCache parmCache(parmSet, solGrid.getBoundingBox());

      // Store the (freq, time) resolution of the solutions.
      vector<double> resolution(2);
      resolution[0] = freqWidth[0];
      resolution[1] = itsTimeIntervalAvg;
      parmDB.setDefaultSteps(resolution);

      // Convert station names from casa::String to std::string.
      ASSERT(getInfo().antennaNames().size() == itsNStation);
      vector<string> stations(itsNStation);
      copy(getInfo().antennaNames().begin(), getInfo().antennaNames().end(),
          stations.begin());

      vector<BBS::Parm> parms;
      for(size_t dr = 0; dr < itsNModel; ++dr) {
        for(size_t st = 0; st < itsNStation; ++st) {
          parms.push_back(BBS::Parm(parmCache, parmSet.addParm(parmDB,
            "DirectionalGain:0:0:Real:" + stations[st] + ":"
            + itsAllSources[dr])));
          parms.push_back(BBS::Parm(parmCache, parmSet.addParm(parmDB,
            "DirectionalGain:0:0:Imag:" + stations[st] + ":"
            + itsAllSources[dr])));

          parms.push_back(BBS::Parm(parmCache, parmSet.addParm(parmDB,
            "DirectionalGain:0:1:Real:" + stations[st] + ":"
            + itsAllSources[dr])));
          parms.push_back(BBS::Parm(parmCache, parmSet.addParm(parmDB,
            "DirectionalGain:0:1:Imag:" + stations[st] + ":"
            + itsAllSources[dr])));

          parms.push_back(BBS::Parm(parmCache, parmSet.addParm(parmDB,
            "DirectionalGain:1:0:Real:" + stations[st] + ":"
            + itsAllSources[dr])));
          parms.push_back(BBS::Parm(parmCache, parmSet.addParm(parmDB,
            "DirectionalGain:1:0:Imag:" + stations[st] + ":"
            + itsAllSources[dr])));

          parms.push_back(BBS::Parm(parmCache, parmSet.addParm(parmDB,
            "DirectionalGain:1:1:Real:" + stations[st] + ":"
            + itsAllSources[dr])));
          parms.push_back(BBS::Parm(parmCache, parmSet.addParm(parmDB,
            "DirectionalGain:1:1:Imag:" + stations[st] + ":"
            + itsAllSources[dr])));
        }
      }

      // Cache parameter values.
      parmCache.cacheValues();

      // Assign solution grid to parameters.
      for(size_t i = 0; i < parms.size(); ++i) {
        parms[i].setSolveGrid(solGrid);
      }

      // Write solutions.
      for(size_t ts = 0; ts < itsNTimeDemix; ++ts) {
        double *unknowns = &(itsUnknowns(IPosition(4, 0, 0, 0, ts)));
//        double *errors = &(itsErrors(IPosition(4, 0, 0, 0, ts)));
        for(size_t i = 0; i < parms.size(); ++i) {
//          parms[i].setCoeff(BBS::Location(0, ts), unknowns + i, 1, errors + i);
          parms[i].setCoeff(BBS::Location(0, ts), unknowns + i, 1);
        }
      }

      // Flush solutions to disk.
      parmCache.flush();
    }


  } //# end namespace DPPP
} //# end namespace LOFAR
