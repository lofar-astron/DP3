//# GainCal.cc: DPPP step class to predict visibilities
//# Copyright (C) 2013
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
//# $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include "Predict.h"

#include <iostream>

#include "../Common/ParameterSet.h"
#include "../Common/ThreadPool.h"
#include "../Common/Timer.h"
#include "../Common/StreamUtil.h"

#include "../ParmDB/ParmDBMeta.h"
#include "../ParmDB/PatchInfo.h"

#include "DPInfo.h"
#include "Exceptions.h"
#include "FlagCounter.h"
#include "Position.h"
#include "ApplyBeam.h"
#include "Simulate.h"
#include "Simulator.h"
#include "Stokes.h"
#include "PointSource.h"
#include "GaussianSource.h"

#include "../ParmDB/SourceDB.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/tables/Tables/RefRows.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/algorithm/string/case_conv.hpp>

using namespace casacore;
using namespace DP3::BBS;

namespace DP3 {
  namespace DPPP {

    Predict::Predict (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix) :
      itsThreadPool(nullptr),
      itsMeasuresMutex(nullptr)
    {
      init(input, parset, prefix, parset.getStringVector(prefix + "sources",
                                                         vector<string>()));
    }

    Predict::Predict (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix,
                      const vector<string>& sourcePatterns) :
      itsThreadPool(nullptr),
      itsMeasuresMutex(nullptr)
    {
      init(input, parset, prefix, sourcePatterns);
    }

    void Predict::init(DPInput* input,
              const ParameterSet& parset,
              const string& prefix, const vector<string>& sourcePatterns) {
      itsInput = input;
      itsName = prefix;
      itsSourceDBName = parset.getString (prefix + "sourcedb");
      setOperation(parset.getString (prefix + "operation", "replace"));
#ifdef HAVE_LOFAR_BEAM
      itsApplyBeam = parset.getBool (prefix + "usebeammodel", false);
#endif
      itsDebugLevel = parset.getInt (prefix + "debuglevel", 0);
      itsPatchList = vector<Patch::ConstPtr> ();

      assert(File(itsSourceDBName).exists());
      BBS::SourceDB sourceDB(BBS::ParmDBMeta("", itsSourceDBName), false);

      // Save directions specifications to pass to applycal
      stringstream ss;
      ss << sourcePatterns;
      itsDirectionsStr = ss.str();

      vector<string> patchNames=makePatchList(sourceDB, sourcePatterns);
      itsPatchList = makePatches (sourceDB, patchNames, patchNames.size());

#ifdef HAVE_LOFAR_BEAM
      if (itsApplyBeam) {
        itsUseChannelFreq=parset.getBool (prefix + "usechannelfreq", true);
        itsOneBeamPerPatch=parset.getBool (prefix + "onebeamperpatch", false);

        string mode=boost::to_lower_copy(parset.getString(prefix + "beammode","default"));
        assert (mode=="default" || mode=="array_factor" || mode=="element");
        if (mode=="default") {
          itsBeamMode=FullBeamCorrection;
        } else if (mode=="array_factor") {
          itsBeamMode=ArrayFactorBeamCorrection;
        } else if (mode=="element") {
          itsBeamMode=ElementBeamCorrection;
        } else {
          throw Exception("Beammode should be DEFAULT, ARRAY_FACTOR or ELEMENT");
        }

        // Rework patch list to contain a patch for every source
        if (!itsOneBeamPerPatch) {
          itsPatchList = makeOnePatchPerComponent(itsPatchList);
        }
      }
#endif

      // If called from h5parmpredict, applycal gets set by that step,
      // so must not be read from parset
      if (parset.isDefined(prefix + "applycal.parmdb") ||
          parset.isDefined(prefix + "applycal.steps")) {
        setApplyCal(input, parset, prefix + "applycal.");
      } else {
        itsDoApplyCal=false;
      }

      itsSourceList = makeSourceList(itsPatchList);

      // Determine whether any sources are polarized. If not, enable Stokes-I-
      // only mode (note that this mode cannot be used with itsApplyBeam)
#ifdef HAVE_LOFAR_BEAM
      if (itsApplyBeam) {
        itsStokesIOnly = false;
      } else {
        itsStokesIOnly = !checkPolarized(sourceDB, patchNames, patchNames.size());
      }
#else
      itsStokesIOnly = !checkPolarized(sourceDB, patchNames, patchNames.size());
#endif
    }

    void Predict::setApplyCal(DPInput* input,
                              const ParameterSet& parset,
                              const string& prefix) {
      itsDoApplyCal=true;
      itsApplyCalStep=ApplyCal(input, parset, prefix, true,
                               itsDirectionsStr);
      assert(!(itsOperation!="replace" &&
               parset.getBool(prefix + "applycal.updateweights", false)));
      itsResultStep=new ResultStep();
      itsApplyCalStep.setNextStep(DPStep::ShPtr(itsResultStep));
    }

    Predict::~Predict()
    {}
    
    void Predict::initializeThreadData()
    {
      const size_t nBl=info().nbaselines();
      const size_t nSt = info().nantenna();
      const size_t nCh = info().nchan();
      const size_t nCr = info().ncorr();
      const size_t nThreads = getInfo().nThreads();
      
      itsUVW.resize(3, nSt);
      itsUVWSplitIndex = nsetupSplitUVW (info().nantenna(), info().getAnt1(),
                                         info().getAnt2());
      
      itsModelVis.resize(nThreads);
      itsModelVisPatch.resize(nThreads);
#ifdef HAVE_LOFAR_BEAM
      itsBeamValues.resize(nThreads);
      itsAntBeamInfo.resize(nThreads);
      // Create the Measure ITRF conversion info given the array position.
      // The time and direction are filled in later.
      itsMeasConverters.resize(nThreads);
      itsMeasFrames.resize(nThreads);
#endif

      for (uint thread=0; thread<nThreads; ++thread) {
        if (itsStokesIOnly) {
          itsModelVis[thread].resize(1,nCh,nBl);
        } else {
          itsModelVis[thread].resize(nCr,nCh,nBl);
        }
#ifdef HAVE_LOFAR_BEAM
        if (itsApplyBeam) {
          itsModelVisPatch[thread].resize(nCr,nCh,nBl);
          itsBeamValues[thread].resize(nSt*nCh);
          itsMeasFrames[thread].set (info().arrayPosCopy());
          itsMeasFrames[thread].set (MEpoch(MVEpoch(info().startTime()/86400),
                                            MEpoch::UTC));
          itsMeasConverters[thread].set (MDirection::J2000,
                                         MDirection::Ref(MDirection::ITRF, itsMeasFrames[thread]));
          itsInput->fillBeamInfo (itsAntBeamInfo[thread], info().antennaNames());
        }
#endif
      }
    }

    void Predict::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();

      const size_t nBl=info().nbaselines();
      for (size_t i=0; i!=nBl; ++i) {
        itsBaselines.push_back (Baseline(info().getAnt1()[i],
                                         info().getAnt2()[i]));
      }

      MDirection dirJ2000(MDirection::Convert(infoIn.phaseCenter(),
                                              MDirection::J2000)());
      Quantum<Vector<Double> > angles = dirJ2000.getAngle();
      itsPhaseRef = Position(angles.getBaseValue()[0],
                             angles.getBaseValue()[1]);

      initializeThreadData();

      if (itsDoApplyCal) {
        info()=itsApplyCalStep.setInfo(info());
      }
    }

    std::pair<double, double> Predict::getFirstDirection() const {
      std::pair<double, double> res;
      res.first = itsPatchList[0]->position()[0];
      res.second = itsPatchList[0]->position()[1];
      return res;
    }

    void Predict::setOperation(const std::string& operation) {
      itsOperation=operation;
      assert(itsOperation=="replace" || itsOperation=="add" ||
             itsOperation=="subtract");
    }

    void Predict::show (std::ostream& os) const
    {
      os << "Predict " << itsName << '\n';
      os << "  sourcedb:           " << itsSourceDBName << '\n';
      os << "   number of patches: " << itsPatchList.size() << '\n';
      os << "   number of sources: " << itsSourceList.size() << '\n';
      os << "   all unpolarized:   " << boolalpha << itsStokesIOnly << '\n';
#ifdef HAVE_LOFAR_BEAM
      os << "  apply beam:         " << boolalpha << itsApplyBeam << '\n';
      if (itsApplyBeam) {
        os << "   mode:              ";
        if (itsBeamMode==FullBeamCorrection)
          os<<"default";
        else if (itsBeamMode==ArrayFactorBeamCorrection)
          os<<"array_factor";
        else os<<"element";
        os << '\n';
        os << "   use channelfreq:   " << boolalpha << itsUseChannelFreq << '\n';
        os << "   one beam per patch:" << boolalpha << itsOneBeamPerPatch << '\n';
      }
#endif
      os << "  operation:          " << itsOperation << '\n';
      os << "  threads:            " << getInfo().nThreads() << '\n';
      if (itsDoApplyCal) {
        itsApplyCalStep.show(os);
      }
    }

    void Predict::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " Predict " << itsName << '\n';
    }

    bool Predict::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      itsTempBuffer.copy (bufin);
      itsInput->fetchUVW(bufin, itsTempBuffer, itsTimer);
      itsInput->fetchWeights(bufin, itsTempBuffer, itsTimer);

      // Determine the various sizes.
      //const size_t nDr = itsPatchList.size();
      const size_t nSt = info().nantenna();
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = info().ncorr();
      const size_t nSamples = nBl * nCh * nCr;

      itsTimerPredict.start();

      nsplitUVW(itsUVWSplitIndex, itsBaselines, itsTempBuffer.getUVW(), itsUVW);

#ifdef HAVE_LOFAR_BEAM
      double time = itsTempBuffer.getTime();
      //Set up directions for beam evaluation
      LOFAR::StationResponse::vector3r_t refdir, tiledir;

      if (itsApplyBeam)
      {
        // Because multiple predict steps might be predicting simultaneously, and
        // Casacore is not thread safe, this needs synchronization.
        std::unique_lock<std::mutex> lock;
        if(itsMeasuresMutex != nullptr)
          lock = std::unique_lock<std::mutex>(*itsMeasuresMutex);
        for (uint thread=0; thread!=getInfo().nThreads(); ++thread) {
          itsMeasFrames[thread].resetEpoch (MEpoch(MVEpoch(time/86400),
                                                   MEpoch::UTC));
          //Do a conversion on all threads, because converters are not
          //thread safe and apparently need to be used at least once
          refdir  = dir2Itrf(info().delayCenter(), itsMeasConverters[thread]);
          tiledir = dir2Itrf(info().tileBeamDir(), itsMeasConverters[thread]);
        }
      }
#endif

      std::unique_ptr<ThreadPool> localThreadPool;
      ThreadPool* pool = itsThreadPool;
      if(pool == nullptr)
      {
        // If no ThreadPool was specified, we create a temporary one just
        // for executation of this part.
        localThreadPool.reset(new ThreadPool(info().nThreads()));
        pool = localThreadPool.get();
      }
      else {
        if(pool->NThreads() != info().nThreads())
          throw std::runtime_error("Thread pool has inconsistent number of threads!");
      }
      std::vector<Simulator> simulators;
      simulators.reserve(pool->NThreads());
      for(size_t thread=0; thread!=pool->NThreads(); ++thread)
      {
        itsModelVis[thread]=dcomplex();
        itsModelVisPatch[thread]=dcomplex();

#ifdef HAVE_LOFAR_BEAM
        //When applying beam, simulate into patch vector
        Cube<dcomplex>& simulatedest=(itsApplyBeam ? itsModelVisPatch[thread]
          : itsModelVis[thread]);
#else
        Cube<dcomplex>& simulatedest=itsModelVis[thread];
#endif
        simulators.emplace_back(itsPhaseRef, nSt, nBl, nCh, itsBaselines,
          info().chanFreqs(), itsUVW, simulatedest,
          itsStokesIOnly);
      }
      std::vector<Patch::ConstPtr> curPatches(pool->NThreads());
      
      pool->For(0, itsSourceList.size(), [&](size_t iter, size_t thread) {
        // Keep on predicting, only apply beam when an entire patch is done
        Patch::ConstPtr& curPatch = curPatches[thread];
#ifdef HAVE_LOFAR_BEAM
        if (itsApplyBeam && curPatch!=itsSourceList[iter].second && curPatch!=nullptr) {
          addBeamToData (curPatch, time, refdir, tiledir, thread, nSamples,
            itsModelVisPatch[thread].data());
        }
#endif
        simulators[thread].simulate(itsSourceList[iter].first);
        curPatch=itsSourceList[iter].second;
      });
#ifdef HAVE_LOFAR_BEAM
      // Apply beam to the last patch
      for(size_t thread=0; thread!=pool->NThreads(); ++thread)
      {
        if (itsApplyBeam && curPatches[thread]!=nullptr) {
          addBeamToData (curPatches[thread], time, refdir, tiledir, thread, nSamples,
            itsModelVisPatch[thread].data());
        }
      }
#endif

      // Add all thread model data to one buffer
      itsTempBuffer.getData()=Complex();
      Complex* tdata=itsTempBuffer.getData().data();
      for (uint thread=0; thread<pool->NThreads(); ++thread) {
        if (itsStokesIOnly) {
          for (uint i=0,j=0;i<nSamples;i+=nCr,j++) {
            tdata[i] += itsModelVis[thread].data()[j];
            tdata[i+nCr-1] += itsModelVis[thread].data()[j];
          }
        } else {
          std::transform(tdata, tdata+nSamples, itsModelVis[thread].data(),
                         tdata, std::plus<dcomplex>());
        }
      }

      // Call ApplyCal step
      if (itsDoApplyCal) {
        if(itsMeasuresMutex == nullptr)
          itsApplyCalStep.process(itsTempBuffer);
        else
          itsApplyCalStep.process(itsTempBuffer, itsMeasuresMutex);
        itsTempBuffer=itsResultStep->get();
        tdata=itsTempBuffer.getData().data();
      }

      // Put predict result from temp buffer into the 'real' buffer
      if (itsOperation=="replace") {
        itsBuffer=itsTempBuffer;
      } else {
        itsBuffer.copy(bufin);
        Complex* data=itsBuffer.getData().data();
        if (itsOperation=="add") {
          std::transform(data, data+nSamples, tdata,
                         data, std::plus<dcomplex>());
        } else if (itsOperation=="subtract") {
          std::transform(data, data+nSamples, tdata,
                         data, std::minus<dcomplex>());
        }
      }

      itsTimerPredict.stop();

      itsTimer.stop();
      getNextStep()->process(itsBuffer);
      return false;
    }

#ifdef HAVE_LOFAR_BEAM
    LOFAR::StationResponse::vector3r_t Predict::dir2Itrf (const MDirection& dir,
                                      MDirection::Convert& measConverter) {
      const MDirection& itrfDir = measConverter(dir);
      const Vector<Double>& itrf = itrfDir.getValue().getValue();
      LOFAR::StationResponse::vector3r_t vec;
      vec[0] = itrf[0];
      vec[1] = itrf[1];
      vec[2] = itrf[2];
      return vec;
    }

    void Predict::addBeamToData (Patch::ConstPtr patch, double time,
                                 const LOFAR::StationResponse::vector3r_t& refdir,
                                 const LOFAR::StationResponse::vector3r_t& tiledir,
                                 uint thread, uint nSamples, dcomplex* data0) {
      //Apply beam for a patch, add result to itsModelVis
      MDirection dir (MVDirection(patch->position()[0],
                                  patch->position()[1]),
                      MDirection::J2000);
      LOFAR::StationResponse::vector3r_t srcdir = dir2Itrf(dir, itsMeasConverters[thread]);

      float* dummyweight = 0;

      ApplyBeam::applyBeam(info(), time, data0, dummyweight, srcdir, refdir,
                           tiledir, itsAntBeamInfo[thread],
                           itsBeamValues[thread], itsUseChannelFreq, false,
                           itsBeamMode, false);

      //Add temporary buffer to itsModelVis
      std::transform(itsModelVis[thread].data(),
                     itsModelVis[thread].data()+nSamples,
                     data0,
                     itsModelVis[thread].data(), std::plus<dcomplex>());
      itsModelVisPatch[thread]=dcomplex();
    }
#endif

    void Predict::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } //# end namespace
}
