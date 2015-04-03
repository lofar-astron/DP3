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

#include <lofar_config.h>
#include <DPPP/Predict.h>

#include <iostream>
//#include <iomanip>
#include <Common/ParameterSet.h>
#include <Common/Timer.h>
#include <Common/OpenMP.h>
#include <ParmDB/ParmDBMeta.h>
#include <ParmDB/PatchInfo.h>
#include <DPPP/DPInfo.h>
#include <DPPP/FlagCounter.h>
#include <DPPP/Position.h>
#include <DPPP/ApplyBeam.h>
#include <DPPP/Simulator.h>

#include <DPPP/Stokes.h>
#include <DPPP/PointSource.h>
#include <DPPP/GaussianSource.h>
#include <ParmDB/SourceDB.h>

#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Quanta/Quantum.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MeasConvert.h>
#include <tables/Tables/RefRows.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using namespace casa;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    Predict::Predict (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix)
      : itsInput         (input),
        itsName          (prefix),
        itsSourceDBName (parset.getString (prefix + "sourcedb")),
        itsApplyBeam     (parset.getBool (prefix + "usebeammodel", false)),
        itsDebugLevel    (parset.getInt (prefix + "debuglevel", 0)),
        itsPatchList     ()
    {
      BBS::SourceDB sourceDB(BBS::ParmDBMeta("", itsSourceDBName), false);

      vector<string> sourcePatterns=parset.getStringVector(prefix + "sources",
                                                           vector<string>());
      vector<string> patchNames=makePatchList(sourceDB, sourcePatterns);
      itsPatchList = makePatches (sourceDB, patchNames, patchNames.size());

      if (itsApplyBeam) {
        itsUseChannelFreq=parset.getBool (prefix + "usechannelfreq", true);
        itsOneBeamPerPatch=parset.getBool (prefix + "onebeamperpatch", true);

        // Rework patch list to contain a patch for every source
        if (!itsOneBeamPerPatch) {
          itsPatchList = makeOnePatchPerComponent(itsPatchList);
        }
      }


      itsSourceList = makeSourceList(itsPatchList);
    }

    Predict::Predict()
    {}

    Predict::~Predict()
    {}

    void Predict::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setNeedWrite();

      uint nBl=info().nbaselines();
      for (uint i=0; i<nBl; ++i) {
        itsBaselines.push_back (Baseline(info().getAnt1()[i],
                                         info().getAnt2()[i]));
      }

      MDirection dirJ2000(MDirection::Convert(infoIn.phaseCenter(),
                                              MDirection::J2000)());
      Quantum<Vector<Double> > angles = dirJ2000.getAngle();
      itsPhaseRef = Position(angles.getBaseValue()[0],
                             angles.getBaseValue()[1]);

      //const size_t nDr = itsPatchList.size();
      const size_t nSt = info().nantenna();
      const size_t nCh = info().nchan();
      const size_t nCr = info().ncorr();

      itsUVW.resize(3,nSt);
      itsUVWSplitIndex = nsetupSplitUVW (info().nantenna(), info().getAnt1(),
                                         info().getAnt2());

      itsModelVis.resize(OpenMP::maxThreads());
      itsModelVisPatch.resize(OpenMP::maxThreads());
      itsBeamValues.resize(OpenMP::maxThreads());

      // Create the Measure ITRF conversion info given the array position.
      // The time and direction are filled in later.
      itsMeasConverters.resize(OpenMP::maxThreads());
      itsMeasFrames.resize(OpenMP::maxThreads());
      itsAntBeamInfo.resize(OpenMP::maxThreads());

      for (uint thread=0;thread<OpenMP::maxThreads();++thread) {
        itsModelVis[thread].resize(nCr,nCh,nBl);
        itsModelVisPatch[thread].resize(nCr,nCh,nBl);
        itsBeamValues[thread].resize(nSt*nCh);
        itsMeasFrames[thread].set (info().arrayPosCopy());
        itsMeasFrames[thread].set (MEpoch(MVEpoch(info().startTime()/86400),
                                          MEpoch::UTC));
        itsMeasConverters[thread].set (MDirection::J2000,
                     MDirection::Ref(MDirection::ITRF, itsMeasFrames[thread]));
        itsInput->fillBeamInfo (itsAntBeamInfo[thread], info().antennaNames());
      }
    }

    void Predict::show (std::ostream& os) const
    {
      os << "Predict " << itsName << endl;
      os << "  sourcedb:           " << itsSourceDBName << endl;
      os << "   number of patches: " << itsPatchList.size() << endl;
      os << "   number of sources: " << itsSourceList.size() << endl;
      os << "  apply beam:         " << boolalpha << itsApplyBeam << endl;
      if (itsApplyBeam) {
        os << "   use channelfreq:   " << boolalpha << itsUseChannelFreq << endl;
        os << "   one beam per patch:" << boolalpha << itsOneBeamPerPatch << endl;
      }
      os << "  threads:            "<<OpenMP::maxThreads()<<endl;
    }

    void Predict::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " Predict " << itsName << endl;
    }

    bool Predict::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      itsBuffer.copy (bufin);
      Complex* data=itsBuffer.getData().data();
      itsInput->fetchUVW(bufin, itsBuffer, itsTimer);
      itsInput->fetchWeights(bufin, itsBuffer, itsTimer);
      itsInput->fetchFullResFlags(bufin, itsBuffer, itsTimer);

      // Determine the various sizes.
      //const size_t nDr = itsPatchList.size();
      const size_t nSt = info().nantenna();
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = 4;
      const size_t nSamples = nBl * nCh * nCr;

      double time = itsBuffer.getTime();

      itsTimerPredict.start();

      nsplitUVW(itsUVWSplitIndex, itsBaselines, itsBuffer.getUVW(), itsUVW);

      //Set up directions for beam evaluation
      StationResponse::vector3r_t refdir, tiledir;

      if (itsApplyBeam) {
        for (uint thread=0;thread<OpenMP::maxThreads();++thread) {
          itsMeasFrames[thread].resetEpoch (MEpoch(MVEpoch(time/86400),
                                                   MEpoch::UTC));
          //Do a conversion on all threads, because converters are not
          //thread safe and apparently need to be used at least once
          refdir  = dir2Itrf(info().delayCenter(),itsMeasConverters[thread]);
          tiledir = dir2Itrf(info().tileBeamDir(),itsMeasConverters[thread]);
        }
      }

#pragma omp parallel
{
      uint thread=OpenMP::threadNum();
      itsModelVis[thread]=dcomplex();
      itsModelVisPatch[thread]=dcomplex();

      //When applying beam, simulate into patch vector,
      Cube<dcomplex>& simulatedest=(itsApplyBeam?itsModelVisPatch[thread]
                                                :itsModelVis[thread]);

      Simulator simulator(itsPhaseRef, nSt, nBl, nCh, itsBaselines,
                          info().chanFreqs(), itsUVW, simulatedest);

      Patch::ConstPtr curPatch;
#pragma omp for
      for (uint i=0;i<itsSourceList.size();++i) {
        if (itsApplyBeam && curPatch!=itsSourceList[i].second && curPatch!=0) {
          addBeamToData (curPatch, time, refdir, tiledir, thread, nSamples,
                         itsModelVisPatch[thread].data());
        }
        simulator.simulate(itsSourceList[i].first);
        curPatch=itsSourceList[i].second;
      }

      if (itsApplyBeam && curPatch!=0) {
        addBeamToData (curPatch, time, refdir, tiledir, thread, nSamples,
                       itsModelVisPatch[thread].data());
      }
}

      //Add all thread model data to one buffer
      itsBuffer.getData()=Complex();
      for (uint thread=0;thread<OpenMP::maxThreads();++thread) {
        std::transform(data, data+nSamples, itsModelVis[thread].data(),
                       data, std::plus<dcomplex>());
      }

      itsTimerPredict.stop();

      itsTimer.stop();
      getNextStep()->process(itsBuffer);
      return false;
    }

    StationResponse::vector3r_t Predict::dir2Itrf (const MDirection& dir,
                                      MDirection::Convert& measConverter) {
      const MDirection& itrfDir = measConverter(dir);
      const Vector<Double>& itrf = itrfDir.getValue().getValue();
      StationResponse::vector3r_t vec;
      vec[0] = itrf[0];
      vec[1] = itrf[1];
      vec[2] = itrf[2];
      return vec;
    }

    void Predict::addBeamToData (Patch::ConstPtr patch, double time,
                                 const StationResponse::vector3r_t& refdir,
                                 const StationResponse::vector3r_t& tiledir,
                                 uint thread, uint nSamples, dcomplex* data0) {
      //Apply beam for a patch, add result to itsModelVis
      MDirection dir (MVDirection(patch->position()[0],
                                  patch->position()[1]),
                      MDirection::J2000);
      StationResponse::vector3r_t srcdir = dir2Itrf(dir,itsMeasConverters[thread]);

      ApplyBeam::applyBeam(info(), time, data0, srcdir, refdir,
                           tiledir, itsAntBeamInfo[thread],
                           itsBeamValues[thread], itsUseChannelFreq, false);

      //Add temporary buffer to itsModelVis
      std::transform(itsModelVis[thread].data(),
                     itsModelVis[thread].data()+nSamples,
                     data0,
                     itsModelVis[thread].data(), std::plus<dcomplex>());
      //threadoutput<<"thread "<<thread<<" has "<<itsModelVis[thread]<<endl;
      itsModelVisPatch[thread]=dcomplex();
    }

    void Predict::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } //# end namespace
}
