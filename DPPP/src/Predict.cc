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
#include <DPPP/Simulator.h>
#include <DPPP/Simulate.h>

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
        itsUseChannelFreq(parset.getBool (prefix + "usechannelfreq", true)),
        itsDebugLevel    (parset.getInt (prefix + "debuglevel", 0)),
        itsPatchList     ()
    {
      BBS::SourceDB sourceDB(BBS::ParmDBMeta("", itsSourceDBName), false);

      vector<string> sourcePatterns=parset.getStringVector(prefix + "sources",
                                                           vector<string>());
      vector<string> patchNames=makePatchList(sourceDB, sourcePatterns);

      itsPatchList = makePatches (sourceDB, patchNames, patchNames.size());

      // Rework patch list to contain a patch for every source
      if (!parset.getBool(prefix + "onebeamperpatch", true)) {
        itsPatchList = makeOnePatchPerComponent(itsPatchList);
      }

      itsSourceList = makeSourceList(itsPatchList);
    }

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
      const size_t nSt = info().antennaUsed().size();
      const size_t nCh = info().nchan();
      const size_t nCr = info().ncorr();

      itsUVW.resize(3,nSt);
      itsUVWSplitIndex = nsetupSplitUVW (info().nantenna(), info().getAnt1(),
                                         info().getAnt2());

      itsModelVis.resize(OpenMP::maxThreads());
      itsModelVisTmp.resize(OpenMP::maxThreads());
      itsBeamValues.resize(OpenMP::maxThreads());

      // Create the Measure ITRF conversion info given the array position.
      // The time and direction are filled in later.
      itsMeasConverters.resize(OpenMP::maxThreads());
      itsMeasFrames.resize(OpenMP::maxThreads());
      itsAntBeamInfo.resize(OpenMP::maxThreads());

      for (uint thread=0;thread<OpenMP::maxThreads();++thread) {
        itsModelVis[thread].resize(nCr,nCh,nBl);
        itsModelVisTmp[thread].resize(nCr,nCh,nBl);
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
      os << "  apply beam:         " << boolalpha << itsApplyBeam << endl;
      os << "   use channelfreq:   " << boolalpha << itsUseChannelFreq << endl;
      os << "Threads: "<<OpenMP::maxThreads()<<endl;
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
      DPBuffer buf(bufin);
      buf.getData().unique();
      Complex* data=buf.getData().data();
      RefRows refRows(buf.getRowNrs());

      buf.setUVW(itsInput->fetchUVW(buf, refRows, itsTimer));
      buf.setWeights(itsInput->fetchWeights(buf, refRows, itsTimer));
      buf.setFullResFlags(itsInput->fetchFullResFlags(buf, refRows, itsTimer));

      // Determine the various sizes.
      //const size_t nDr = itsPatchList.size();
      const size_t nSt = info().antennaUsed().size();
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = 4;
      const size_t nSamples = nBl * nCh * nCr;

      double time = buf.getTime();

      itsTimerPredict.start();

      nsplitUVW(itsUVWSplitIndex, itsBaselines, buf.getUVW(), itsUVW);

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
      StationResponse::vector3r_t srcdir;
      uint thread=OpenMP::threadNum();
      itsModelVis[thread]=dcomplex();
      itsModelVisTmp[thread]=dcomplex();
      Simulator simulator(itsPhaseRef, nSt, nBl, nCh, itsBaselines,
                          info().chanFreqs(), itsUVW, itsModelVisTmp[thread]);

      Patch::ConstPtr curPatch;
#pragma omp for
      for (uint i=0;i<itsSourceList.size();++i) {
        if (curPatch!=itsSourceList[i].second && curPatch!=0) {
          //Apply beam for curPatch, copy itsModelVisTmp to itsModelVis
          MDirection dir (MVDirection(curPatch->position()[0],
                                      curPatch->position()[1]),
                          MDirection::J2000);
          srcdir = dir2Itrf(dir,itsMeasConverters[thread]);
          if (itsApplyBeam) {
            applyBeam(info().chanFreqs(), time, itsModelVisTmp[thread],
                      srcdir, refdir, tiledir, itsAntBeamInfo[thread],
                      itsBeamValues[thread], itsUseChannelFreq);
          }

          //Add itsModelVisTmp to itsModelVis
          std::transform(itsModelVis[thread].data(),
                         itsModelVis[thread].data()+nSamples,
                         itsModelVisTmp[thread].data(),
                         itsModelVis[thread].data(), std::plus<dcomplex>());
          //threadoutput<<"thread "<<thread<<" has "<<itsModelVis[thread]<<endl;
          itsModelVisTmp[thread]=dcomplex();
        }
        simulator.simulate(itsSourceList[i].first);
        curPatch=itsSourceList[i].second;
      }

      if (curPatch!=0) {
        //Apply beam for curPatch, copy itsModelVisTmp to itsModelVis
        MDirection dir (MVDirection(curPatch->position()[0],
                                    curPatch->position()[1]),
                        MDirection::J2000);
        srcdir = dir2Itrf(dir,itsMeasConverters[thread]);
        if (itsApplyBeam) {
          applyBeam(info().chanFreqs(), time, itsModelVisTmp[thread],
                    srcdir, refdir, tiledir, itsAntBeamInfo[thread],
                    itsBeamValues[thread], itsUseChannelFreq);
        }

        //Add itsModelVisTmp to itsModelVis
        std::transform(itsModelVis[thread].data(),
                       itsModelVis[thread].data()+nSamples,
                       itsModelVisTmp[thread].data(),
                       itsModelVis[thread].data(), std::plus<dcomplex>());
        //threadoutput<<"thread "<<thread<<" has "<<itsModelVis[thread]<<endl;
        itsModelVisTmp[thread]=dcomplex();
      }
}

      //Add all thread model data to one buffer
      buf.getData()=Complex();
      for (uint thread=0;thread<OpenMP::maxThreads();++thread) {
        std::transform(data, data+nSamples, itsModelVis[thread].data(),
                       data, std::plus<dcomplex>());
      }

      itsTimerPredict.stop();

      itsTimer.stop();
      getNextStep()->process(buf);
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

    void Predict::applyBeam (const Vector<double>& chanFreqs, double time,
                             Cube<dcomplex>& data0,
                             const StationResponse::vector3r_t& srcdir,
                             const StationResponse::vector3r_t& refdir,
                             const StationResponse::vector3r_t& tiledir,
                             const vector<StationResponse::Station::Ptr>&
                               antBeamInfo,
                             vector<StationResponse::matrix22c_t>& beamValues,
                             bool useChannelFreq)  {
      // Get the beam values for each station.
      uint nCh = chanFreqs.size();
      uint nSt   = beamValues.size()/nCh;
      uint nBl   = info().nbaselines();

      if (!useChannelFreq) {
        for (size_t st=0; st<nSt; ++st) {
          antBeamInfo[st]->response (nCh, time, chanFreqs.cbegin(),
                                     srcdir, info().refFreq(), refdir,
                                     tiledir, &(beamValues[nCh*st]));
        }
      }

      // Apply the beam values of both stations to the predicted data.
      dcomplex tmp[4];
      for (size_t ch=0; ch<nCh; ++ch) {
        if (useChannelFreq) {
          for (size_t st=0; st<nSt; ++st) {
            antBeamInfo[st]->response (nCh, time, chanFreqs.cbegin(),
                                       srcdir, chanFreqs[ch], refdir,
                                       tiledir, &(beamValues[nCh*st]));
          }
        }
        for (size_t bl=0; bl<nBl; ++bl) {
          dcomplex* data=&data0(0,ch,bl);
          StationResponse::matrix22c_t *left =
              &(beamValues[nCh * info().getAnt1()[bl]]);
          StationResponse::matrix22c_t *right=
              &(beamValues[nCh * info().getAnt2()[bl]]);
          dcomplex l[] = {left[ch][0][0], left[ch][0][1],
                          left[ch][1][0], left[ch][1][1]};
          // Form transposed conjugate of right.
          dcomplex r[] = {conj(right[ch][0][0]), conj(right[ch][1][0]),
                          conj(right[ch][0][1]), conj(right[ch][1][1])};
          // left*data
          tmp[0] = l[0] * data[0] + l[1] * data[2];
          tmp[1] = l[0] * data[1] + l[1] * data[3];
          tmp[2] = l[2] * data[0] + l[3] * data[2];
          tmp[3] = l[2] * data[1] + l[3] * data[3];
          // data*conj(right)
          data[0] = tmp[0] * r[0] + tmp[1] * r[2];
          data[1] = tmp[0] * r[1] + tmp[1] * r[3];
          data[2] = tmp[2] * r[0] + tmp[3] * r[2];
          data[3] = tmp[2] * r[1] + tmp[3] * r[3];
        }
      }
    }


    void Predict::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } //# end namespace
}
