//# GainCal.cc: DPPP step class to Predict1 visibilities
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
#include <DPPP/Predict1-poging1.h>
#include <DPPP/Simulate.h>
#include <DPPP/CursorUtilCasa.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/SourceDBUtil.h>
#include <DPPP/MSReader.h>
#include <ParmDB/SourceDB.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <Common/LofarLogger.h>
#include <Common/OpenMP.h>

#include <fstream>
#include <ctime>

#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCDirection.h>
#include <casa/OS/File.h>

#include <vector>
#include <algorithm>

#include <iostream>
#include <iomanip>

using namespace casa;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    Predict1::Predict1 (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix)
      : itsInput         (input),
        itsName          (prefix),
        itsSourceDBName  (""),
        itsApplyBeam     (parset.getBool (prefix + "usebeammodel", false)),
        itsOneBeamPerPatch  (parset.getBool (prefix + "onebeamperpatch", true)),
        itsUseChannelFreq(parset.getBool (prefix + "usechannelfreq", true)),
        itsDebugLevel    (parset.getInt (prefix + "debuglevel", 0)),
        itsBaselines     (),
        itsThreadStorage (),
        itsPatchList     ()
    {
      itsSourceDBName = parset.getString (prefix + "sourcedb","");
      BBS::SourceDB sourceDB(BBS::ParmDBMeta("", itsSourceDBName), false);

      vector<PatchInfo> patchInfo=sourceDB.getPatchInfo();
      vector<string> patchNames;

      vector<string> sourcePatterns=parset.getStringVector(prefix + "sources",
                                                           vector<string>());
      patchNames=makePatchList(sourceDB, sourcePatterns);

      itsPatchList = makePatches (sourceDB, patchNames, patchNames.size());
      if (!itsOneBeamPerPatch) {
        itsPatchList = makeOnePatchPerComponent(itsPatchList);
      }
    }

    Predict1::~Predict1()
    {}

    void Predict1::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setNeedWrite();

      uint nBl=info().nbaselines();
      for (uint i=0; i<nBl; ++i) {
        itsBaselines.push_back (Baseline(info().getAnt1()[i],
                                         info().getAnt2()[i]));
      }

      MDirection dirJ2000(MDirection::Convert(infoIn.phaseCenterCopy(),
                                              MDirection::J2000)());
      Quantum<Vector<Double> > angles = dirJ2000.getAngle();
      itsPhaseRef = Position(angles.getBaseValue()[0],
                             angles.getBaseValue()[1]);

      const size_t nDr = itsPatchList.size();
      const size_t nSt = info().antennaUsed().size();
      const size_t nCh = info().nchan();

      const size_t nThread=1;//OpenMP::maxThreads();

      itsThreadStorage.resize(nThread);
      for(vector<ThreadPrivateStorage>::iterator it = itsThreadStorage.begin(),
          end = itsThreadStorage.end(); it != end; ++it)
      {
        initThreadPrivateStorage(*it, nDr, nSt, nBl, nCh, nCh);
      }

      itsInput->fillBeamInfo (itsAntBeamInfo, info().antennaNames());
    }

    StationResponse::vector3r_t Predict1::dir2Itrf (const MDirection& dir,
                                                   MDirection::Convert& converter) const
    {
      const MDirection& itrfDir = converter(dir);
      const Vector<Double>& itrf = itrfDir.getValue().getValue();
      StationResponse::vector3r_t vec;
      vec[0] = itrf[0];
      vec[1] = itrf[1];
      vec[2] = itrf[2];
      return vec;
    }

    void Predict1::show (std::ostream& os) const
    {
      os << "Precit " << itsName << endl;
      os << "  sourcedb:           " << itsSourceDBName << endl;
      os << "   number of patches: " << itsPatchList.size() << endl;
      os << "  apply beam:         " << boolalpha << itsApplyBeam << endl;
      os << "   beam per patch:    " << boolalpha << itsOneBeamPerPatch << endl;
      os << "   use channelfreq:   " << boolalpha << itsUseChannelFreq << endl;
    }

    void Predict1::showTimings (std::ostream& os, double duration) const
    {
      double totaltime=itsTimer.getElapsed();
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " Predict1 " << itsName << endl;
    }

    bool Predict1::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      DPBuffer buf(bufin);
      buf.getData().unique();
      RefRows refRows(buf.getRowNrs());

      buf.setUVW(itsInput->fetchUVW(buf, refRows, itsTimer));
      buf.setWeights(itsInput->fetchWeights(buf, refRows, itsTimer));
      buf.setFullResFlags(itsInput->fetchFullResFlags(buf, refRows, itsTimer));

      // Determine the various sizes.
      const size_t nDr = itsPatchList.size();
      const size_t nSt = info().antennaUsed().size();
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = 4;
      const size_t nSamples = nBl * nCh * nCr;
      // Define various cursors to iterate through arrays.
      const_cursor<double> cr_freq = casa_const_cursor(info().chanFreqs());
      const_cursor<Baseline> cr_baseline(&(itsBaselines[0]));

      const size_t thread = 0;//OpenMP::threadNum();

      Complex* data=buf.getData().data();
      float* weight = buf.getWeights().data();
      const Bool* flag=buf.getFlags().data();

      // Simulate.
      //
      // Model visibilities for each direction of interest will be computed
      // and stored.

      itsTimerPredict1.start();

      ThreadPrivateStorage &storage = itsThreadStorage[thread];
      double time = buf.getTime();

      size_t stride_uvw[2] = {1, 3};
      cursor<double> cr_uvw_split(&(storage.uvw[0]), 2, stride_uvw);

      size_t stride_model[3] = {1, nCr, nCr * nCh};
      fill(storage.model.begin(), storage.model.end(), dcomplex());

      const_cursor<double> cr_uvw = casa_const_cursor(buf.getUVW());
      splitUVW(nSt, nBl, cr_baseline, cr_uvw, cr_uvw_split);
      cursor<dcomplex> cr_model(&(storage.model_patch[0]), 3, stride_model);

      // Convert the directions to ITRF for the given time.
      storage.measFrame.resetEpoch (MEpoch(MVEpoch(time/86400), MEpoch::UTC));
      StationResponse::vector3r_t refdir = dir2Itrf(info().delayCenter(),storage.measConverter);
      StationResponse::vector3r_t tiledir = dir2Itrf(info().tileBeamDir(),storage.measConverter);

      //#pragma omp parallel for
      for(size_t dr = 0; dr < nDr; ++dr) {
        fill(storage.model_patch.begin(), storage.model_patch.end(), dcomplex());

        simulate(itsPhaseRef, itsPatchList[dr], nSt, nBl, nCh, cr_baseline,
                 cr_freq, cr_uvw_split, cr_model);

        /*applyBeam(time, itsPatchList[dr]->position(), itsApplyBeam,
                  info().chanFreqs(), &(itsThreadStorage[thread].model_patch[0]),
                  refdir, tiledir, &(itsThreadStorage[thread].beamvalues[0]),
                  storage.measConverter);*/

        for (size_t i=0; i<itsThreadStorage[thread].model_patch.size();++i) {
          itsThreadStorage[thread].model[i]+=
              itsThreadStorage[thread].model_patch[i];
        }
      }

      itsTimerPredict1.stop();
      //copy result of model to data
      if (itsOperation=="Predict1") {
        copy(storage.model.begin(),storage.model.begin()+nSamples,data);
      }

      itsTimer.stop();
      getNextStep()->process(buf);
      return false;
    }


    void GainCal::applyBeam (double time, const Position& pos, bool apply,
                             const Vector<double>& chanFreqs, dcomplex* data0,
                             StationResponse::vector3r_t& refdir,
                             StationResponse::vector3r_t& tiledir,
                             StationResponse::matrix22c_t* beamvalues,
                             casa::MDirection::Convert& converter)
    {
      if (! apply) {
        return;
      }

      MDirection dir (MVDirection(pos[0], pos[1]), MDirection::J2000);
      StationResponse::vector3r_t srcdir = dir2Itrf(dir,converter);
      // Get the beam values for each station.
      uint nchan = chanFreqs.size();
      uint nSt   = info().antennaUsed().size();
      uint nBl   = info().nbaselines();

      if (!itsUseChannelFreq) {
        for (size_t st=0; st<nSt; ++st) {
          itsAntBeamInfo[st]->response (nchan, time, chanFreqs.cbegin(),
                                        srcdir, info().refFreq(), refdir,
                                        tiledir, &(beamvalues[nchan*st]));
        }
      }

      // Apply the beam values of both stations to the Predict1ed data.
      dcomplex tmp[4];
      for (size_t ch=0; ch<nchan; ++ch) {
        if (itsUseChannelFreq) {
          for (size_t st=0; st<nSt; ++st) {
            itsAntBeamInfo[st]->response (nchan, time, chanFreqs.cbegin(),
                                          srcdir, chanFreqs[ch], refdir,
                                          tiledir, &(beamvalues[nchan*st]));
          }
        }
        for (size_t bl=0; bl<nBl; ++bl) {
          dcomplex* data=data0+bl*4*nchan + ch*4; //TODO
          StationResponse::matrix22c_t *left =
              &(beamvalues[nchan * info().getAnt1()[bl]]);
          StationResponse::matrix22c_t *right=
              &(beamvalues[nchan * info().getAnt2()[bl]]);
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

    void Predict1::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }


  } //# end namespace
}
