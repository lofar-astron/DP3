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
#include <ParmDB/SourceDB.h>
#include <DPPP/CursorUtilCasa.h>
#include <DPPP/DPInfo.h>
#include <DPPP/FlagCounter.h>
#include <DPPP/Position.h>
#include <DPPP/Simulator.h>
#include <DPPP/Simulate.h>

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

      MDirection dirJ2000(MDirection::Convert(infoIn.phaseCenterCopy(),
                                              MDirection::J2000)());
      Quantum<Vector<Double> > angles = dirJ2000.getAngle();
      itsPhaseRef = Position(angles.getBaseValue()[0],
                             angles.getBaseValue()[1]);

      //const size_t nDr = itsPatchList.size();
      //const size_t nSt = info().antennaUsed().size();
      const size_t nCh = info().nchan();
      const size_t nCr = info().ncorr();

      itsUVWSplitIndex = nsetupSplitUVW (info().nantenna(), info().getAnt1(),
                                         info().getAnt2());

      itsModelVis.resize(nCr,nCh,nBl);
      //
      //itsInput->fillBeamInfo (itsAntBeamInfo, info().antennaNames());
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
      const size_t nDr = itsPatchList.size();
      const size_t nSt = info().antennaUsed().size();
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = 4;
      const size_t nSamples = nBl * nCh * nCr;

      itsTimerPredict.start();

      itsUVW.resize(3,nSt);
      nsplitUVW(itsUVWSplitIndex, itsBaselines, buf.getUVW(), itsUVW);

      itsModelVis=dcomplex();

      Simulator simulator(itsPhaseRef, nSt, nBl, nCh, itsBaselines,
                          info().chanFreqs(), itsUVW, itsModelVis);

      for (size_t dr=0; dr<nDr; dr++) {
        for(size_t i = 0; i < itsPatchList[dr]->nComponents(); ++i) {
          simulator.simulate(itsPatchList[dr]->component(i));
        }
      }

      copy(itsModelVis.data(),itsModelVis.data()+nSamples,data);

      itsTimerPredict.stop();

      itsTimer.stop();
      getNextStep()->process(buf);
      return false;
    }

    void Predict::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }


  } //# end namespace
}
