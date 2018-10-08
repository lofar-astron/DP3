//# ApplyBeam.cc: DPPP step class to ApplyBeam visibilities
//# Copyright (C) 2015
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

#include <iostream>

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"
#include "../Common/StringUtil.h"

#include "ApplyBeam.h"
#include "DPInfo.h"
#include "Exceptions.h"
#include "FlagCounter.h"
#include "Position.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/algorithm/string/case_conv.hpp>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    ApplyBeam::ApplyBeam(DPInput* input, const ParameterSet& parset,
                         const string& prefix, bool substep)
        :
          itsInput(input),
          itsName(prefix),
          itsUpdateWeights(parset.getBool(prefix + "updateweights", false)),
          itsUseChannelFreq(parset.getBool(prefix + "usechannelfreq", true)),
          itsDebugLevel(parset.getInt(prefix + "debuglevel", 0))
    {
      // only read 'invert' parset key if it is a separate step
      // if applybeam is called from gaincal/predict, the invert key should always be false
      if (substep) {
        itsInvert=false;
      } else {
        itsInvert=parset.getBool(prefix + "invert", true);
      }
      string mode=boost::to_lower_copy(parset.getString(prefix + "beammode","default"));
      assert (mode=="default" || mode=="array_factor" || mode=="element");
      if (mode=="default") {
        itsMode=DEFAULT;
      } else if (mode=="array_factor") {
        itsMode=ARRAY_FACTOR;
      } else if (mode=="element") {
        itsMode=ELEMENT;
      } else {
        throw Exception("Beammode should be DEFAULT, ARRAY_FACTOR or ELEMENT");
      }
    }

    ApplyBeam::ApplyBeam()
    {
    }

    ApplyBeam::~ApplyBeam()
    {
    }

    void ApplyBeam::updateInfo(const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();
      if (itsUpdateWeights) {
        info().setWriteWeights();
      }

      MDirection dirJ2000(
          MDirection::Convert(infoIn.phaseCenter(), MDirection::J2000)());
      Quantum<Vector<Double> > angles = dirJ2000.getAngle();
      itsPhaseRef = Position(angles.getBaseValue()[0],
                             angles.getBaseValue()[1]);

      const size_t nSt = info().nantenna();
      const size_t nCh = info().nchan();

      itsBeamValues.resize(NThreads());

      // Create the Measure ITRF conversion info given the array position.
      // The time and direction are filled in later.
      itsMeasConverters.resize(NThreads());
      itsMeasFrames.resize(NThreads());
      itsAntBeamInfo.resize(NThreads());

      for (size_t thread = 0; thread < NThreads(); ++thread) {
        itsBeamValues[thread].resize(nSt * nCh);
        itsMeasFrames[thread].set(info().arrayPosCopy());
        itsMeasFrames[thread].set(
            MEpoch(MVEpoch(info().startTime() / 86400), MEpoch::UTC));
        itsMeasConverters[thread].set(
            MDirection::J2000,
            MDirection::Ref(MDirection::ITRF, itsMeasFrames[thread]));
        itsInput->fillBeamInfo(itsAntBeamInfo[thread], info().antennaNames());
      }
    }

    void ApplyBeam::show(std::ostream& os) const
    {
      os << "ApplyBeam " << itsName << endl;
      os << "  mode:              ";
      if (itsMode==DEFAULT)
        os<<"default";
      else if (itsMode==ARRAY_FACTOR)
        os<<"array_factor";
      else os<<"element";
      os << endl;
      os << "  use channelfreq:   " << boolalpha << itsUseChannelFreq << endl;
      os << "  invert:            " << boolalpha << itsInvert << endl;
      os << "  update weights:    " << boolalpha << itsUpdateWeights << endl;
    }

    void ApplyBeam::showTimings(std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
      os << " ApplyBeam " << itsName << endl;
    }

    bool ApplyBeam::processMultithreaded(const DPBuffer& bufin, size_t thread)
    {
      itsTimer.start();
      itsBuffer.copy (bufin);
      Complex* data=itsBuffer.getData().data();

      if (itsUpdateWeights) {
        itsInput->fetchWeights (bufin, itsBuffer, itsTimer);
      }
      float* weight = itsBuffer.getWeights().data();

      double time = itsBuffer.getTime();

      //Set up directions for beam evaluation
      LOFAR::StationResponse::vector3r_t refdir, tiledir;

      /**
       * I'm not sure this is correct the way it is. These loops
       * seem to initialize variables that are never used in a
       * multi-threaded way, and if they were used from multiple
       * threads, it would imply process() is called multiple times,
       * and hence this initialization is already subject to a race
       * condition... ???
       * itsMeasFrames seems not to be actually used.
       * André, 2018-10-07
       */
      for (size_t threadIter = 0; threadIter < NThreads(); ++threadIter) {
        itsMeasFrames[threadIter].resetEpoch(
            MEpoch(MVEpoch(time / 86400), MEpoch::UTC));
        //Do a conversion on all threads, because converters are not
        //thread safe and apparently need to be used at least once
        refdir = dir2Itrf(info().delayCenter(), itsMeasConverters[threadIter]);
        tiledir = dir2Itrf(info().tileBeamDir(), itsMeasConverters[threadIter]);
      }

      LOFAR::StationResponse::vector3r_t srcdir = refdir;
      applyBeam(info(), time, data, weight, srcdir, refdir, tiledir,
                itsAntBeamInfo[thread], itsBeamValues[thread],
                itsUseChannelFreq, itsInvert, itsMode, itsUpdateWeights);

      itsTimer.stop();
      getNextStep()->process(itsBuffer);
      return false;
    }

    LOFAR::StationResponse::vector3r_t ApplyBeam::dir2Itrf(
        const MDirection& dir, MDirection::Convert& measConverter)
    {
      const MDirection& itrfDir = measConverter(dir);
      const Vector<Double>& itrf = itrfDir.getValue().getValue();
      LOFAR::StationResponse::vector3r_t vec;
      vec[0] = itrf[0];
      vec[1] = itrf[1];
      vec[2] = itrf[2];
      return vec;
    }

    void ApplyBeam::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } //# end namespace
}
