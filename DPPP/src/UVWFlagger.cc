//# UVWFlagger.cc: DPPP step class to flag data on channel, baseline, or time
//# Copyright (C) 2010
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
#include <DPPP/UVWFlagger.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/AverageInfo.h>
#include <Common/ParameterSet.h>
#include <Common/StreamUtil.h>
#include <Common/LofarLogger.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Utilities/GenSort.h>
#include <iostream>
#include <algorithm>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    UVWFlagger::UVWFlagger (DPInput* input,
                            const ParameterSet& parset, const string& prefix)
      : itsInput  (input),
        itsName   (prefix),
        itsNTimes (0)
    {
      itsRangeUVm = fillUVW (parset, prefix, "uvm", true);
      itsRangeUm  = fillUVW (parset, prefix, "um", false);
      itsRangeVm  = fillUVW (parset, prefix, "vm", false);
      itsRangeWm  = fillUVW (parset, prefix, "wm", false);
      itsRangeUVl = fillUVW (parset, prefix, "uvlambda", true);
      itsRangeUl  = fillUVW (parset, prefix, "ulambda", false);
      itsRangeVl  = fillUVW (parset, prefix, "vlambda", false);
      itsRangeWl  = fillUVW (parset, prefix, "wlambda", false);
      ASSERTSTR (itsRangeUVm.size() + itsRangeUVl.size() +
                 itsRangeUm.size() + itsRangeVm.size() + itsRangeWm.size() +
                 itsRangeUl.size() + itsRangeVl.size() + itsRangeWl.size() > 0,
                 "One or more ranges in UVWFlagger has to be filled in");
      itsCenter = parset.getStringVector (prefix+"phasecenter",
                                          vector<string>());
      if (! itsCenter.empty()) {
        handleCenter();
      }
    }

    UVWFlagger::~UVWFlagger()
    {}

    void UVWFlagger::show (std::ostream& os) const
    {
      os << "UVWFlagger " << itsName << std::endl;
    }

    void UVWFlagger::showCounts (std::ostream& os) const
    {
      os << endl << "Flag statistics of UVWFlagger " << itsName;
      os << endl << "=============================" << endl;
      itsFlagCounter.showBaseline (os, itsInput->getAnt1(),
                                   itsInput->getAnt2(), itsNTimes);
      itsFlagCounter.showChannel  (os, itsNTimes);
    }

    void UVWFlagger::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " UVWFlagger " << itsName << endl;
    }

    void UVWFlagger::updateAverageInfo (AverageInfo& info)
    {
      // Convert the given frequencies to possibly averaged frequencies.
      averageFreqs (info.startChan(), info.nchanAvg());
    }

    bool UVWFlagger::process (const DPBuffer& buf)
    {
      DPBuffer out(buf);
      // The flags will be changed, so make sure we have a unique array.
      out.getFlags().unique();
      // Let the next step do its processing.
      itsNTimes++;
      getNextStep()->process (out);
      return true;
    }

    void UVWFlagger::finish()
    {
      // Let the next step finish its processing.
      getNextStep()->finish();
    }

    void UVWFlagger::flagUV (const Matrix<double>& uvw,
                             Cube<bool>& flags)
    {
      const IPosition& shape = flags.shape();
      uint nr = shape[0] * shape[1];
      uint nrbl = shape[2];
      const double* uvwPtr = uvw.data();
      bool* flagPtr = flags.data();
      for (uint i=0; i<nrbl; ++i) {
        // UV-distance is sqrt(u^2 + v^2).
        // The sqrt is not needed because minuv and maxuv are squared.
        double uvdist = uvwPtr[0] * uvwPtr[0] + uvwPtr[1] * uvwPtr[1];
        uvwPtr  += 3;
        flagPtr += nr;
      }
    }

    vector<double> UVWFlagger::fillUVW (const ParameterSet& parset,
                                        const string& prefix,
                                        const string& name,
                                        bool square)
    {
      // Get possible range, minimum, and maximum.
      vector<string> uvs = parset.getStringVector (prefix + name + "range",
                                                   vector<string>());
      double minuv = parset.getDouble (prefix + name + "min", 0.);
      double maxuv = parset.getDouble (prefix + name + "max", 0.);
      // Process the ranges.
      vector<double> vals;
      vals.reserve (2*uvs.size());
      for (vector<string>::const_iterator str = uvs.begin();
           str != uvs.end(); ++str) {
        // Each range can be given as st..end or val+-halfwidth.
        // Find the .. or +- token.
        bool usepm = false;
        string::size_type pos;
        pos = str->find ("..");
        if (pos == string::npos) {
          usepm = true;
          pos = str->find ("+-");
          ASSERTSTR (pos != string::npos, "UVWFlagger " << name << "range '"
                     << *str << "' should be range using .. or +-");
        }
        string str1 = str->substr (0, pos);
        string str2 = str->substr (pos+2);
        vals.push_back (strToDouble(str1));
        vals.push_back (strToDouble(str2));
      }
      // If minimum or maximum is given, add them as a range as well.
      if (minuv > 0) {
        vals.push_back (0.);
        vals.push_back (minuv);
      }
      if (maxuv > 0) {
        vals.push_back (maxuv);
        vals.push_back (1e15);
      }
      if (square) {
        for (vector<double>::iterator iter = vals.begin();
             iter != vals.end(); ++iter) {
          *iter = *iter * *iter;
        }
      }
      return vals;
    }

    void UVWFlagger::handleCenter()
    {
      // The phase center can be given as one, two, or three values.
      // I.e., as source name, ra,dec or ra,dec,frame.
      ASSERTSTR (itsCenter.size() < 4,
                 "Up to 3 values can be given in UVWFlagger phasecenter");
      MDirection phaseCenter;
      if (itsCenter.size() == 1) {
        string type = toLower(itsCenter[0]);
        if (type == "sun") {
          phaseCenter = MDirection(MDirection::SUN);
        } else {
          ASSERTSTR (false, itsCenter[0] << " is an unknown source");
        }
      } else {
        Quantity q0, q1;
        ASSERTSTR (MVAngle::read (q0, itsCenter[0]),
                   itsCenter[0] << " is an invalid RA or longitude");
        ASSERTSTR (MVAngle::read (q1, itsCenter[1]),
                   itsCenter[1] << " is an invalid DEC or latitude");
        MDirection::Types type = MDirection::J2000;
        if (itsCenter.size() > 2) {
        }
        phaseCenter = MDirection(MVDirection(q0, q1), type);
      }
      // Create the UVW calculator.
      itsUVWCalc = UVWCalculator (phaseCenter, itsInput->antennaPos());
    }

    void UVWFlagger::averageFreqs (uint startChan, uint nchanAvg)
    {
      uint nchan = (itsFreqs.size() - startChan) / nchanAvg;
      Vector<double> freqs(nchan); 
      for (uint i=0; i<nchan; ++i) {
        freqs[i] = 0.5 * (itsFreqs[startChan + i*nchanAvg] +
                          itsFreqs[startChan + (i+1)*nchanAvg-1]);
      }
      // Replace the original freqs by the averaged ones.
      itsFreqs.reference (freqs);
    }

  } //# end namespace
}
