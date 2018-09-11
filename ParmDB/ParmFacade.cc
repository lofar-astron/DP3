//# ParmFacade.cc: Object access the parameter database
//#
//# Copyright (C) 2006
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
//# $Id: ParmFacade.cc 27639 2013-12-04 08:02:12Z diepen $

#include "ParmFacade.h"
#include "ParmFacadeLocal.h"
#include "ParmFacadeDistr.h"

#include <casacore/tables/Tables/Table.h>

using namespace std;
using namespace casacore;

// Create tParmFacade.in_mep with parmdb using:
//   create tablename='tParmFacade.in_mep'
//   add parm1 domain=[1,5,4,10],values=2
//   add parm2 domain=[1,5,4,10],values=[2,0.1],nx=2
//   add parm3 type='expression',expression='parm1*parm2'

namespace DP3 {
  namespace BBS {

    ParmFacade::ParmFacade (const string& tableName, bool create)
    {
      // If create, only a local ParmDB can be done.
      // If it is an existing table, open it directly.
      // Otherwise it is a distributed ParmDB.
      if (create) {
        itsRep = ParmFacadeRep::ShPtr(new ParmFacadeLocal(tableName, create));
      } else if (Table::isReadable(tableName)) {
        itsRep = ParmFacadeRep::ShPtr(new ParmFacadeLocal(tableName));
      } else {
        //itsRep = ParmFacadeRep::ShPtr(new ParmFacadeDistr(tableName));
				throw std::runtime_error("distributed parm facade not available");
      }
    }

    ParmFacade::~ParmFacade()
    {}

    Record ParmFacade::getValues (const string& parmNamePattern,
                                  double freqv1, double freqv2,
                                  double timev1, double timev2,
                                  bool asStartEnd,
                                  bool includeDefaults)
    {
      double sfreq = freqv1;
      double efreq = freqv2;
      double stime = timev1;
      double etime = timev2;
      if (!asStartEnd) {
        sfreq = freqv1 - freqv2/2;
        efreq = sfreq  + freqv2;
        stime = timev1 - timev2/2;
        etime = stime  + timev2;
      }
      vector<double> rng = getRange (parmNamePattern);
      // No values if the range is null.
      if (rng[0] == 0  &&  rng[1] == 0) {
        return Record();
      }
      if (sfreq < rng[0]) sfreq = rng[0];
      if (efreq > rng[1]) efreq = rng[1];
      if (stime < rng[2]) stime = rng[2];
      if (etime > rng[3]) etime = rng[3];
      return itsRep->getValues (parmNamePattern, sfreq, efreq, 0,
                                stime, etime, 0, true, includeDefaults);
    }

    // Get the parameter values for the given parameters and domain.
    map<string, vector<double> >
    ParmFacade::getValuesMap (const string& parmNamePattern,
                              double freqv1, double freqv2, double freqStep,
                              double timev1, double timev2, double timeStep,
                              bool asStartEnd, bool includeDefaults)
    {
      return record2Map (getValues (parmNamePattern,
                                    freqv1, freqv2, freqStep,
                                    timev1, timev2, timeStep,
                                    asStartEnd, includeDefaults));
    }

    map<string,vector<double> >
    ParmFacade::record2Map (const Record& rec) const
    {
      map<string, vector<double> > out;
      // Copy all values from the record to the map.
      for (uint i=0; i<rec.nfields(); ++i) {
        const String& name = rec.name(i);
        Record subrec = rec.subRecord(i);
        // First make empty vector; thereafter copy values to it.
        vector<double>& vec = out[name];
        assert (vec.size() == 0);
        // Get result and put in vector.
        const Array<double>& arr = subrec.asArrayDouble("values");
        vec.assign (arr.begin(), arr.end());
      }
      return out;
    }

  } // namespace ParmDB
} // namespace LOFAR
