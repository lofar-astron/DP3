//# ParSet.cc: Wrapper around ParaMeterSet to keep track of parameters asked for
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

#include <DPPP/ParSet.h>
#include <Common/ParameterRecord.h>
#include <set>

namespace LOFAR {
  namespace DPPP {

    ParSet::ParSet (const ParameterSet& parset )
      : itsParSet (parset)
    {}

    bool   ParSet::getBool  (const string& aKey) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getBool (aKey);
    }

    bool   ParSet::getBool  (const string& aKey, bool aValue) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getBool (aKey, aValue);
    }

    int    ParSet::getInt   (const string& aKey) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getInt (aKey);
    }

    int    ParSet::getInt   (const string& aKey, int aValue) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getInt (aKey, aValue);
    }

    uint   ParSet::getUint  (const string& aKey) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getUint (aKey);
    }

    uint   ParSet::getUint  (const string& aKey, uint aValue) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getUint (aKey, aValue);
    }

    float  ParSet::getFloat (const string& aKey) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getFloat (aKey);
    }

    float  ParSet::getFloat (const string& aKey, float aValue) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getFloat (aKey, aValue);
    }

    double ParSet::getDouble(const string& aKey) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getDouble (aKey);
    }

    double ParSet::getDouble(const string& aKey, double aValue) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getDouble (aKey, aValue);
    }

    string ParSet::getString(const string& aKey) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getString (aKey);
    }

    string ParSet::getString(const string& aKey, const string& aValue) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getString (aKey, aValue);
    }

    vector<bool>   ParSet::getBoolVector  (const string& aKey,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getBoolVector (aKey, expandable);
    }

    vector<bool>   ParSet::getBoolVector  (const string& aKey,
                                           const vector<bool>& aValue,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getBoolVector (aKey, aValue, expandable);
    }

    vector<int>    ParSet::getIntVector   (const string& aKey,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getIntVector (aKey, expandable);
    }

    vector<int>    ParSet::getIntVector   (const string& aKey,
                                           const vector<int>& aValue,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getIntVector (aKey, aValue, expandable);
    }

    vector<uint>   ParSet::getUintVector  (const string& aKey,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getUintVector (aKey, expandable);
    }

    vector<uint>   ParSet::getUintVector  (const string& aKey,
                                           const vector<uint>& aValue,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getUintVector (aKey, aValue, expandable);
    }

    vector<float>  ParSet::getFloatVector (const string& aKey,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getFloatVector (aKey, expandable);
    }

    vector<float>  ParSet::getFloatVector (const string& aKey,
                                           const vector<float>& aValue,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getFloatVector (aKey, aValue, expandable);
    }

    vector<double> ParSet::getDoubleVector(const string& aKey,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getDoubleVector (aKey, expandable);
    }

    vector<double> ParSet::getDoubleVector(const string& aKey,
                                           const vector<double>& aValue,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getDoubleVector (aKey, aValue, expandable);
    }

    vector<string> ParSet::getStringVector(const string& aKey,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getStringVector (aKey, expandable);
    }

    vector<string> ParSet::getStringVector(const string& aKey,
                                           const vector<string>& aValue,
                                           bool expandable) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getStringVector (aKey, aValue, expandable);
    }

    ParameterRecord ParSet::getRecord(const string& aKey) const
    {
      itsAskedParms.insert (aKey);
      return itsParSet.getRecord (aKey);
    }

    vector<string> ParSet::unusedKeys() const
    {
      vector<string> vec;
      for (ParameterSet::const_iterator iter = itsParSet.begin();
           iter != itsParSet.end(); ++iter) {
        if (itsAskedParms.find (iter->first) == itsAskedParms.end()) {
          vec.push_back (iter->first);
        }
      }
      return vec;
    }


  } //# end namespace
}
