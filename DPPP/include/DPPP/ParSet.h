//# ParSet.h: Wrapper around ParaMeterSet to keep track of parameters asked for
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

#ifndef DPPP_PARSET_H
#define DPPP_PARSET_H

// @file
// @brief Wrapper around ParaMeterSet to keep track of parameters asked for

#include <Common/ParameterSet.h>
#include <Common/lofar_string.h>
#include <Common/lofar_vector.h>
#include <Common/lofar_set.h>

namespace LOFAR {

  namespace DPPP {

    // @ingroup NDPPP

    // This class keeps track of the parameters asked for.
    // In this way it is possible to give a warning or error that the
    // ParameterSet has more parameters than asked for. It is an aid in
    // detecting misspelled parameters names.

    class ParSet
    {
    public:
      ParSet (const ParameterSet&);

      // Ask a parameter with a possible default value.
      // @{
      bool   getBool  (const string& aKey) const;
      bool   getBool  (const string& aKey, bool aValue) const;
      int    getInt   (const string& aKey) const;
      int    getInt   (const string& aKey, int aValue) const;
      uint   getUint  (const string& aKey) const;
      uint   getUint  (const string& aKey, uint aValue) const;
      float  getFloat (const string& aKey) const;
      float  getFloat (const string& aKey, float aValue) const;
      double getDouble(const string& aKey) const;
      double getDouble(const string& aKey, double aValue) const;
      string getString(const string& aKey) const;
      string getString(const string& aKey, const string& aValue) const;
      vector<bool>   getBoolVector  (const string& aKey,
                                     bool expandable = false) const;
      vector<bool>   getBoolVector  (const string& aKey,
                                     const vector<bool>& aValue,
                                     bool expandable = false) const;
      vector<int>    getIntVector   (const string& aKey,
                                     bool expandable = false) const;
      vector<int>    getIntVector   (const string& aKey,
                                     const vector<int>& aValue,
                                     bool expandable = false) const;
      vector<uint>   getUintVector  (const string& aKey,
                                     bool expandable = false) const;
      vector<uint>   getUintVector  (const string& aKey,
                                     const vector<uint>& aValue,
                                     bool expandable = false) const;
      vector<float>  getFloatVector (const string& aKey,
                                     bool expandable = false) const;
      vector<float>  getFloatVector (const string& aKey,
                                     const vector<float>& aValue,
                                     bool expandable = false) const;
      vector<double> getDoubleVector(const string& aKey,
                                     bool expandable = false) const;
      vector<double> getDoubleVector(const string& aKey,
                                     const vector<double>& aValue,
                                     bool expandable = false) const;
      vector<string> getStringVector(const string& aKey,
                                     bool expandable = false) const;
      vector<string> getStringVector(const string& aKey,
                                     const vector<string>& aValue,
                                     bool expandable = false) const;
      ParameterRecord getRecord (const string& aKey) const;

      // Get all unused parameters.
      vector<string> unusedKeys() const;

      // Get the underlying ParameterSet.
      const ParameterSet& parameterSet() const
        { return itsParSet; }

    private:
      ParameterSet        itsParSet;
      mutable set<string> itsAskedParms;
    };

  } //# end namespace
}

#endif
