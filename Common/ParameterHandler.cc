//#  ParameterHandler.cc: 
//#
//#  Copyright (C) 2007
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
//#  $Id: ParameterHandler.cc 16886 2010-12-08 10:43:17Z diepen $

#include "ParameterHandler.h"

using namespace std;

namespace DP3 { namespace CEP {

  ParameterHandler::ParameterHandler (const ParameterSet& parSet)
    : itsParms (parSet)
  {}

  string ParameterHandler::getString (const string& parm,
				      const string& defVal) const
  {
    if (itsParms.isDefined(parm)) {
      return itsParms.getString (parm);
    }
    return defVal;
  }

  double ParameterHandler::getDouble (const string& parm,
				      double defVal) const
  {
    if (itsParms.isDefined(parm)) {
      return itsParms.getDouble (parm);
    }
    return defVal;
  }

  unsigned ParameterHandler::getUint (const string& parm,
				      unsigned defVal) const
  {
    if (itsParms.isDefined(parm)) {
      return itsParms.getUint32 (parm);
    }
    return defVal;
  }

  bool ParameterHandler::getBool (const string& parm,
				  bool defVal) const
  {
    if (itsParms.isDefined(parm)) {
      return itsParms.getBool (parm);
    }
    return defVal;
  }

  vector<string> ParameterHandler::getStringVector
  (const string& parm, const vector<string>& defVal) const
  {
    if (itsParms.isDefined(parm)) {
      return itsParms.getStringVector (parm);
    }
    return defVal;
  }

  void ParameterHandler::fillString (const string& parm,
				     string& value) const
  {
    if (itsParms.isDefined(parm)) {
      value = itsParms.getString (parm);
    }
  }

  void ParameterHandler::fillDouble (const string& parm,
				     double& value) const
  {
    if (itsParms.isDefined(parm)) {
      value = itsParms.getDouble (parm);
    }
  }

  void ParameterHandler::fillUint (const string& parm,
				   unsigned& value) const
  {
    if (itsParms.isDefined(parm)) {
      value = itsParms.getUint32 (parm);
    }
  }

  void ParameterHandler::fillBool (const string& parm,
				   bool& value) const
  {
    if (itsParms.isDefined(parm)) {
      value = itsParms.getBool (parm);
    }
  }

  void ParameterHandler::fillStringVector (const string& parm,
					   vector<string>& value) const
  {
    if (itsParms.isDefined(parm)) {
      value = itsParms.getStringVector (parm);
    }
  }


  DP3::BlobOStream& operator<< (DP3::BlobOStream& bs, const ParameterSet& m)
  {
    bs.putStart ("ParameterSet", 1);
    bs << static_cast<uint32_t>(m.size());
    for (ParameterSet::const_iterator it=m.begin();
         it!=m.end();
         ++it) {
      bs << it->first << it->second.get();
    }
    bs.putEnd();
    return bs;
  }

  DP3::BlobIStream& operator>> (DP3::BlobIStream& bs, ParameterSet& m)
  {
    bs.getStart ("ParameterSet");
    m.clear();
    uint32_t size;
    bs >> size;
    std::string k,v;
    for (uint32_t i=0; i<size; ++i) {
      bs >> k >> v;
      m.add (k, v);
    }
    bs.getEnd();
    return bs;
  }

}} // end namespaces
