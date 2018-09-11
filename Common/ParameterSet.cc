//# ParameterSet.cc: Implements a map of Key-Value pairs.
//#
//# Copyright (C) 2002-2003
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
//# $Id: ParameterSet.cc 23459 2013-01-08 08:49:52Z diepen $

//# Always #include <lofar_config.h> first!

#include "ParameterSet.h"
#include "ParameterRecord.h"

namespace DP3 {

//-------------------------- creation and destroy ---------------------------

ParameterSet* globalParameterSet()
{
  static ParameterSet ps;

  return &ps;
}

//-------------------------- creation and destroy ---------------------------
//
// Default constructor
//
ParameterSet::ParameterSet(KeyCompare::Mode	mode)
  : itsSet (new ParameterSetImpl(mode))
{}

ParameterSet::ParameterSet(bool caseInsensitive)
  : itsSet (new ParameterSetImpl(caseInsensitive ?
                                 KeyCompare::NOCASE : KeyCompare::NORMAL))
{}

//
// Construction by reading a parameter file.
//
ParameterSet::ParameterSet(const std::string& theFilename, bool caseInsensitive)
  : itsSet (new ParameterSetImpl(theFilename, caseInsensitive ?
                                 KeyCompare::NOCASE : KeyCompare::NORMAL))
{}

//
// Construction by reading a parameter file.
//
ParameterSet::ParameterSet(const std::string&	theFilename,
			   KeyCompare::Mode	mode)
  : itsSet (new ParameterSetImpl(theFilename, mode))
{}

ParameterSet::ParameterSet(const char*	theFilename,
			   KeyCompare::Mode	mode)
  : itsSet (new ParameterSetImpl(std::string(theFilename), mode))
{}

//
// Copying is allowed.
//
ParameterSet::ParameterSet(const ParameterSet& that)
  : itsSet (that.itsSet)
{}

//
// operator= copying
//
ParameterSet& 
ParameterSet::operator=(const ParameterSet& that)
{
  if (this != &that) {
    itsSet = that.itsSet;
  }
  return (*this);
}

//
//	Destructor
//
ParameterSet::~ParameterSet()
{}

ParameterRecord ParameterSet::getRecord (const std::string& aKey) const
{
  return get(aKey).getRecord();
}

//
// operator<<
//
std::ostream&	operator<< (std::ostream& os, const ParameterSet &thePS)
{
  os << *thePS.itsSet;
  return os;
}

} // namespace LOFAR
