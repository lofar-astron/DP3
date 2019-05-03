//# ParmDBBlob.cc: Dummy class to hold parmaeter values
//#
//# Copyright (C) 2012
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
//# $Id: ParmDBBlob.cc 21598 2012-07-16 08:07:34Z diepen $

#include "ParmDBBlob.h"

using namespace std;

namespace DP3 {
namespace BBS {

  ParmDBBlob::ParmDBBlob (const string&, bool)
  {}

  ParmDBBlob::~ParmDBBlob()
  {}

  void ParmDBBlob::flush (bool)
  {}

  void ParmDBBlob::lock (bool)
  {}

  void ParmDBBlob::unlock()
  {}

  void ParmDBBlob::clearTables()
  {}

  void ParmDBBlob::setDefaultSteps (const vector<double>&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  int ParmDBBlob::getNameId (const std::string&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  Box ParmDBBlob::getRange (const string&) const
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  Box ParmDBBlob::getRange (const std::vector<std::string>&) const
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  void ParmDBBlob::getValues (vector<ParmValueSet>&,
                              const vector<unsigned int>&,
                              const vector<ParmId>&,
                              const Box&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  void ParmDBBlob::getDefValues (ParmMap&,
                                 const string&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  void ParmDBBlob::putValues (const string&, int&, ParmValueSet&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  void ParmDBBlob::putDefValue (const string&, const ParmValueSet&,
                                bool)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  void ParmDBBlob::deleteValues (const string&,
                                 const Box&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  void ParmDBBlob::deleteDefValues (const string&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  vector<string> ParmDBBlob::getNames (const string&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

  void ParmDBBlob::fillDefMap (ParmMap&)
  {
    throw std::runtime_error("ParmDBBlob not implemented");
  }

} // namespace BBS
} // namespace LOFAR
