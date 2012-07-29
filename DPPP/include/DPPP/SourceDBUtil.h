//# SourceDBUtil.h: Helper functions to extract patch and source information
//# from a SourceDB.
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
//# $Id$

#ifndef DPPP_SOURCEDBUTIL_H
#define DPPP_SOURCEDBUTIL_H

// \file
// Helper functions to extract patch and source information from a SourceDB.

#include <DPPP/Patch.h>
#include <Common/lofar_string.h>
#include <Common/LofarTypes.h>

namespace LOFAR
{
namespace BBS
{
class SourceDB;
}

namespace DPPP
{

// \addtogroup NDPPP
// @{
  vector<Patch::ConstPtr> makePatches(BBS::SourceDB &sourceDB,
                                      const vector<string> &patchNames,
                                      uint nModel);
// @}

} //# namespace DPPP
} //# namespace LOFAR

#endif
