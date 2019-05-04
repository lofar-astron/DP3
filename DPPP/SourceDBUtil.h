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

#include "Patch.h"

#include <vector>

namespace DP3
{
namespace BBS
{
class SourceDB;
}

namespace DPPP
{

// \addtogroup NDPPP
// @{
  std::vector<Patch::ConstPtr> makePatches(BBS::SourceDB &sourceDB,
                                      const std::vector<std::string> &patchNames,
                                      unsigned int nModel);

  // Create a source list (with patch name) from a patchlist
  // Needed for efficient multithreading
  std::vector<std::pair<ModelComponent::ConstPtr,Patch::ConstPtr> >
  makeSourceList (const std::vector<Patch::ConstPtr>& patchList);

  // From a given PatchList, create a new one with one patch per component
  std::vector<Patch::ConstPtr> makeOnePatchPerComponent(
      const std::vector<Patch::ConstPtr>&);

  std::vector<std::string>  makePatchList(BBS::SourceDB &sourceDB,
                                std::vector<std::string> patterns);

  bool checkPolarized(BBS::SourceDB &sourceDB,
                      const std::vector<std::string> &patchNames,
                      unsigned int nModel);

// @}

} //# namespace DPPP
} //# namespace LOFAR

#endif
