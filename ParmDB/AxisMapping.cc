//# AxisMapping.cc: Map the cells of one axis to another
//#
//# Copyright (C) 2008
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
//# $Id: AxisMapping.cc 14038 2009-09-17 13:59:12Z diepen $

#include "AxisMapping.h"

namespace DP3 {
namespace BBS {

  // Map one axis onto another.
  AxisMapping::AxisMapping (const Axis& from, const Axis& to)
  {
    size_t nr = from.size();
    size_t maxto = to.size() - 1;
    itsMapping.reserve (nr);
    itsCenters.reserve (nr);
    itsBorders.reserve (nr);
    size_t inx=0;
    for (size_t i=0; i<nr; ++i) {
      size_t inxn = to.find (from.center(i), true, inx).first;
      if (inxn > maxto) inxn = maxto;
      if (inxn != inx  &&  i > 0) {
        itsBorders.push_back (i);
      }
      inx = inxn;
      itsMapping.push_back (inx);
      itsCenters.push_back ((from.center(i) - to.lower(inx)) / to.width(inx));
    }
    itsBorders.push_back (nr);
  }


  const AxisMapping& AxisMappingCache::makeMapping (const Axis& from,
                                                    const Axis& to)
  {
    std::pair<std::map<AxisKey,AxisMapping>::iterator, bool> result = 
      itsCache.insert (std::make_pair(AxisKey(from.getId(), to.getId()),
                                 AxisMapping(from, to)));
    // Make sure it got inserted.
    assert (result.second);
    // Return the new mapping.
    return result.first->second;
  }


  Location GridMapping::findLocation (AxisMappingCache& cache,
                                      const Location& location,
                                      const Grid& src,
                                      const Grid& dest)
  {
    const AxisMapping& mapx = cache.get (*src.getAxis(0), *dest.getAxis(0));
    const AxisMapping& mapy = cache.get (*src.getAxis(1), *dest.getAxis(1));
    return Location (mapx[location.first], mapy[location.second]);
  }

  Location GridMapping::findLocation (AxisMappingCache& cache,
                                      unsigned int cellId,
                                      const Grid& src,
                                      const Grid& dest)
  {
    return findLocation (cache, src.getCellLocation(cellId), src, dest);
  }

  unsigned int GridMapping::findCellId (AxisMappingCache& cache,
                                const Location& location,
                                const Grid& src,
                                const Grid& dest)
  {
    return dest.getCellId (findLocation (cache, location, src, dest));
  }

  unsigned int GridMapping::findcellId (AxisMappingCache& cache,
                                unsigned int cellId,
                                const Grid& src,
                                const Grid& dest)
  {
    return dest.getCellId (findLocation (cache, cellId, src, dest));
  }


} // namespace BBS
} // namespace LOFAR
