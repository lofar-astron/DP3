// AxisMapping.cc: Map the cells of one axis to another
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "AxisMapping.h"

namespace dp3 {
namespace parmdb {

// Map one axis onto another.
AxisMapping::AxisMapping(const Axis& from, const Axis& to) {
  size_t nr = from.size();
  size_t maxto = to.size() - 1;
  itsMapping.reserve(nr);
  itsCenters.reserve(nr);
  itsBorders.reserve(nr);
  size_t inx = 0;
  for (size_t i = 0; i < nr; ++i) {
    size_t inxn = to.find(from.center(i), true, inx).first;
    if (inxn > maxto) inxn = maxto;
    if (inxn != inx && i > 0) {
      itsBorders.push_back(i);
    }
    inx = inxn;
    itsMapping.push_back(inx);
    itsCenters.push_back((from.center(i) - to.lower(inx)) / to.width(inx));
  }
  itsBorders.push_back(nr);
}

const AxisMapping& AxisMappingCache::makeMapping(const Axis& from,
                                                 const Axis& to) {
  std::pair<std::map<AxisKey, AxisMapping>::iterator, bool> result =
      itsCache.insert(std::make_pair(AxisKey(from.getId(), to.getId()),
                                     AxisMapping(from, to)));
  // Make sure it got inserted.
  assert(result.second);
  // Return the new mapping.
  return result.first->second;
}

Grid::Location GridMapping::findLocation(AxisMappingCache& cache,
                                         const Grid::Location& location,
                                         const Grid& src, const Grid& dest) {
  const AxisMapping& mapx = cache.get(*src.getAxis(0), *dest.getAxis(0));
  const AxisMapping& mapy = cache.get(*src.getAxis(1), *dest.getAxis(1));
  return Grid::Location(mapx[location.first], mapy[location.second]);
}

Grid::Location GridMapping::findLocation(AxisMappingCache& cache,
                                         unsigned int cellId, const Grid& src,
                                         const Grid& dest) {
  return findLocation(cache, src.getCellLocation(cellId), src, dest);
}

unsigned int GridMapping::findCellId(AxisMappingCache& cache,
                                     const Grid::Location& location,
                                     const Grid& src, const Grid& dest) {
  return dest.getCellId(findLocation(cache, location, src, dest));
}

unsigned int GridMapping::findcellId(AxisMappingCache& cache,
                                     unsigned int cellId, const Grid& src,
                                     const Grid& dest) {
  return dest.getCellId(findLocation(cache, cellId, src, dest));
}

}  // namespace parmdb
}  // namespace dp3
