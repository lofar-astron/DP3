// AxisMapping.h: Map the cells of one axis to another
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Map the cells of one axis to another
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_AXISMAPPING_H
#define LOFAR_PARMDB_AXISMAPPING_H

#include "Grid.h"

#include <map>

namespace dp3 {
namespace parmdb {

/// @ingroup ParmDB
/// @brief Map the cells of one axis to another

/// @{
/// This class defines the mapping of one axis to another.
/// It is meant for mapping the grid axes of a predict to the axes
/// of the domain grid, so it has to be calculated only once per predict.
class AxisMapping {
 public:
  /// Create the mapping.
  AxisMapping(const Axis& from, const Axis& to);

  /// Iterator to get the next interval.
  ///@{
  typedef std::vector<int>::const_iterator const_iterator;
  const_iterator begin() const { return itsMapping.begin(); }
  const_iterator end() const { return itsMapping.end(); }
  ///@}

  /// Get the number of elements.
  size_t size() const { return itsMapping.size(); }

  /// Get the to interval for the i-th from interval.
  int operator[](int i) const { return itsMapping[i]; }

  /// Get a pointer to the scaled center values.
  /// The center value of each interval in the from axis is scaled for
  /// its interval in the to axis.
  const double* getScaledCenters() const { return &(itsCenters[0]); }

  /// Get the borders telling which from-cells map to the same to-cell.
  /// For example: borders (5,8,12) mean that from-cells 0-4 map to the same
  /// to-cell, and so do from-cells 5-7 and 8-11.
  const std::vector<int>& getBorders() const { return itsBorders; }

 private:
  std::vector<int> itsMapping;     ///< cellnr in to for from-cell i
  std::vector<double> itsCenters;  ///< center of from-cell i scaled to to-cell
  std::vector<int> itsBorders;     ///< last from-cell mapped to same to-cell i
};

/// @brief This class caches Grid objects.
///
/// It is used to achieve that parameters
/// with equal domains share the same Axis objects. In this way
/// the AxisMapping objects can also be shared which can improve the
/// performance.
//
/// The cache consists of two parts:
/// <ol>
///  <li> A map of axis ID to Axis objects is used to know which Axis
/// objects are available and to find them by Axis ID.
///  <li> Another map is used to find the possible pair of Axis objects
/// that form a Grid with a given hash value. Note that different Grids
/// might map to the same hash value, although in practice that will
/// hardly ever occur.
/// </ol>
class AxisCache {
 public:
  /// Get a Grid object for the given domains.
  /// It adds Axis objects to the cache if they are new.
  Grid getGrid(const std::vector<Box>& domains);

  /// Test if an Axis object equal to the given one occurs in the cache.
  /// If so, return that object. Otherwise add the Axis to the cache.
  Axis::ShPtr getAxis(const Axis& axis);

  /// Clear the cache.
  void clear() {
    itsGridCache.clear();
    itsAxisCache.clear();
  }

 private:
  std::map<int64_t, std::vector<std::pair<int, int>>> itsGridCache;
  std::map<int, Axis::ShPtr> itsAxisCache;
};

/// @brief This class caches axis mappings. It uses the unique id of the
/// from-axis and to-axis as the key in the cache. <br>It is meant to avoid
/// creating many identical mappings, because the same axes are used for many
/// parameters.
class AxisMappingCache {
 private:
  /// Define the key consisting of both id-s.
  struct AxisKey {
    AxisKey(unsigned int fromId, unsigned int toId)
        : itsFrom(fromId), itsTo(toId) {}
    unsigned int itsFrom;
    unsigned int itsTo;

    bool operator==(const AxisKey that) const {
      return itsFrom == that.itsFrom && itsTo == that.itsTo;
    }
    bool operator!=(const AxisKey that) const {
      return itsFrom != that.itsFrom || itsTo != that.itsTo;
    }
    bool operator<(const AxisKey that) const {
      return itsFrom < that.itsFrom ||
             (itsFrom == that.itsFrom && itsTo < that.itsTo);
    }
  };

 public:
  /// Get the number of elements.
  size_t size() const { return itsCache.size(); }

  /// Clear the cache.
  void clear() { itsCache.clear(); }

  /// Find the possible mapping of axis 'from' to axis 'to'.
  /// If no existing, create it.
  const AxisMapping& get(const Axis& from, const Axis& to) {
    std::map<AxisKey, AxisMapping>::const_iterator iter =
        itsCache.find(AxisKey(from.getId(), to.getId()));
    return (iter == itsCache.end() ? makeMapping(from, to) : iter->second);
  }

 private:
  /// Create a mapping of axis from to axis to and add it to the cache.
  const AxisMapping& makeMapping(const Axis& from, const Axis& to);

  std::map<AxisKey, AxisMapping> itsCache;
};

class GridMapping {
 public:
  /// Find the location in grid 'dest', given the location in grid 'src'.
  static Grid::Location findLocation(AxisMappingCache& cache,
                                     const Grid::Location& location,
                                     const Grid& src, const Grid& dest);

  /// Find the location in grid 'dest', given the cellId in grid 'src'.
  static Grid::Location findLocation(AxisMappingCache& cache,
                                     unsigned int cellId, const Grid& src,
                                     const Grid& dest);

  /// Find the cellId in grid 'dest', given the location in grid 'src'.
  static unsigned int findCellId(AxisMappingCache& cache,
                                 const Grid::Location& location,
                                 const Grid& src, const Grid& dest);

  /// Find the cellId in grid 'dest', given the cellId in grid 'src'.
  static unsigned int findcellId(AxisMappingCache& cache, unsigned int cellId,
                                 const Grid& src, const Grid& dest);
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
