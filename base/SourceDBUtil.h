// SourceDBUtil.h: Helper functions to extract patch and source information
// from a SourceDB.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// \file
/// Helper functions to extract patch and source information from a SourceDB.

#ifndef DPPP_SOURCEDBUTIL_H
#define DPPP_SOURCEDBUTIL_H

#include "Patch.h"
#include "../parmdb/SourceDB.h"
#include "../parmdb/SourceDBSkymodel.h"

#include <string>
#include <variant>
#include <vector>

namespace dp3 {
namespace base {

std::vector<std::shared_ptr<Patch>> makePatches(
    parmdb::SourceDB &sourceDB, const std::vector<std::string> &patchNames,
    unsigned int nModel);

std::vector<std::shared_ptr<Patch>> MakePatches(
    const parmdb::SourceDBSkymodel &source_db,
    const std::vector<std::string> &patch_names);

/// Create a source list (with patch name) from a patchlist
/// Needed for efficient multithreading
std::vector<std::pair<std::shared_ptr<ModelComponent>, std::shared_ptr<Patch>>>
makeSourceList(std::vector<std::shared_ptr<Patch>> &patchList);

/// From a given PatchList, create a new one with one patch per component
std::vector<std::shared_ptr<Patch>> makeOnePatchPerComponent(
    const std::vector<std::shared_ptr<Patch>> &);

std::vector<std::shared_ptr<Patch>> clusterProximateSources(
    const std::vector<std::shared_ptr<Patch>> &patchList,
    double proximityLimit);

std::vector<std::string> makePatchList(parmdb::SourceDB &sourceDB,
                                       std::vector<std::string> patterns);

std::vector<std::string> MakePatchList(
    const parmdb::SourceDBSkymodel &source_db,
    const std::vector<std::string> &patterns);

/**
 * Creates a list of directions, using packed directions and/or a SourceDB file.
 *
 * @param packed_directions Packed direction lists. Each string should contain a
 * list of directions, e.g., "[direction1, direction2, direction3]".
 * @param source_db_filename A filename of a Source DB file.
 * @return A list of patches.
 */
std::vector<std::vector<std::string>> MakeDirectionList(
    const std::vector<std::string> &packed_directions,
    const std::string &source_db_filename);

bool checkPolarized(parmdb::SourceDB &sourceDB,
                    const std::vector<std::string> &patchNames,
                    unsigned int nModel);

bool CheckPolarized(const parmdb::SourceDBSkymodel &source_db,
                    const std::vector<std::string> &patch_names);

/// Check whether any source in a sourcedb has absolute orientation
bool CheckAnyOrientationIsAbsolute(const parmdb::SourceDBSkymodel &source_db,
                                   const std::vector<std::string> &patch_names);

/// Check whether any source in a skymodel-sourcedb has absolute orientation
bool CheckAnyOrientationIsAbsolute(const parmdb::SourceDBSkymodel &source_db,
                                   const std::vector<std::string> &patch_names);

/// A SourceDB abstraction layer.
///
/// A SourceDB can be read as a SourceDB database or a .skymodel file. This
/// class abstracts the processing of the data regardless of the data source
/// used.
class SourceDB {
 public:
  /// The method used to filter the supplied patches.
  enum class FilterMode {
    /// Filter as a pattern.
    ///
    /// When used the @p patch_names_ is initialised with the found patches
    /// sorted in alphabetic order.
    kPattern,
    /// Filter as a value.
    ///
    /// When used the @p patch_names_ is initialised filter as list of values.
    /// They are stored in the same order as supplied to the constructor.
    kValue
  };

  /// Constructor
  ///
  /// @param source_db_name The name of the source DB to create. The name can
  ///                       either be the name of a binary textual source DB.
  ///                       The entension of the name determines which is used:
  ///                       * If .txt and .skymodel textual
  ///                       * else binary
  /// @param filter         The list of patches to filter. The interpretation
  ///                       of the filter depends on the @a filter_mode.
  /// @param filter_mode    Determines how the @a filter is applied.
  explicit SourceDB(const std::string &source_db_name,
                    const std::vector<std::string> &filter,
                    FilterMode filter_mode);

  std::vector<std::shared_ptr<Patch>> MakePatchList();

  bool CheckPolarized();

  bool CheckAnyOrientationIsAbsolute();

 private:
  std::vector<std::string> patch_names_;
  std::variant<std::monostate, parmdb::SourceDB, parmdb::SourceDBSkymodel>
      source_db_;

  template <class T>
  T &Get() {
    return std::get<T>(source_db_);
  }

  template <class T>
  bool HoldsAlternative() const {
    return std::holds_alternative<T>(source_db_);
  }

  /// Helper for the constructor.
  ///
  /// @pre @code HasSkymodelExtension(source_db_name) @endcode
  ///
  /// The arguments are forwarded from the constructor.
  void InitialiseUsingSkymodel(const std::string &source_db_name,
                               const std::vector<std::string> &filter,
                               FilterMode filter_mode);

  /// Helper for the constructor.
  ///
  /// @pre @code !HasSkymodelExtension(source_db_name) @endcode
  ///
  /// The arguments are forwarded from the constructor.
  void InitialiseUsingSourceDb(const std::string &source_db_name,
                               const std::vector<std::string> &filter,
                               FilterMode filter_mode);
};

}  // namespace base
}  // namespace dp3

#endif
