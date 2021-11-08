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

#include <boost/variant.hpp>

#include <string>
#include <vector>

namespace dp3 {
namespace base {

std::vector<Patch::ConstPtr> makePatches(
    parmdb::SourceDB &sourceDB, const std::vector<std::string> &patchNames,
    unsigned int nModel);

std::vector<Patch::ConstPtr> MakePatches(
    const parmdb::SourceDBSkymodel &source_db,
    const std::vector<std::string> &patch_names);

/// Create a source list (with patch name) from a patchlist
/// Needed for efficient multithreading
std::vector<std::pair<ModelComponent::ConstPtr, Patch::ConstPtr>>
makeSourceList(const std::vector<Patch::ConstPtr> &patchList);

/// From a given PatchList, create a new one with one patch per component
std::vector<Patch::ConstPtr> makeOnePatchPerComponent(
    const std::vector<Patch::ConstPtr> &);

std::vector<Patch::ConstPtr> clusterProximateSources(
    const std::vector<Patch::ConstPtr> &patchList, double proximityLimit);

std::vector<std::string> makePatchList(parmdb::SourceDB &sourceDB,
                                       std::vector<std::string> patterns);

std::vector<std::string> MakePatchList(
    const parmdb::SourceDBSkymodel &source_db,
    const std::vector<std::string> &patterns);

bool checkPolarized(parmdb::SourceDB &sourceDB,
                    const std::vector<std::string> &patchNames,
                    unsigned int nModel);

bool CheckPolarized(const parmdb::SourceDBSkymodel &source_db,
                    const std::vector<std::string> &patch_names);

/// A SourceDB abstraction layer.
///
/// A SourceDB can be read as a SourceDB database or a .skymodel file. This
/// class abstracts the processing of the data regardless of the data source
/// used.
class SourceDB {
 public:
  explicit SourceDB(const std::string &source_db_name,
                    const std::vector<std::string> &source_patterns);

  std::vector<Patch::ConstPtr> MakePatchList();

  bool CheckPolarized();

 private:
  /// The type of an "empty" variant.
  ///
  /// The original version of this code used boost::variant2, which is modelled
  /// after C++17's variant. However that's not available at Ubuntu 18.04. The
  /// naming used in the code still models the C++17 naming to make porting to
  /// std::variant easier.
  using monostate = int;
  std::vector<std::string> patch_names_;
  boost::variant<monostate, parmdb::SourceDB, parmdb::SourceDBSkymodel>
      source_db_;

  template <class T>
  T &Get() {
    return boost::get<T>(source_db_);
  }

  template <class T>
  bool HoldsAlternative() const {
    return boost::get<T>(&source_db_);
  }
};

}  // namespace base
}  // namespace dp3

#endif
