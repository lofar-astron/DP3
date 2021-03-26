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

#include <vector>

namespace dp3 {
namespace parmdb {
class SourceDB;
}

namespace base {

std::vector<Patch::ConstPtr> makePatches(
    parmdb::SourceDB &sourceDB, const std::vector<std::string> &patchNames,
    unsigned int nModel);

/// Create a source list (with patch name) from a patchlist
/// Needed for efficient multithreading
std::vector<std::pair<ModelComponent::ConstPtr, Patch::ConstPtr> >
makeSourceList(const std::vector<Patch::ConstPtr> &patchList);

/// From a given PatchList, create a new one with one patch per component
std::vector<Patch::ConstPtr> makeOnePatchPerComponent(
    const std::vector<Patch::ConstPtr> &);

std::vector<Patch::ConstPtr> clusterProximateSources(
    const std::vector<Patch::ConstPtr> &patchList, double proximityLimit);

std::vector<std::string> makePatchList(parmdb::SourceDB &sourceDB,
                                       std::vector<std::string> patterns);

bool checkPolarized(parmdb::SourceDB &sourceDB,
                    const std::vector<std::string> &patchNames,
                    unsigned int nModel);

}  // namespace base
}  // namespace dp3

#endif
