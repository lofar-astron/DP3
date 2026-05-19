// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// \file
/// Helper functions to extract patch and source information from a SourceDB.

#ifndef DP3_SKY_MODEL_SKY_MODEL_FUNCTIONS_H_
#define DP3_SKY_MODEL_SKY_MODEL_FUNCTIONS_H_

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "sky_model/SkyModel.h"

#include "Patch.h"

namespace dp3::sky_model {

std::vector<std::shared_ptr<Patch>> MakePatches(
    const sky_model::SkyModel &source_db,
    const std::vector<std::string> &patch_names);

/// Create a source list (with patch name) from a patchlist
/// Needed for efficient multithreading
std::vector<
    std::pair<std::shared_ptr<base::ModelComponent>, std::shared_ptr<Patch>>>
makeSourceList(std::vector<std::shared_ptr<Patch>> &patchList);

/**
 * Makes sure that every Patch's stored index is set to the index the patch has
 * in the list. @see Patch::Index().
 */
void SetPatchIndices(std::vector<std::shared_ptr<Patch>> &patch_list);

/// From a given patch list, create a new one with one patch per component
std::vector<std::shared_ptr<Patch>> makeOnePatchPerComponent(
    const std::vector<std::shared_ptr<Patch>> &);

std::vector<std::shared_ptr<Patch>> clusterProximateSources(
    const std::vector<std::shared_ptr<Patch>> &patchList,
    double proximityLimit);

std::vector<std::string> MakePatchList(
    const sky_model::SkyModel &source_db,
    const std::vector<std::string> &patterns);

/**
 * Creates a list of directions, using packed directions and/or a SourceDB file.
 *
 * @param packed_directions Packed direction lists. Each string should contain a
 * list of directions, e.g., "[direction1, direction2, direction3]".
 * @param sky_model_filename A filename of a sky model file.
 * @return A list of patches.
 */
std::vector<std::vector<std::string>> MakeDirectionList(
    const std::vector<std::string> &packed_directions,
    const std::string &sky_model_filename);

bool CheckPolarized(const sky_model::SkyModel &source_db,
                    const std::vector<std::string> &patch_names);

/// Check whether any source in a skymodel-sourcedb has absolute orientation
bool CheckAnyOrientationIsAbsolute(const sky_model::SkyModel &sky_model,
                                   const std::vector<std::string> &patch_names);

}  // namespace dp3::sky_model

#endif
