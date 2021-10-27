// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SourceDBSkymodel.h"

#include <casacore/casa/Utilities/Regex.h>

namespace dp3 {
namespace parmdb {

static void ValidateUniqueName(const std::string& source_name,
                               const std::vector<SourceData>& sourcees) {
  if (std::any_of(sourcees.begin(), sourcees.end(),
                  [&](const SourceData& source) {
                    return source.getInfo().getName() == source_name;
                  }))
    throw std::runtime_error("Source " + source_name + " already exists");
}

void SourceDBSkymodel::addSource(const SourceInfo& source_info,
                                 const string& patch_name, int cat_type,
                                 double apparent_brightness,
                                 const ParmMap& default_parameters, double ra,
                                 double dec, bool check) {
  if (check) {
    ValidateUniqueName(source_info.getName(), sources_);
    ValidatePatchName(patch_name);
  }

  addPatch(patch_name, cat_type, apparent_brightness, ra, dec, false);
  addSource(source_info, patch_name, default_parameters, ra, dec, false);
}

void SourceDBSkymodel::ValidatePatchName(const string& patch_name) const {
  if (patches_lut_.find(patch_name) != patches_lut_.end())
    throw std::runtime_error("Patch " + patch_name + " already exists");
}

unsigned SourceDBSkymodel::GetPatchRowId(const string& patch_name) const {
  auto iter = patches_lut_.find(patch_name);
  if (iter == patches_lut_.end())
    throw std::runtime_error("Patch " + patch_name + " does not exist");
  return iter->second;
}

static SourceData MakeSource(const SourceInfo& source_info,
                             const std::string& patch_name,
                             const ParmMap& default_parameters, double ra,
                             double dec) {
  SourceData result{source_info, patch_name, ra, dec};
  result.setParms(default_parameters);
  return result;
}

void SourceDBSkymodel::addSource(const SourceInfo& source_info,
                                 const string& patch_name,
                                 const ParmMap& default_parameters, double ra,
                                 double dec, bool check) {
  if (check) {
    ValidateUniqueName(source_info.getName(), sources_);
  }
  sources_.push_back(
      MakeSource(source_info, patch_name, default_parameters, ra, dec));
}

const PatchInfo& SourceDBSkymodel::GetPatch(
    const std::string& patch_name) const {
  return patches_[GetPatchRowId(patch_name)];
}

unsigned SourceDBSkymodel::addPatch(const std::string& patch_name, int cat_type,
                                    double apparent_brightness, double ra,
                                    double dec, bool check) {
  if (check) {
    ValidatePatchName(patch_name);
  }
  const unsigned patch_id = patches_.size();
  patches_lut_[patch_name] = patch_id;
  patches_.emplace_back(patch_name, ra, dec, cat_type, apparent_brightness);
  return patch_id;
}

void SourceDBSkymodel::updatePatch(unsigned patch_id,
                                   double apparent_brightness, double ra,
                                   double dec) {
  assert(patch_id < patches_.size() && "Patch index out of bounds.");
  PatchInfo& info = patches_[patch_id];

  info.setApparentBrightness(apparent_brightness);
  info.setRa(ra);
  info.setDec(dec);
}

std::vector<std::string> KeyToVector(
    const std::map<std::string, unsigned>& lut) {
  std::vector<std::string> result;
  std::transform(
      lut.begin(), lut.end(), std::back_inserter(result),
      [](const std::pair<std::string, unsigned>& pair) { return pair.first; });
  return result;
}

std::vector<std::string> SourceDBSkymodel::FindPatches(
    const std::string& pattern) const {
  assert(!pattern.empty() && "The pattern should contain data.");
  if (pattern == "*") {
    return KeyToVector(patches_lut_);
  }

  const casacore::Regex regex{casacore::Regex::fromPattern(pattern)};
  std::vector<std::string> result;
  for (const auto& lut : patches_lut_) {
    const auto& key = lut.first;
    if (regex.match(key.data(), key.size()) == key.size()) {
      result.push_back(key);
    }
  }

  return result;
}

}  // namespace parmdb
}  // namespace dp3
