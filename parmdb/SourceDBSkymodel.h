// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_PARMDB_SOURCEDBSKYMODEL_H
#define DP3_PARMDB_SOURCEDBSKYMODEL_H

#include "SourceDB.h"

namespace dp3 {
namespace parmdb {

/// @ingroup ParmDB
/// @{

class SourceDBSkymodel final : public SourceDBBase {
 public:
  void addSource(const SourceInfo& source_info, const std::string& patch_name,
                 int cat_type, double apparent_brightness,
                 const ParmMap& default_parameters, double ra, double dec,
                 bool check) override;

  void addSource(const SourceInfo& source_info, const std::string& patch_name,
                 const ParmMap& default_parameters, double ra, double dec,
                 bool check) override;

  /// Add a patch and return its patch_id.
  ///
  /// Optionally it is checked if the patch already exists.
  unsigned addPatch(const std::string& patch_name, int cat_type,
                    double apparent_brightness, double ra, double dec,
                    bool check) override;

  void updatePatch(unsigned patch_id, double apparent_brightness, double ra,
                   double dec) override;

  std::vector<std::string> FindPatches(const std::string& pattern) const;
  const std::vector<PatchInfo>& GetPatches() const { return patches_; }

  const PatchInfo& GetPatch(const std::string& patch_name) const;
  const std::vector<SourceData>& GetSources() const { return sources_; }

 private:
  std::vector<PatchInfo> patches_;

  std::map<std::string, unsigned> patches_lut_;

  void ValidatePatchName(const string& patch_name) const;
  unsigned GetPatchRowId(const string& patch_name) const;

  std::vector<SourceData> sources_;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
