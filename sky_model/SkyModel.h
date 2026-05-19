// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_PARMDB_SKY_MODEL_H
#define DP3_PARMDB_SKY_MODEL_H

#include <string>

#include "../parmdb/ParmValue.h"
#include "PatchInfo.h"
#include "Source.h"

namespace dp3::sky_model {

class SkyModel {
 public:
  void addSource(
      const SourceInfo& source_info, const std::string& patch_name,
      int cat_type, double apparent_brightness,
      const std::map<std::string, parmdb::ParmValue>& default_parameters,
      double ra, double dec, bool check);

  void addSource(
      const SourceInfo& source_info, const std::string& patch_name,
      const std::map<std::string, parmdb::ParmValue>& default_parameters,
      double ra, double dec, bool check);

  /// Add a patch and return its patch_id.
  ///
  /// Optionally it is checked if the patch already exists.
  unsigned addPatch(const std::string& patch_name, int cat_type,
                    double apparent_brightness, double ra, double dec,
                    bool check);

  void updatePatch(unsigned patch_id, double apparent_brightness, double ra,
                   double dec);

  std::vector<std::string> FindPatches(const std::string& pattern) const;
  std::vector<std::string> GetPatchNames() const;
  const std::vector<PatchInfo>& GetPatches() const { return patches_; }

  const PatchInfo& GetPatch(const std::string& patch_name) const;
  const std::vector<Source>& GetSources() const { return sources_; }

 private:
  std::vector<PatchInfo> patches_;

  std::map<std::string, unsigned> patches_lut_;

  void ValidatePatchName(const std::string& patch_name) const;
  unsigned GetPatchRowId(const std::string& patch_name) const;

  std::vector<Source> sources_;
};

}  // namespace dp3::sky_model

#endif
