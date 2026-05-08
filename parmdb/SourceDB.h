// SourceDB.h: Base class for a table holding sources and their parameters
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Base class for a table holding sources and their parameters
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_SOURCEDB_H
#define LOFAR_PARMDB_SOURCEDB_H

#include "SourceData.h"
#include "PatchInfo.h"
#include "ParmDBMeta.h"
#include "ParmDB.h"

#include <filesystem>

namespace dp3 {
namespace parmdb {

class ParmMap;

class SourceDBBase {
 public:
  virtual ~SourceDBBase() = default;

  /// Add a source to a patch.
  /// Its ra and dec and default parameters will be stored as default
  /// values in the associated ParmDB tables. The names of the parameters
  /// will be succeeded by a colon and the source name.
  /// The map should contain the parameters belonging to the source type.
  /// Not all parameters need to be present. The ParmDB classes will
  /// use a default of 0 for missing ones.
  ///@{
  virtual void addSource(const SourceInfo& source_info,
                         const std::string& patch_name, int cat_type,
                         double apparent_brightness,
                         const ParmMap& default_parameters, double ra,
                         double dec, bool check) = 0;

  virtual void addSource(const SourceInfo& source_info,
                         const std::string& patch_name,
                         const ParmMap& default_parameters, double ra,
                         double dec, bool check) = 0;
  ///@}

  /// Add a patch and return its patch_id.
  ///
  /// Optionally it is checked if the patch already exists.
  virtual unsigned addPatch(const std::string& patch_name, int cat_type,
                            double apparent_brightness, double ra, double dec,
                            bool check) = 0;

  virtual void updatePatch(unsigned patch_id, double apparent_brightness,
                           double ra, double dec) = 0;
};

}  // namespace parmdb
}  // namespace dp3

#endif
