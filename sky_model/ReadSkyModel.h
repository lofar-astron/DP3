// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief  Class holding all functions needed to convert a .sky_model text file
/// into a .sourcedb directory
/// @author Chiara Salvoni

#ifndef DP3_SKY_MODEL_READ_SKY_MODEL_H_
#define DP3_SKY_MODEL_READ_SKY_MODEL_H_

#include "SkyModel.h"

#include "common/StringTools.h"
#include "common/StreamUtil.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace dp3::sky_model {

std::string ReadFormat(const std::string& file, const std::string& cat_file);

SkyModel ReadSkyModel(const std::string& filename, const std::string& format);

inline SkyModel ReadSkyModel(const std::string& filename) {
  return ReadSkyModel(filename, sky_model::ReadFormat("", filename));
}

}  // namespace dp3::sky_model

#endif
