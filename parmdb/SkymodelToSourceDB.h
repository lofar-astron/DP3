// Skymodel.h:
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief  Class holding all functions needed to convert a .skymodel text file
/// into a .sourcedb directory
/// @author Chiara Salvoni

#ifndef DP3_SKYMODEL_TO_SOURCEDB_H
#define DP3_SKYMODEL_TO_SOURCEDB_H

#include "SourceDB.h"
#include "SourceDBSkymodel.h"

#include "common/StringTools.h"
#include "common/StreamUtil.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace dp3::parmdb::skymodel_to_source_db {

SourceDBSkymodel MakeSourceDBSkymodel(const std::string& filename,
                                      const std::string& format);

std::string ReadFormat(const std::string& file, const std::string& cat_file);

}  // namespace dp3::parmdb::skymodel_to_source_db

#endif
