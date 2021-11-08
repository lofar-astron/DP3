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

#include "../common/StringTools.h"
#include "../common/StreamUtil.h"

#include <boost/algorithm/string.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Inputs/Input.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Utilities/Regex.h>

struct SearchInfo {
  double ra;
  double dec;
  double sinDec;
  double cosDec;
  double cosRadius;
  double raStart;
  double raEnd;
  double decStart;
  double decEnd;
  bool search;  // false no search info given
  bool asCone;  // true is search in a cone, otherwise a box
};

namespace dp3 {

namespace parmdb {

namespace skymodel_to_source_db {

SourceDB MakeSourceDb(const std::string& in, const std::string& out,
                      const std::string& outType, const std::string& format,
                      const std::string& prefix, const std::string& suffix,
                      bool append, bool average, bool check,
                      const SearchInfo& search_info);

SourceDBSkymodel MakeSourceDBSkymodel(const std::string& filename,
                                      const std::string& format);

SearchInfo GetSearchInfo(const std::string& center, const std::string& radius,
                         const std::string& width);

std::string ReadFormat(std::string file, const std::string& cat_file);

}  // namespace skymodel_to_source_db
}  // namespace parmdb
}  // namespace dp3

#endif
