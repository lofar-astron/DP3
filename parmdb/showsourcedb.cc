// showsourcedb.cc: Show contents of a SourceDB catalog
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// This program shows the contents of a SourceDB catalogs.
// It can show only patches or sources or both.
//
// The program can be run as:
//    showsourcedb  in=inname mode=all|patch|source
// in            name of the SourceDB.
// mode          all      = show patches and sources
//               patch    = only show patches
//               source   = only show sources
//               skymodel = show the output in a skymodel format

#include "SourceDB.h"

#include "../common/StringTools.h"
#include "../common/StreamUtil.h"
#include "../base/Version.h"

#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Inputs/Input.h>

#include <boost/algorithm/string/case_conv.hpp>

#include <vector>

#include <aocommon/logger.h>

using namespace casacore;
using aocommon::Logger;

void show(const string& name, const string& mode, const string& patt) {
  // Open the input SourceDB.
  dp3::parmdb::SourceDB in(dp3::parmdb::ParmDBMeta(string(), name), false,
                           false);
  // Read all patches from the SourceDB and write them.
  vector<dp3::parmdb::PatchInfo> patch(in.getPatchInfo(-1, patt));
  for (size_t i = 0; i < patch.size(); ++i) {
    if (mode != "source") {
      Logger::Info << patch[i] << '\n';
    }
    if (mode != "patch") {
      vector<dp3::parmdb::SourceData> sources(
          in.getPatchSourceData(patch[i].getName()));
      std::ostringstream stream;
      for (dp3::parmdb::SourceData& s : sources) s.print(stream);
      Logger::Info << stream.str();
    }
  }
}

static void showSkymodel(const std::string& name, const std::string& patches) {
  dp3::parmdb::SourceDB source_db(dp3::parmdb::ParmDBMeta("", name), false,
                                  false);
  Logger::Info << "# (Name, Type, Patch, Ra, Dec, I, ReferenceFrequency, "
                  "SpectralIndex='[]', LogarithmicSI, MajorAxis, MinorAxis, "
                  "Orientation, OrientationIsAbsolute) = format\n";

  std::ostringstream stream;
  for (const auto& patch : source_db.getPatchInfo(-1, patches)) {
    dp3::parmdb::toSkymodel(stream, patch);

    for (const dp3::parmdb::SourceData& source :
         source_db.getPatchSourceData(patch.getName())) {
      dp3::parmdb::toSkymodel(stream, source);
    }
  }
  Logger::Info << stream.str();
}

int main(int argc, char* argv[]) {
  try {
    // Define the input parameters.
    Input inputs(1);
    inputs.version(DP3Version::AsString());
    inputs.create("in", "", "Input SourceDB", "string");
    inputs.create("mode", "all",
                  "patch=show all patches, "
                  "source=show all sources, "
                  "all=show patches and sources, "
                  "skymodel=show data as skymodel file",
                  "string");
    inputs.create("patches", "*", "Pattern for names of patches to show",
                  "string");
    // Read and check the input parameters.
    inputs.readArguments(argc, argv);
    string in = inputs.getString("in");
    if (in.empty()) throw std::runtime_error("no input sourcedb name given");
    string mode = boost::to_lower_copy(inputs.getString("mode"));
    if (mode != "patch" && mode != "source" && mode != "all" &&
        mode != "skymodel")
      throw std::runtime_error("incorrect mode given");
    string patt = inputs.getString("patches");
    if (mode == "skymodel")
      showSkymodel(in, patt);
    else
      show(in, mode, patt);
  } catch (AipsError& x) {
    Logger::Error << "Caught AIPS error: " << x.what() << '\n';
    return 1;
  } catch (std::exception& x) {
    Logger::Error << "Caught LOFAR exception: " << x.what() << '\n';
    return 1;
  }
  return 0;
}
