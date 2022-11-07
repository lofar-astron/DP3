// InputStep.cc: Abstract base class for a Step generating input
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "InputStep.h"
#include "MultiMSReader.h"
#include "MSReader.h"
#include "MSBDAReader.h"

#include "../base/MS.h"
#include "../common/ParameterSet.h"

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/Utilities/Copy.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using casacore::Cube;
using casacore::IPosition;
using casacore::Matrix;
using casacore::RefRows;

namespace dp3 {
namespace steps {

InputStep::~InputStep() {}

std::string InputStep::msName() const { return std::string(); }

const casacore::Table& InputStep::table() const {
  throw std::runtime_error("InputStep::table not implemented");
}

double InputStep::firstTime() const {
  throw std::runtime_error("InputStep::firstTime not implemented");
}

double InputStep::lastTime() const {
  throw std::runtime_error("InputStep::lastTime not implemented");
}

bool InputStep::HasBda(const casacore::MeasurementSet& ms) {
  return ms.keywordSet().isDefined(base::DP3MS::kBDAFactorsTable) &&
         (ms.keywordSet().asTable(base::DP3MS::kBDAFactorsTable).nrow() > 0);
}

std::unique_ptr<InputStep> InputStep::CreateReader(
    const common::ParameterSet& parset) {
  // Get input and output MS name.
  // Those parameters were always called msin and msout.
  // However, SAS/MAC cannot handle a parameter and a group with the same
  // name, hence one can also use msin.name and msout.name.
  std::vector<std::string> inNames =
      parset.getStringVector("msin.name", std::vector<std::string>());
  if (inNames.empty()) {
    inNames = parset.getStringVector("msin");
  }
  if (inNames.empty())
    throw std::runtime_error("No input MeasurementSets given");
  // Find all file names matching a possibly wildcarded input name.
  // This is only possible if a single name is given.
  if (inNames.size() == 1) {
    if (inNames[0].find_first_of("*?{['") != std::string::npos) {
      std::vector<std::string> names;
      names.reserve(80);
      casacore::Path path(inNames[0]);
      casacore::String dirName(path.dirName());
      casacore::Directory dir(dirName);
      // Use the basename as the file name pattern.
      casacore::DirectoryIterator dirIter(
          dir, casacore::Regex(casacore::Regex::fromPattern(path.baseName())));
      while (!dirIter.pastEnd()) {
        names.push_back(dirName + '/' + dirIter.name());
        dirIter++;
      }
      if (names.empty())
        throw std::runtime_error("No datasets found matching msin " +
                                 inNames[0]);
      inNames = names;
    }
  }

  if (inNames.size() == 1) {
    if (!casacore::Table::isReadable(inNames.front())) {
      throw std::invalid_argument("No such MS: " + inNames.front());
    }
    const casacore::MeasurementSet ms(inNames.front(),
                                      casacore::TableLock::AutoNoReadLocking);
    if (HasBda(ms)) {
      return std::make_unique<MSBDAReader>(ms, parset, "msin.");
    } else {
      return std::make_unique<MSReader>(ms, parset, "msin.");
    }
  } else {
    // MultiMSReader checks that all MS's have regular (non-BDA) data.
    return std::make_unique<MultiMSReader>(inNames, parset, "msin.");
  }
}

}  // namespace steps
}  // namespace dp3
