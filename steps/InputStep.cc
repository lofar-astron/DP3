// InputStep.cc: Abstract base class for a Step generating input
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "InputStep.h"
#include "MultiMSReader.h"
#include "MSReader.h"
#include "MSBDAReader.h"

#include "../base/Exceptions.h"
#include "../base/MS.h"
#include "../common/ParameterSet.h"

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/Utilities/Copy.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <boost/make_unique.hpp>

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

const Cube<bool>& InputStep::fetchFullResFlags(const base::DPBuffer& bufin,
                                               base::DPBuffer& bufout,
                                               common::NSTimer& timer,
                                               bool merge) {
  // If already defined in the buffer, return those fullRes flags.
  if (!bufin.getFullResFlags().empty()) {
    return bufin.getFullResFlags();
  }
  // No fullRes flags in buffer, so get them from the input.
  timer.stop();
  bool fnd = getFullResFlags(bufin.getRowNrs(), bufout);
  timer.start();
  Cube<bool>& fullResFlags = bufout.getFullResFlags();
  if (!fnd) {
    // No fullRes flags in input; form them from the flags in the buffer.
    // Only use the XX flags; no averaging done, thus navgtime=1.
    // (If any averaging was done, the flags would be in the buffer).
    IPosition shp(bufin.getFlags().shape());
    if (fullResFlags.shape()[0] != shp[1] || fullResFlags.shape()[1] != 1 ||
        fullResFlags.shape()[2] != shp[2])
      throw std::runtime_error("Invalid shape of full res flags");
    casacore::objcopy(fullResFlags.data(), bufin.getFlags().data(),
                      fullResFlags.size(), 1, shp[0]);  // only copy XX.
    return fullResFlags;
  }
  // There are fullRes flags.
  // If needed, merge them with the buffer's flags.
  if (merge) {
    DPBuffer::mergeFullResFlags(fullResFlags, bufin.getFlags());
  }
  return fullResFlags;
}

const Cube<float>& InputStep::fetchWeights(const DPBuffer& bufin,
                                           DPBuffer& bufout,
                                           common::NSTimer& timer) {
  // If already defined in the buffer, return those weights.
  if (!bufin.getWeights().empty()) {
    return bufin.getWeights();
  }
  // No weights in buffer, so get them from the input.
  // It might need the data and flags in the buffer.
  timer.stop();
  getWeights(bufin.getRowNrs(), bufout);
  timer.start();
  return bufout.getWeights();
}

const Matrix<double>& InputStep::fetchUVW(const DPBuffer& bufin,
                                          DPBuffer& bufout,
                                          common::NSTimer& timer) {
  // If already defined in the buffer, return those UVW.
  if (!bufin.getUVW().empty()) {
    return bufin.getUVW();
  }
  // No UVW in buffer, so get them from the input.
  timer.stop();
  getUVW(bufin.getRowNrs(), bufin.getTime(), bufout);
  timer.start();
  return bufout.getUVW();
}

void InputStep::getUVW(const RefRows&, double, DPBuffer&) {
  throw Exception("InputStep::getUVW not implemented");
}

void InputStep::getWeights(const RefRows&, DPBuffer&) {
  throw Exception("InputStep::getWeights not implemented");
}

bool InputStep::getFullResFlags(const RefRows&, DPBuffer&) {
  throw Exception("InputStep::getFullResFlags not implemented");
}

std::unique_ptr<everybeam::telescope::Telescope> InputStep::GetTelescope(
    const everybeam::ElementResponseModel /*element_response_model*/,
    bool /*use_channel_frequency*/) const {
  throw Exception("InputStep::GetTelescope not implemented");
}

std::vector<size_t> InputStep::SelectStationIndices(
    const everybeam::telescope::Telescope* telescope,
    const casacore::Vector<casacore::String>& station_names) {
  const everybeam::telescope::PhasedArray* phased_array =
      dynamic_cast<const everybeam::telescope::PhasedArray*>(telescope);
  if (phased_array == nullptr) {
    throw Exception(
        "Currently, only PhasedArray telescopes accepted as input, i.e. "
        "OSKAR or LOFAR. Support for other telescope may become available "
        "soon.");
  }

  // Copy only those stations for which the name matches.
  // Note: the order of the station names in both vectors match,
  // thus avoiding a nested loop.
  std::vector<size_t> station_to_msindex;
  station_to_msindex.reserve(station_names.size());
  size_t station_idx = 0;
  for (size_t i = 0; i < phased_array->GetNrStations(); ++i) {
    if (station_idx < station_names.size() &&
        casacore::String(phased_array->GetStation(i)->GetName()) ==
            station_names[station_idx]) {
      station_to_msindex.push_back(i);
      station_idx++;
    }
  }

  if (station_idx != station_names.size()) {
    throw Exception(
        "InputStep::SelectStationIndices -"
        " some stations miss the beam info");
  }
  return station_to_msindex;
}

void InputStep::setReadVisData(bool) {
  throw Exception("InputStep::setReadVisData not implemented");
}

const casacore::Table& InputStep::table() const {
  throw Exception("InputStep::table not implemented");
}

double InputStep::firstTime() const {
  throw Exception("InputStep::firstTime not implemented");
}

double InputStep::lastTime() const {
  throw Exception("InputStep::lastTime not implemented");
}

unsigned int InputStep::spectralWindow() const {
  throw Exception("InputStep::spectralWindow not implemented");
}

unsigned int InputStep::nchanAvgFullRes() const {
  throw Exception("InputStep::nchanAvgFullRes not implemented");
}

unsigned int InputStep::ntimeAvgFullRes() const {
  throw Exception("InputStep::ntimeAvgFullRes not implemented");
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
  if (inNames.empty()) throw Exception("No input MeasurementSets given");
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
        throw Exception("No datasets found matching msin " + inNames[0]);
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
      return boost::make_unique<MSBDAReader>(ms, parset, "msin.");
    } else {
      return boost::make_unique<MSReader>(ms, parset, "msin.");
    }
  } else {
    // MultiMSReader checks that all MS's have regular (non-BDA) data.
    return boost::make_unique<MultiMSReader>(inNames, parset, "msin.");
  }
}

}  // namespace steps
}  // namespace dp3
