// InputStep.cc: Abstract base class for a Step generating input
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../base/Exceptions.h"

#include "InputStep.h"
#include "MultiMSReader.h"
#include "MSReader.h"
#include "MSBDAReader.h"

#include "../common/ParameterSet.h"

#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/casa/Utilities/Copy.h>

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

void InputStep::getModelData(const RefRows&, Cube<casacore::Complex>&) {
  throw Exception("InputStep::getModelData not implemented");
}

void InputStep::fillBeamInfo(std::vector<std::shared_ptr<everybeam::Station>>&,
                             const casacore::Vector<casacore::String>&,
                             const everybeam::ElementResponseModel) const {
  throw Exception("InputStep::fillBeamInfo not implemented");
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

std::unique_ptr<InputStep> InputStep::CreateReader(
    const common::ParameterSet& parset, const std::string& prefix) {
  // Get input and output MS name.
  // Those parameters were always called msin and msout.
  // However, SAS/MAC cannot handle a parameter and a group with the same
  // name, hence one can also use msin.name and msout.name.
  std::vector<std::string> inNames =
      parset.getStringVector("msin.name", std::vector<std::string>());
  if (inNames.empty()) {
    inNames = parset.getStringVector("msin");
  }
  if (inNames.size() == 0) throw Exception("No input MeasurementSets given");
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

  if (!parset.getBool(prefix + "bda", false)) {
    // Get the steps.
    // Currently the input MS must be given.
    // In the future it might be possible to have a simulation step instead.
    // Create MSReader step if input ms given.
    if (inNames.size() == 1) {
      return boost::make_unique<MSReader>(inNames[0], parset, "msin.");
    } else {
      return boost::make_unique<MultiMSReader>(inNames, parset, "msin.");
    }
  } else {
    if (inNames.size() > 1) {
      throw std::invalid_argument(
          "DP3 does not support multiple in MS for BDA data.");
    }
    return boost::make_unique<MSBDAReader>(inNames[0], parset, "msin.");
  }
}

}  // namespace steps
}  // namespace dp3
