// ParmDBBlob.cc: Dummy class to hold parmaeter values
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParmDBBlob.h"

using namespace std;

namespace dp3 {
namespace parmdb {

ParmDBBlob::ParmDBBlob(const string&, bool) {}

ParmDBBlob::~ParmDBBlob() {}

void ParmDBBlob::flush(bool) {}

void ParmDBBlob::lock(bool) {}

void ParmDBBlob::unlock() {}

void ParmDBBlob::clearTables() {}

void ParmDBBlob::setDefaultSteps(const vector<double>&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

int ParmDBBlob::getNameId(const std::string&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

Box ParmDBBlob::getRange(const string&) const {
  throw std::runtime_error("ParmDBBlob not implemented");
}

Box ParmDBBlob::getRange(const std::vector<std::string>&) const {
  throw std::runtime_error("ParmDBBlob not implemented");
}

void ParmDBBlob::getValues(vector<ParmValueSet>&, const vector<unsigned int>&,
                           const vector<ParmId>&, const Box&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

void ParmDBBlob::getDefValues(ParmMap&, const std::string&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

void ParmDBBlob::putValues(const std::string&, int&, ParmValueSet&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

void ParmDBBlob::putDefValue(const std::string&, const ParmValueSet&, bool) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

void ParmDBBlob::deleteValues(const std::string&, const Box&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

void ParmDBBlob::deleteDefValues(const std::string&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

vector<string> ParmDBBlob::getNames(const std::string&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

void ParmDBBlob::fillDefMap(ParmMap&) {
  throw std::runtime_error("ParmDBBlob not implemented");
}

}  // namespace parmdb
}  // namespace dp3
