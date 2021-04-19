// ParameterSet.cc: Implements a map of Key-Value pairs.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Always #include <lofar_config.h> first!

#include "ParameterSet.h"
#include "ParameterRecord.h"

namespace dp3 {
namespace common {

//-------------------------- creation and destroy ---------------------------
//
// Default constructor
//
ParameterSet::ParameterSet(KeyCompare::Mode mode)
    : itsSet(new ParameterSetImpl(mode)) {}

ParameterSet::ParameterSet(bool caseInsensitive)
    : itsSet(new ParameterSetImpl(caseInsensitive ? KeyCompare::NOCASE
                                                  : KeyCompare::NORMAL)) {}

//
// Construction by reading a parameter file.
//
ParameterSet::ParameterSet(const std::string& theFilename, bool caseInsensitive)
    : itsSet(new ParameterSetImpl(theFilename, caseInsensitive
                                                   ? KeyCompare::NOCASE
                                                   : KeyCompare::NORMAL)) {}

//
// Construction by reading a parameter file.
//
ParameterSet::ParameterSet(const std::string& theFilename,
                           KeyCompare::Mode mode)
    : itsSet(new ParameterSetImpl(theFilename, mode)) {}

ParameterSet::ParameterSet(const char* theFilename, KeyCompare::Mode mode)
    : itsSet(new ParameterSetImpl(std::string(theFilename), mode)) {}

//
// Copying is allowed.
//
ParameterSet::ParameterSet(const ParameterSet& that) : itsSet(that.itsSet) {}

//
// operator= copying
//
ParameterSet& ParameterSet::operator=(const ParameterSet& that) {
  if (this != &that) {
    itsSet = that.itsSet;
  }
  return (*this);
}

//
//  Destructor
//
ParameterSet::~ParameterSet() {}

ParameterRecord ParameterSet::getRecord(const std::string& aKey) const {
  return get(aKey).getRecord();
}

//
// operator<<
//
std::ostream& operator<<(std::ostream& os, const ParameterSet& thePS) {
  os << *thePS.itsSet;
  return os;
}

}  // namespace common
}  // namespace dp3
