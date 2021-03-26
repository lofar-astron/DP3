// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParameterHandler.h"

using namespace std;

namespace dp3 {
namespace common {

ParameterHandler::ParameterHandler(const ParameterSet& parSet)
    : itsParms(parSet) {}

string ParameterHandler::getString(const string& parm,
                                   const string& defVal) const {
  if (itsParms.isDefined(parm)) {
    return itsParms.getString(parm);
  }
  return defVal;
}

double ParameterHandler::getDouble(const string& parm, double defVal) const {
  if (itsParms.isDefined(parm)) {
    return itsParms.getDouble(parm);
  }
  return defVal;
}

unsigned ParameterHandler::getUint(const string& parm, unsigned defVal) const {
  if (itsParms.isDefined(parm)) {
    return itsParms.getUint32(parm);
  }
  return defVal;
}

bool ParameterHandler::getBool(const string& parm, bool defVal) const {
  if (itsParms.isDefined(parm)) {
    return itsParms.getBool(parm);
  }
  return defVal;
}

vector<string> ParameterHandler::getStringVector(
    const string& parm, const vector<string>& defVal) const {
  if (itsParms.isDefined(parm)) {
    return itsParms.getStringVector(parm);
  }
  return defVal;
}

void ParameterHandler::fillString(const string& parm, string& value) const {
  if (itsParms.isDefined(parm)) {
    value = itsParms.getString(parm);
  }
}

void ParameterHandler::fillDouble(const string& parm, double& value) const {
  if (itsParms.isDefined(parm)) {
    value = itsParms.getDouble(parm);
  }
}

void ParameterHandler::fillUint(const string& parm, unsigned& value) const {
  if (itsParms.isDefined(parm)) {
    value = itsParms.getUint32(parm);
  }
}

void ParameterHandler::fillBool(const string& parm, bool& value) const {
  if (itsParms.isDefined(parm)) {
    value = itsParms.getBool(parm);
  }
}

void ParameterHandler::fillStringVector(const string& parm,
                                        vector<string>& value) const {
  if (itsParms.isDefined(parm)) {
    value = itsParms.getStringVector(parm);
  }
}

blob::BlobOStream& operator<<(blob::BlobOStream& bs, const ParameterSet& m) {
  bs.putStart("ParameterSet", 1);
  bs << static_cast<uint32_t>(m.size());
  for (ParameterSet::const_iterator it = m.begin(); it != m.end(); ++it) {
    bs << it->first << it->second.get();
  }
  bs.putEnd();
  return bs;
}

blob::BlobIStream& operator>>(blob::BlobIStream& bs, ParameterSet& m) {
  bs.getStart("ParameterSet");
  m.clear();
  uint32_t size;
  bs >> size;
  std::string k, v;
  for (uint32_t i = 0; i < size; ++i) {
    bs >> k >> v;
    m.add(k, v);
  }
  bs.getEnd();
  return bs;
}

}  // namespace common
}  // namespace dp3
