// ParmDB.cc: Object to hold parameters in a table.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParmDB.h"
#include "ParmDBCasa.h"
#include "ParmDBBlob.h"

#include <casacore/casa/Utilities/Regex.h>

#ifdef AIPS_NO_TEMPLATE_SRC
#include <casacore/casa/Utilities/GenSort.cc>  // for automatic template
#else
#include <casacore/casa/Utilities/GenSort.h>
#endif

using namespace std;
using namespace casacore;

namespace dp3 {
namespace parmdb {

map<string, int> ParmDB::theirDBNames;
vector<ParmDBRep*> ParmDB::theirParmDBs;

ParmDBRep::ParmDBRep()
    : itsCount(0), itsSeqNr(-1), itsDefFilled(false), itsDefSteps(2) {
  itsDefSteps[0] = 1000;  // 1 KHz
  itsDefSteps[1] = 5;     // 5 sec
}

ParmDBRep::~ParmDBRep() {}

void ParmDBRep::flush(bool) {}

void ParmDBRep::lock(bool) {}

void ParmDBRep::unlock() {}

ParmValueSet ParmDBRep::getDefValue(const std::string& parmName,
                                    const ParmValue& defaultValue) {
  // Fill the map with default values if not done yet.
  if (!itsDefFilled) {
    fillDefMap(itsDefValues);
    itsDefFilled = true;
  }
  // Try to find the default value.
  // The parameter name consists of parts (separated by colons), so the
  // parameters are categorised in that way.
  // An initial value can be defined for the full name or for a higher
  // category.
  // So look up until found or until no more parts are left.
  string name = parmName;
  while (true) {
    ParmMap::const_iterator pos = itsDefValues.find(name);
    if (pos != itsDefValues.end()) {
      return pos->second;
    }
    std::string::size_type idx = name.rfind(':');
    // Exit loop if no more name parts.
    if (idx == std::string::npos) {
      break;
    }
    // Remove last part and try again.
    name = name.substr(0, idx);
  }
  // Nothing found; return the default ParmValue.
  return ParmValueSet(defaultValue);
}

void ParmDBRep::getValuesPattern(ParmMap& result,
                                 const std::string& parmNamePattern,
                                 const Box& domain) {
  vector<std::string> parmNames = getNames(parmNamePattern);
  vector<unsigned int> nameIds;
  nameIds.reserve(parmNames.size());
  for (unsigned int i = 0; i < parmNames.size(); ++i) {
    int id = getNameId(parmNames[i]);
    if (id >= 0) {
      nameIds.push_back(id);
    }
  }
  vector<ParmId> parmIds;
  parmIds.reserve(nameIds.size());
  for (unsigned int i = 0; i < nameIds.size(); ++i) {
    parmIds.push_back(i);
  }
  vector<ParmValueSet> sets(nameIds.size());
  getValues(sets, nameIds, parmIds, domain);
  for (unsigned int i = 0; i < sets.size(); ++i) {
    if (sets[i].size() > 0) {
      result.define(parmNames[i], sets[i]);
    }
  }
}

ParmDB::ParmDB(const ParmDBMeta& ptm, bool forceNew) {
  // Attach to existing one if already opened.
  map<string, int>::iterator pos = theirDBNames.find(ptm.getTableName());
  if (pos != theirDBNames.end()) {
    itsRep = theirParmDBs[pos->second];
    itsRep->link();
    return;
  }
  // Open the correct ParmDB.
  if (ptm.getType() == "casa") {
    itsRep = new ParmDBCasa(ptm.getTableName(), forceNew);
  } else if (ptm.getType() == "blob") {
    itsRep = new ParmDBBlob(ptm.getTableName(), forceNew);
    ///  } else if (ptm.getType() == "bdb") {
    /// itsRep = new ParmDBBDB (ptm, forceNew);
  } else if (ptm.getType() == "postgres") {
    //#if defined(HAVE_PGSQL)
#if 0
      itsRep = new ParmDBPostgres(ptm.getDBName(),
                                  ptm.getUserName(),
                                  ptm.getDBPwd(),
                                  ptm.getHostName(),
                                  "");
#else
    throw std::runtime_error("unsupported parmTableType: " + ptm.getType());
#endif
  } else {
    throw std::runtime_error("unknown parmTableType: " + ptm.getType());
  }
  itsRep->link();
  itsRep->setParmDBMeta(ptm);
  // Get the sequence number of the ParmDBs opened.
  unsigned int dbnr = theirParmDBs.size();
  if (dbnr == theirDBNames.size()) {
    theirParmDBs.push_back(itsRep);
  } else {
    // Some entry has been deleted; reuse it.
    for (dbnr = 0; dbnr < theirParmDBs.size(); ++dbnr) {
      if (theirParmDBs[dbnr] == nullptr) {
        theirParmDBs[dbnr] = itsRep;
        break;
      }
    }
  }
  itsRep->setParmDBSeqNr(dbnr);
  theirDBNames.insert(make_pair(ptm.getTableName(), dbnr));
}

ParmDB::ParmDB(ParmDBRep* rep) : itsRep(rep) { itsRep->link(); }

ParmDB::ParmDB(const ParmDB& that) : itsRep(that.itsRep) { itsRep->link(); }

ParmDB& ParmDB::operator=(const ParmDB& that) {
  if (this != &that) {
    decrCount();
    itsRep = that.itsRep;
    itsRep->link();
  }
  return *this;
}

void ParmDB::decrCount() {
  if (itsRep->unlink() == 0) {
    string tabName = itsRep->getParmDBMeta().getTableName();
    map<string, int>::iterator pos = theirDBNames.find(tabName);
    if (pos == theirDBNames.end())
      throw std::runtime_error("~ParmDB " + tabName + " not found in map");
    assert(theirParmDBs[pos->second] == itsRep);
    theirParmDBs[pos->second] = nullptr;
    theirDBNames.erase(pos);
    delete itsRep;
    itsRep = nullptr;
  }
}

ParmDB ParmDB::getParmDB(unsigned int index) {
  if (index >= theirParmDBs.size() || theirParmDBs[index] == nullptr)
    throw std::runtime_error("ParmDB index " + std::to_string(index) +
                             " is unknown");
  return ParmDB(theirParmDBs[index]);
}

}  // namespace parmdb
}  // namespace dp3
