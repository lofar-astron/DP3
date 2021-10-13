// SourceDBBlob.cc: Class for a Blob file holding sources and their parameters
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SourceDBBlob.h"
#include "ParmMap.h"

#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Utilities/Regex.h>
#include <casacore/casa/Utilities/Sort.h>

#include <iostream>

using namespace casacore;
using namespace std;

namespace dp3 {
namespace parmdb {

SourceDBBlob::SourceDBBlob(const ParmDBMeta& pdm, bool forceNew)
    : SourceDBRep(pdm, forceNew), itsCanWrite(true), itsEndPos(0) {
  if (!forceNew) {
    // First open as readonly to see if it exists.
    itsFile.open(pdm.getTableName().c_str(), ios::in | ios::binary);
    if (!itsFile) {
      forceNew = true;
    } else {
      // Get eof position.
      itsFile.seekg(0, ios::end);
      itsEndPos = itsFile.tellg();
      itsFile.close();
      // See if it can be opened as read/write.
      // If not, open it again as readonly.
      itsFile.open(pdm.getTableName().c_str(),
                   ios::in | ios::out | ios::binary);
      if (!itsFile) {
        itsFile.open(pdm.getTableName().c_str(), ios::in | ios::binary);
        itsCanWrite = false;
      }
      if (!itsFile)
        throw std::runtime_error("Error opening file " + pdm.getTableName());
    }
  }
  if (forceNew) {
    itsFile.open(pdm.getTableName().c_str(),
                 ios::out | ios::in | ios::trunc | ios::binary);
    if (!itsFile)
      throw std::runtime_error("SourceDB blob file " + pdm.getTableName() +
                               " cannot be created");
  }
  itsBufIn =
      std::shared_ptr<blob::BlobIBufStream>(new blob::BlobIBufStream(itsFile));
  itsBufOut =
      std::shared_ptr<blob::BlobOBufStream>(new blob::BlobOBufStream(itsFile));
  itsBlobIn =
      std::shared_ptr<blob::BlobIStream>(new blob::BlobIStream(*itsBufIn));
  itsBlobOut =
      std::shared_ptr<blob::BlobOStream>(new blob::BlobOStream(*itsBufOut));
}

SourceDBBlob::~SourceDBBlob() {}

void SourceDBBlob::lock(bool) {}

void SourceDBBlob::unlock() {}

void SourceDBBlob::clearTables() {
  // Recreate the file.
}

void SourceDBBlob::checkDuplicates() {}

vector<string> SourceDBBlob::findDuplicatePatches() { return vector<string>(); }

vector<string> SourceDBBlob::findDuplicateSources() { return vector<string>(); }

bool SourceDBBlob::patchExists(const string&) { return false; }

bool SourceDBBlob::sourceExists(const string&) { return false; }

unsigned int SourceDBBlob::addPatch(const string& patchName, int catType,
                                    double apparentBrightness, double ra,
                                    double dec, bool) {
  // Write at the end of the file.
  if (!itsCanWrite)
    throw std::runtime_error("SourceDBBlob: file is not writable");
  itsFile.seekp(0, ios::end);
  unsigned int filePos = itsFile.tellp();
  *itsBlobOut << PatchInfo(patchName, ra, dec, catType, apparentBrightness);
  itsEndPos = itsFile.tellp();
  return filePos;
}

void SourceDBBlob::updatePatch(unsigned int filePos, double apparentBrightness,
                               double ra, double dec) {
  // First read the blob to get name and category.
  itsFile.seekp(filePos, ios::beg);
  PatchInfo info;
  *itsBlobIn >> info;
  // Write the calculated position and flux of the patches.
  info.setRa(ra);
  info.setDec(dec);
  info.setApparentBrightness(apparentBrightness);
  itsFile.seekp(filePos, ios::beg);
  *itsBlobOut << info;
}

void SourceDBBlob::addSource(const SourceInfo& sourceInfo,
                             const string& patchName,
                             const ParmMap& defaultParameters, double ra,
                             double dec, bool) {
  if (!itsCanWrite)
    throw std::runtime_error("SourceDBBlob: file is not writable");
  itsFile.seekp(0, ios::end);
  SourceData src(sourceInfo, patchName, ra, dec);
  src.setParms(defaultParameters);
  src.writeSource(*itsBlobOut);
  itsEndPos = itsFile.tellp();
}

void SourceDBBlob::addSource(const SourceData& source, bool) {
  source.writeSource(*itsBlobOut);
}

void SourceDBBlob::addSource(const SourceInfo& sourceInfo,
                             const string& patchName, int catType,
                             double apparentBrightness,
                             const ParmMap& defaultParameters, double ra,
                             double dec, bool check) {
  addPatch(patchName, catType, apparentBrightness, ra, dec, check);
  addSource(sourceInfo, patchName, defaultParameters, ra, dec, check);
}

void SourceDBBlob::deleteSources(const string&) {
  throw std::runtime_error("SourceDBBlob::deleteSources not possible");
}

vector<string> SourceDBBlob::getPatches(int category, const string& pattern,
                                        double minBrightness,
                                        double maxBrightness) {
  // If not done yet, read all data from the file.
  readAll();
  Regex regex;
  if (pattern.size() > 0) {
    regex = Regex::fromPattern(pattern);
  }
  // Fill the patch names selecting only the required ones.
  std::vector<std::string> names;
  std::vector<Int> categories;
  std::vector<double> brightness;
  names.reserve(itsPatches.size());
  categories.reserve(itsPatches.size());
  brightness.reserve(itsPatches.size());
  for (map<std::string, PatchInfo>::const_iterator iter = itsPatches.begin();
       iter != itsPatches.end(); ++iter) {
    if ((category < 0 || iter->second.getCategory() == category) &&
        (minBrightness < 0 ||
         iter->second.apparentBrightness() >= minBrightness) &&
        (maxBrightness < 0 ||
         iter->second.apparentBrightness() <= maxBrightness) &&
        (pattern.size() == 0 || String(iter->first).matches(regex))) {
      names.push_back(iter->first);
      categories.push_back(iter->second.getCategory());
      brightness.push_back(iter->second.apparentBrightness());
    }
  }
  // Sort in order of category, brightness, name.
  vector<string> nmout;
  if (!names.empty()) {
    Sort sort;
    sort.sortKey(&(categories[0]), TpInt);
    sort.sortKey(&(brightness[0]), TpDouble, 0, Sort::Descending);
    sort.sortKey(&(names[0]), TpString);
    Vector<uInt> index(names.size());
    sort.sort(index, names.size());
    nmout.reserve(names.size());
    for (unsigned int i = 0; i < names.size(); ++i) {
      nmout.push_back(names[index[i]]);
    }
  }
  return nmout;
}

vector<PatchInfo> SourceDBBlob::getPatchInfo(int category,
                                             const string& pattern,
                                             double minBrightness,
                                             double maxBrightness) {
  vector<string> names =
      getPatches(category, pattern, minBrightness, maxBrightness);
  vector<PatchInfo> info;
  info.reserve(names.size());
  for (vector<string>::const_iterator iter = names.begin(); iter != names.end();
       ++iter) {
    info.push_back(itsPatches.find(*iter)->second);
  }
  return info;
}

vector<SourceInfo> SourceDBBlob::getPatchSources(const string& patchName) {
  // If not done yet, read all data from the file.
  readAll();
  vector<SourceInfo> info;
  map<string, vector<SourceData>>::const_iterator iter =
      itsSources.find(patchName);
  if (iter != itsSources.end()) {
    const vector<SourceData>& sources = iter->second;
    info.reserve(sources.size());
    for (vector<SourceData>::const_iterator srciter = sources.begin();
         srciter != sources.end(); ++srciter) {
      info.push_back(srciter->getInfo());
    }
  }
  return info;
}

vector<SourceData> SourceDBBlob::getPatchSourceData(const string& patchName) {
  // If not done yet, read all data from the file.
  readAll();
  return itsSources[patchName];
}

SourceInfo SourceDBBlob::getSource(const string&) {
  throw std::runtime_error("SourceDBBlob::getSource not implemented");
}

vector<SourceInfo> SourceDBBlob::getSources(const string&) {
  throw std::runtime_error("SourceDBBlob::getSources not implemented");
  return vector<SourceInfo>();
}

bool SourceDBBlob::atEnd() { return itsFile.tellg() >= itsEndPos; }

void SourceDBBlob::rewind() { itsFile.seekg(0, ios::beg); }

void SourceDBBlob::getNextSource(SourceData& src) {
  while (true) {
    string type = itsBlobIn->getNextType();
    if (type == "source") {
      break;
    }
    // Skip the patch info.
    PatchInfo info;
    *itsBlobIn >> info;
  }
  src.readSource(*itsBlobIn);
}

void SourceDBBlob::readAll() {
  // Only read if not done yet.
  if (itsSources.size() > 0) {
    return;
  }
  // Keep current position.
  int64_t pos = itsFile.tellg();
  // Read from beginning till end.
  rewind();
  while (itsFile.tellg() < itsEndPos) {
    if (itsBlobIn->getNextType() == "patch") {
      PatchInfo info;
      *itsBlobIn >> info;
      itsPatches.insert(make_pair(info.getName(), info));
    } else {
      SourceData info;
      info.readSource(*itsBlobIn);
      itsSources[info.getPatchName()].push_back(info);
    }
  }
  // Reset original position.
  itsFile.seekp(pos);
}

}  // namespace parmdb
}  // namespace dp3
