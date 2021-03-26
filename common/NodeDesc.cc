// NodeDesc.cc: Description of a node
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen <diepen AT astron nl>

#include "NodeDesc.h"

#include "StreamUtil.h"
#include "StringTools.h"

#include <ostream>

#include <boost/algorithm/string/case_conv.hpp>

using namespace std;

namespace dp3 {
namespace common {

NodeDesc::NodeDesc(const ParameterSet& parset) {
  itsName = parset.getString("NodeName");
  string type(boost::to_lower_copy(parset.getString("NodeType", "Compute")));
  if (type == "compute") {
    itsType = Compute;
  } else if (type == "storage") {
    itsType = Storage;
  } else if (type == "head") {
    itsType = Head;
  } else {
    itsType = Any;
  }
  itsMounts = parset.getStringVector("NodeMountPoints", true);
  itsFileSys = parset.getStringVector("NodeFileSys", itsMounts, true);
  if (itsFileSys.size() != itsMounts.size())
    throw std::runtime_error(
        "NodeMountPoints and NodeFileSys should have the same size");
  for (unsigned int i = 0; i < itsMounts.size(); ++i) {
    if (itsFileSys[i].size() == 0)
      throw std::runtime_error("NodeFileSys elements can't be empty");
    if (itsMounts[i].size() == 0 || itsMounts[i][0] != '/')
      throw std::runtime_error(
          "NodeMountPoints elements can't be empty and need to start with /");
  }
}

void NodeDesc::addFileSys(const string& fsName, const string& mountPoint) {
  assert(fsName.size() > 0);
  string mp(mountPoint);
  if (mp.size() > 5 && mp.substr(0, 5) == "/auto") {
    mp = mp.substr(5);
  }
  assert(mp.size() > 0 && mp[0] == '/');
  itsFileSys.push_back(fsName);
  itsMounts.push_back(mp);
}

string NodeDesc::findFileSys(const string& fileName) const {
  // The file name must be absolute.
  assert(fileName.size() > 1 && fileName[0] == '/');
  // Determine the max nr of parts in the mount point.
  // Remember the root filesys (a single /).
  int nrp = 0;
  int rootfs = -1;
  for (unsigned int i = 0; i < itsMounts.size(); ++i) {
    int nr = 0;
    const string& str = itsMounts[i];
    // A single / counts as no part.
    if (str.size() == 1) {
      rootfs = i;
    } else {
      for (unsigned int j = 0; j < str.size(); ++j) {
        if (str[j] == '/') {
          ++nr;
        }
      }
    }
    if (nr > nrp) {
      nrp = nr;
    }
  }
  // Find the slashes in the file name for each part.
  vector<int> pos(nrp, -1);
  int nr = 0;
  for (unsigned int i = 1; i < fileName.size() && nr < nrp; ++i) {
    if (fileName[i] == '/') {
      pos[nr++] = i;
    }
  }
  // Now compare if it matches the file name.
  // Start with the longest possible string.
  for (int p = nr - 1; p >= 0; --p) {
    string filePart = fileName.substr(0, pos[p]);
    for (unsigned int i = 0; i < itsMounts.size(); ++i) {
      if (filePart == itsMounts[i]) {
        return itsFileSys[i];
      }
    }
  }
  // No match, so return root file system if there.
  // Otherwise return empty string.
  if (rootfs >= 0) {
    return itsFileSys[rootfs];
  }
  return "";
}

void NodeDesc::write(ostream& os, const string& prefix) const {
  string type = "Any";
  if (itsType == Compute) {
    type = "Compute";
  } else if (itsType == Storage) {
    type = "Storage";
  } else if (itsType == Head) {
    type = "Head";
  }
  os << prefix << "NodeName = " << itsName << endl;
  os << prefix << "NodeType = " << type << endl;
  os << prefix << "NodeFileSys     = " << itsFileSys << endl;
  os << prefix << "NodeMountPoints = " << itsMounts << endl;
}

}  // namespace common
}  // namespace dp3
