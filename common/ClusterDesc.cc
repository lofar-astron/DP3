// ClusterDesc.cc: Description of a cluster
//
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen <diepen AT astron nl>

#include "ClusterDesc.h"

#include <casacore/casa/OS/Path.h>

#include "ParameterSet.h"

using namespace std;
using namespace casacore;

namespace dp3 {
namespace common {

ClusterDesc::ClusterDesc(const std::string& parsetName) { init(parsetName); }

void ClusterDesc::init(const std::string& parsetName) {
  // Get absolute file name (it expands possible ~ and $ in the file name).
  String fullName = Path(parsetName).absoluteName();
  ParameterSet parset(fullName);
  itsName = parset.getString("ClusterName");
  if (parset.isDefined("Node0.NodeName")) {
    // The cluster can be heterogeneous and is described in detail.
    getHetCluster(parset);
  } else if (parset.isDefined("SubClusters")) {
    // Get subclusters; use parent's directory as default directory.
    getSubClusters(parset.getStringVector("SubClusters", true),
                   Path(fullName).dirName());
  } else {
    // The cluster is homogeneous and is described in a concise way.
    getHomCluster(parset);
  }
}

void ClusterDesc::getHetCluster(const ParameterSet& parset) {
  // Iterate sequentially until no node with that number is found.
  int i = 0;
  while (true) {
    ostringstream prefix;
    prefix << "Node" << i << '.';
    if (!parset.isDefined(prefix.str() + "NodeName")) {
      break;
    }
    ParameterSet subset = parset.makeSubset(prefix.str());
    NodeDesc node(subset);
    addNode(node);
    ++i;
  }
}

void ClusterDesc::getHomCluster(const ParameterSet& parset) {
  vector<std::string> defVal;
  // Add the different kind of nodes.
  addNodes(parset.makeSubset("Compute."), NodeDesc::Compute);
  addNodes(parset.makeSubset("Storage."), NodeDesc::Storage);
  addNodes(parset.makeSubset("Head."), NodeDesc::Head);
  if (parset.isDefined("Nodes")) {
    addNodes(parset, NodeDesc::Any);
  }
}

void ClusterDesc::addNodes(const ParameterSet& parset,
                           NodeDesc::NodeType type) {
  if (parset.empty()) {
    return;
  }
  vector<std::string> defVal;
  vector<std::string> names = parset.getStringVector("Nodes", true);
  vector<std::string> localDisks =
      parset.getStringVector("LocalDisks", defVal, true);
  vector<std::string> remoteDisks =
      parset.getStringVector("RemoteDisks", defVal, true);
  vector<std::string> remoteFilesys =
      parset.getStringVector("RemoteFileSys", defVal, true);
  if (remoteFilesys.empty()) {
    remoteFilesys = remoteDisks;
  } else {
    if (remoteFilesys.size() != remoteDisks.size())
      throw std::runtime_error(
          "RemoteFileSys must be empty or have same length as RemoteDisks");
  }
  for (unsigned int i = 0; i < names.size(); ++i) {
    vector<std::string> rdisks, rfilesys, lfilesys, ldisks;
    const vector<std::string>* ldiskp = &localDisks;
    const vector<std::string>* rdiskp = &remoteDisks;
    const vector<std::string>* rfsysp = &remoteFilesys;
    // It is possible to override disk specifications on a per node basis.
    string key = "LocalDisks." + names[i];
    if (parset.isDefined(key)) {
      ldisks = parset.getStringVector(key, true);
      ldiskp = &ldisks;
    }
    key = string("RemoteDisks." + names[i]);
    if (parset.isDefined(key)) {
      rdisks = parset.getStringVector(key, true);
      rdiskp = &rdisks;
      rfilesys =
          parset.getStringVector("RemoteFileSys." + names[i], defVal, true);
      if (rfilesys.empty()) {
        rfilesys = rdisks;
      } else {
        if (rfilesys.size() != rdisks.size())
          throw std::runtime_error(
              "RemoteFileSys," + names[i] +
              " must be empty or have same length as RemoteDisks." + names[i]);
      }
      rfsysp = &rfilesys;
    }
    NodeDesc node;
    node.setName(names[i]);
    node.setType(type);
    for (unsigned int j = 0; j < rdiskp->size(); ++j) {
      node.addFileSys((*rfsysp)[j], (*rdiskp)[j]);
    }
    // Add node name to local filesys to make it unique.
    for (unsigned int j = 0; j < ldiskp->size(); ++j) {
      node.addFileSys(names[i] + ':' + (*ldiskp)[j], (*ldiskp)[j]);
    }
    addNode(node);
  }
}

void ClusterDesc::getSubClusters(const vector<std::string>& parsetNames,
                                 const std::string& defaultDir) {
  for (unsigned int i = 0; i < parsetNames.size(); ++i) {
    // Expand possible ~ and $.
    string name = Path(parsetNames[i]).expandedName();
    // Add directory of parent parset if name is not absolute.
    if (name[0] != '/') {
      name = defaultDir + '/' + name;
    }
    ClusterDesc cdesc(name);
    const vector<NodeDesc>& nodes = cdesc.getNodes();
    for (unsigned int j = 0; j < nodes.size(); ++j) {
      // The same nodes can occur in multiple subclusters.
      addNode(nodes[j], true);
    }
  }
}

void ClusterDesc::write(ostream& os) const {
  os << "ClusterName = " << itsName << endl;
  os << "NNodes = " << itsNodes.size() << endl;
  for (unsigned i = 0; i < itsNodes.size(); ++i) {
    ostringstream prefix;
    prefix << "Node" << i << '.';
    itsNodes[i].write(os, prefix.str());
  }
}

void ClusterDesc::addNode(const NodeDesc& node, bool canExist) {
  // Check if node name does not exist yet.
  map<std::string, int>::const_iterator loc = itsNodeMap.find(node.getName());
  if (loc != itsNodeMap.end()) {
    if (!canExist) {
      throw std::runtime_error("Node name " + node.getName() +
                               " multiply specified in clusterdesc " + itsName);
    }
  } else {
    int inx = itsNodes.size();
    itsNodeMap[node.getName()] = inx;
    itsNodes.push_back(node);
    add2Map(inx);
  }
}

void ClusterDesc::add2Map(int nodeIndex) {
  const NodeDesc& node = itsNodes[nodeIndex];
  for (vector<std::string>::const_iterator iter = node.getFileSys().begin();
       iter != node.getFileSys().end(); ++iter) {
    vector<int>& vec = itsFS2Nodes[*iter];
    vec.push_back(nodeIndex);
  }
}

const NodeDesc& ClusterDesc::getNode(const std::string& nodeName) const {
  map<std::string, int>::const_iterator loc = itsNodeMap.find(nodeName);
  if (loc == itsNodeMap.end()) {
    throw std::runtime_error("Node name " + nodeName +
                             " not found in clusterdesc " + itsName);
  }
  return itsNodes[loc->second];
}

//   string ClusterDesc::findNode (const std::string& fileSystem,
// 				const map<std::string,int>& done) const
//   {
//     map<std::string,vector<std::string> >::const_iterator iter =
//                                              itsFS2Nodes.find(fileSystem);
//     if (iter == itsFS2Nodes.end()) {
//       return "";
//     }
//     const vector<std::string>& nodes = iter->second;
//     for (unsigned i=0; i<nodes.size(); ++i) {
//       if (done.find(nodes[i]) == done.end()) {
// 	return nodes[i];
//       }
//     }
//     return "";
//   }

}  // namespace common
}  // namespace dp3
