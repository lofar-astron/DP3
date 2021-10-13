// ClusterDesc.h:  Description of a cluster and the nodes in it
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Description of a cluster and the nodes in it.
/// @author Ger van Diepen diepen AT astron nl

#ifndef LOFAR_LMWCOMMON_CLUSTERDESC_H
#define LOFAR_LMWCOMMON_CLUSTERDESC_H

#include "NodeDesc.h"

#include <string>
#include <vector>
#include <iosfwd>

namespace dp3 {
namespace common {
class ParameterSet;
}
}  // namespace dp3

namespace dp3 {
namespace common {

/// @ingroup LMWCommon
/// @brief Description of a cluster and the nodes in it.

/// This class holds the basic description of a cluster.
/// It defines which nodes are part of the cluster and which file systems
/// each node has access to.
/// If a data set is distributed over many file systems, the cluster
/// description tells which node can handle a data set part on a particular
/// file system.
///
/// Currently the information is made persistent in a LOFAR .parset file.
/// In the future it needs to use the Centrol Processor Resource Manager.

class ClusterDesc {
 public:
  /// Construct an empty object.
  ClusterDesc() {}

  /// Construct from the given parameterset.
  explicit ClusterDesc(const std::string& parsetName);

  /// Set cluster name.
  void setName(const std::string& name) { itsName = name; }

  /// Add a node description.
  /// A node with an already existing name is not added.
  /// If \p canExist=false, an exception is thrown if existing.
  void addNode(const NodeDesc& node, bool canExist = false);

  /// Write it in parset format.
  void write(std::ostream& os) const;

  /// Get the cluster name.
  const std::string& getName() const { return itsName; }

  /// Get a specific node. An exception is thrown if not found.
  const NodeDesc& getNode(const std::string& nodeName) const;

  /// Get all nodes.
  const std::vector<NodeDesc>& getNodes() const { return itsNodes; }

  /// Get the map of file system to node index.
  const std::map<std::string, std::vector<int>>& getMap() const {
    return itsFS2Nodes;
  }

 private:
  /// Fill the object from the given parset file.
  void init(const std::string& parsetName);

  /// Get the description of a homogeneous cluster.
  void getHomCluster(const ParameterSet& parset);

  /// Add nodes for a homogeneous cluster.
  void addNodes(const ParameterSet& parset, NodeDesc::NodeType type);

  /// Get the description of a heterogeneous cluster.
  void getHetCluster(const ParameterSet& parset);

  /// Fill the object from the subcluster definitions.
  /// Use the given directory for relative clusterdesc names.
  void getSubClusters(const std::vector<std::string>& parsetNames,
                      const std::string& defaultDir);

  /// Add entries to the mapping of FileSys to Nodes.
  void add2Map(int nodeIndex);

  std::string itsName;
  std::vector<NodeDesc> itsNodes;
  std::map<std::string, int> itsNodeMap;
  std::map<std::string, std::vector<int>> itsFS2Nodes;
};

}  // namespace common
}  // namespace dp3

#endif
