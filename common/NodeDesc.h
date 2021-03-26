// NodeDesc.h: Description of a node in a cluster
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Description of a node in a cluster.
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_LMWCOMMON_NODEDESC_H
#define LOFAR_LMWCOMMON_NODEDESC_H

#include "ParameterHandler.h"

#include <string>
#include <vector>
#include <iosfwd>

namespace dp3 {
namespace common {

/// @ingroup LMWCommon
/// @brief Description of a node in a cluster.

/// This class holds the basic description of a node.
/// It tells the name of the node and which file systems it has access to.
///
/// Currently the information is made persistent in a LOFAR .parset file.
/// In the future it needs to use the Central Processor Resource Manager.

class NodeDesc {
 public:
  /// Define the node types.
  enum NodeType { Compute, Storage, Head, Any };

  /// Construct an empty object.
  /// By default its type is Any.
  NodeDesc() : itsType(Any) {}

  /// Construct from the given parameterset.
  explicit NodeDesc(const ParameterSet&);

  /// Set node name.
  void setName(const std::string& name) { itsName = name; }

  /// Set node type.
  void setType(NodeType type) { itsType = type; }

  /// Add a file system the node has access to.
  /// A possible leading /auto is removed from the mountPoint.
  void addFileSys(const std::string& fsName, const std::string& mountPoint);

  /// Write it in parset format.
  void write(std::ostream& os, const std::string& prefix) const;

  /// Get the name.
  const std::string& getName() const { return itsName; }

  /// Get the type.
  NodeType getType() const { return itsType; }

  /// Get the file systems it has access to.
  const std::vector<std::string>& getFileSys() const { return itsFileSys; }

  /// Get the mount points of the file systems.
  const std::vector<std::string>& getMountPoints() const { return itsMounts; }

  /// Find the file system a file is on.
  /// The file must be given with its absolute file name.
  /// It does it by comparing the mount points with the leading part
  /// of the file name.
  std::string findFileSys(const std::string& fileName) const;

 private:
  std::string itsName;  ///< full name of the node
  NodeType itsType;
  std::vector<std::string> itsFileSys;  ///< names of file systems
  std::vector<std::string> itsMounts;   ///< and their mount points
};

}  // namespace common
}  // namespace dp3

#endif
