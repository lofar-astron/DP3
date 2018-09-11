//# ClusterDesc.h:  Description of a cluster and the nodes in it
//#
//# Copyright (C) 2005
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: ClusterDesc.h 16886 2010-12-08 10:43:17Z diepen $
//#

#ifndef LOFAR_LMWCOMMON_CLUSTERDESC_H
#define LOFAR_LMWCOMMON_CLUSTERDESC_H

// @file
// @brief Description of a cluster and the nodes in it.
// @author Ger van Diepen <diepen AT astron nl>

//# Includes
#include "NodeDesc.h"

#include <string>
#include <vector>
#include <iosfwd>

//# Forard Declarations;
namespace DP3 {
  class ParameterSet;
}

namespace DP3 { namespace CEP {

  // @ingroup LMWCommon
  // @brief Description of a cluster and the nodes in it.

  // This class holds the basic description of a cluster.
  // It defines which nodes are part of the cluster and which file systems
  // each node has access to.
  // If a data set is distributed over many file systems, the cluster
  // description tells which node can handle a data set part on a particular
  // file system.
  //
  // Currently the information is made persistent in a LOFAR .parset file.
  // In the future it needs to use the Centrol Processor Resource Manager.

  class ClusterDesc
  {
  public:
    // Construct an empty object.
    ClusterDesc()
      {}

    // Construct from the given parameterset.
    explicit ClusterDesc (const std::string& parsetName);

    // Set cluster name.
    void setName (const std::string& name)
      { itsName = name; }

    // Add a node description.
    // A node with an already existing name is not added.
    // If <src>canExist=false</src>, an exception is thrown if existing.
    void addNode (const NodeDesc& node, bool canExist=false);

    // Write it in parset format.
    void write (std::ostream& os) const;

    // Get the cluster name.
    const std::string& getName() const
      { return itsName; }

    // Get a specific node. An exception is thrown if not found.
    const NodeDesc& getNode (const std::string& nodeName) const;

    // Get all nodes.
    const std::vector<NodeDesc>& getNodes() const
      { return itsNodes; }

    // Get the map of file system to node index.
    const std::map<std::string, std::vector<int> >& getMap() const
      { return itsFS2Nodes; }

  private:
    // Fill the object from the given parset file.
    void init (const std::string& parsetName);

    // Get the description of a homogeneous cluster.
    void getHomCluster (const ParameterSet& parset);

    // Add nodes for a homogeneous cluster.
    void addNodes (const ParameterSet& parset,
                   NodeDesc::NodeType type);

    // Get the description of a heterogeneous cluster.
    void getHetCluster (const ParameterSet& parset);

    // Fill the object from the subcluster definitions.
    // Use the given directory for relative clusterdesc names.
    void getSubClusters (const std::vector<std::string>& parsetNames,
                         const std::string& defaultDir);

    // Add entries to the mapping of FileSys to Nodes.
    void add2Map (int nodeIndex);

    std::string itsName;
    std::vector<NodeDesc>      itsNodes;
    std::map<std::string, int> itsNodeMap;
    std::map<std::string, std::vector<int> > itsFS2Nodes;
  };
    
}} // end namespaces

#endif
