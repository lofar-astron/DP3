//# NodeDesc.h: Description of a node in a cluster
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
//# $Id: NodeDesc.h 16886 2010-12-08 10:43:17Z diepen $

#ifndef LOFAR_LMWCOMMON_NODEDESC_H
#define LOFAR_LMWCOMMON_NODEDESC_H

// @file
// @brief Description of a node in a cluster.
// @author Ger van Diepen (diepen AT astron nl)

//# Includes
#include "ParameterHandler.h"

#include <string>
#include <vector>
#include <iosfwd>

namespace DP3 { namespace CEP {

  // @ingroup LMWCommon
  // @brief Description of a node in a cluster.

  // This class holds the basic description of a node.
  // It tells the name of the node and which file systems it has access to.
  //
  // Currently the information is made persistent in a LOFAR .parset file.
  // In the future it needs to use the Central Processor Resource Manager.

  class NodeDesc
  {
  public:
    // Define the node types.
    enum NodeType {
      Compute,
      Storage,
      Head,
      Any
    };

    // Construct an empty object.
    // By default its type is Any.
    NodeDesc()
      : itsType(Any) {}

    // Construct from the given parameterset.
    explicit NodeDesc (const ParameterSet&);
 
    // Set node name.
    void setName (const std::string& name)
      { itsName = name; }

    // Set node type.
    void setType (NodeType type)
      { itsType = type; }

   // Add a file system the node has access to.
    // A possible leading /auto is removed from the mountPoint.
    void addFileSys (const std::string& fsName, const std::string& mountPoint);

    // Write it in parset format.
    void write (std::ostream& os, const std::string& prefix) const;

    // Get the name.
    const std::string& getName() const
      { return itsName; }

    // Get the type.
    NodeType getType() const
      { return itsType; }

    // Get the file systems it has access to.
    const std::vector<std::string>& getFileSys() const
      { return itsFileSys; }

    // Get the mount points of the file systems.
    const std::vector<std::string>& getMountPoints() const
      { return itsMounts; }

    // Find the file system a file is on.
    // The file must be given with its absolute file name.
    // It does it by comparing the mount points with the leading part
    // of the file name.
    std::string findFileSys (const std::string& fileName) const;

  private:
    std::string itsName;                  //# full name of the node
    NodeType    itsType;
    std::vector<std::string> itsFileSys;  //# names of file systems
    std::vector<std::string> itsMounts;   //# and their mount points
  };
    
}} //# end namespaces

#endif
