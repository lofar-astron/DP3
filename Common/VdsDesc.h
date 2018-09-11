//# VdsDesc.h: Describe an entire visibility data set
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
//# $Id: VdsDesc.h 16886 2010-12-08 10:43:17Z diepen $

#ifndef LOFAR_LMWCOMMON_VDSDESC_H
#define LOFAR_LMWCOMMON_VDSDESC_H

// @file
// @brief Describe an entire visibility data set
// @author Ger van Diepen (diepen AT astron nl)

//# Includes
#include "VdsPartDesc.h"
#include "ParameterHandler.h"

#include <casacore/casa/Utilities/Regex.h>

namespace DP3 { namespace CEP {

  // @ingroup LMWCommon
  // @brief Describe an entire visibility data set

  // This class holds the description of an entire visibility data set (VDS).
  // In VdsPartDesc objects it describes the parts it consists of and
  // on which file systems they are located.
  // A VdsPartDesc object is also used to describe the entire VDS.
  // Furthermore it contains the names of all antennae, which can be used
  // to map the antenna name to the antenna number when a selection on
  // antenna names is done.
  //
  // Currently the information is made persistent in a LOFAR .parset file.
  // In the future it needs to use the Centrol Processor Resource Manager.

  class VdsDesc
  {
  public:
    // Construct with a description of the entire visibility data set.
    // The description can be empty and set later using setDesc.
    explicit VdsDesc (const VdsPartDesc& = VdsPartDesc());

    // Construct from the given parameterset.
    // @{
    explicit VdsDesc (const std::string& parsetName);
    explicit VdsDesc (const ParameterSet& parset)
      { init (parset); }
    // @}

    // Add a part.
    void addPart (const VdsPartDesc& part)
      { itsParts.push_back (part); }

    // Get the description of the parts.
    const std::vector<VdsPartDesc>& getParts() const
      { return itsParts; }

    // Get the description of the VDS.
    const VdsPartDesc& getDesc() const
      { return itsDesc; }

    // Set the description of the VDS.
    // Usually the description is set in the constructor, but this offers
    // another way of doing it.
    void setDesc (const VdsPartDesc& desc)
      { itsDesc = desc; }

    // Write it in parset format.
    void write (std::ostream& os) const;

  private:
    // Fill the object from the given parset file.
    void init (const ParameterSet& parset);

    VdsPartDesc              itsDesc;
    std::vector<VdsPartDesc> itsParts;
  };

}} //# end namespaces

#endif
