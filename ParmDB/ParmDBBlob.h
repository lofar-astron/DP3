//# ParmDBBlob.h: Dummy class to hold parameters in a Blob
//#
//# Copyright (C) 2012
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
//# $Id: ParmDBBlob.h 21598 2012-07-16 08:07:34Z diepen $

// @file
// @brief Dummy class to hold parameters in a Blob
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMDBBLOB_H
#define LOFAR_PARMDB_PARMDBBLOB_H

//# Includes
#include "ParmDB.h"


namespace DP3 {
namespace BBS {


  // @ingroup ParmDB
  // @{

  // @brief Dummy class to hold parameters in a Blob

  // This class is only meant to ensure that Blobs can be used for SOurceDB.
  // It throws an exception as soon as it is really used.
  class ParmDBBlob : public ParmDBRep
  {
  public:
    explicit ParmDBBlob (const std::string& tableName, bool forceNew=false);

    virtual ~ParmDBBlob();

    // Flush possible changes to disk.
    // It does not do anything.
    virtual void flush (bool fsync);

    // Writelock and unlock the table.
    // They do not do anything.
    // <group>
    virtual void lock (bool lockForWrite);
    virtual void unlock();
    // </group>

    // Get the domain range (time,freq) of the given parameters in the table.
    // They throw a "not implemented" exception.
    // <group>
    virtual Box getRange (const std::string& parmNamePattern) const;
    virtual Box getRange (const std::vector<std::string>& parmNames) const;
    // </group>

    // Set the default step values.
    // It throws a "not implemented" exception.
    virtual void setDefaultSteps (const std::vector<double>&);

    // Get the parameter values for the given parameters and domain.
    // It throws a "not implemented" exception.
    virtual void getValues (std::vector<ParmValueSet>& values,
                            const std::vector<uint>& nameIds,
                            const std::vector<ParmId>& parmIds,
                            const Box& domain);

    // Put the values for the given parameter name and id.
    // It throws a "not implemented" exception.
    virtual void putValues (const string& parmName, int& nameId,
                            ParmValueSet& values);

    // Delete the value records for the given parameters and domain.
    // It throws a "not implemented" exception.
    virtual void deleteValues (const std::string& parmNamePattern,
                               const Box& domain);

    // Get the default value for the given parameters.
    // It throws a "not implemented" exception.
    virtual void getDefValues (ParmMap& result,
                               const std::string& parmNamePattern);

    // Put the default value.
    // It throws a "not implemented" exception.
    virtual void putDefValue (const string& name, const ParmValueSet& value,
                              bool check=true);

    // Delete the default value records for the given parameters.
    // It throws a "not implemented" exception.
    virtual void deleteDefValues (const std::string& parmNamePattern);

    // Get the names of all parms matching the given (filename like) pattern.
    // It throws a "not implemented" exception.
    virtual std::vector<std::string> getNames (const std::string& pattern);

    // Get the id of a parameter.
    // It throws a "not implemented" exception.
    virtual int getNameId (const std::string& parmName);

    // Clear database or table.
    // It does not do anything.
    virtual void clearTables();

    // Fill the map with default values.
    virtual void fillDefMap (ParmMap& defMap);
  };

  // @}

} // namespace BBS
} // namespace LOFAR

#endif
