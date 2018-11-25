//# ParmDB.h: Base class for a table holding parameters
//#
//# Copyright (C) 2002
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
//# $Id: ParmDB.h 20771 2012-04-19 12:04:48Z diepen $

// @file
// @brief Base class for a table holding parameters
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMDB_H
#define LOFAR_PARMDB_PARMDB_H

//# Includes
#include "ParmMap.h"
#include "ParmDBMeta.h"
#include "ParmSet.h"

#include <vector>
#include <map>


namespace DP3 {
namespace BBS {

  // @ingroup ParmDB
  // @{

  // @brief Abstract base class for a table holding parameters.
  class ParmDBRep
  {
  public:
    ParmDBRep();

    virtual ~ParmDBRep();

    // Link to the DBRep by incrementing the count.
    void link()
      { itsCount++; }

    // Unlink by decrementing the count.
    int unlink()
      { return --itsCount; }

    // Flush possible changes to disk.
    // <br>If <src>fsync=True</src> the file contents are fsync-ed to disk,
    // to ensure that the system buffers are actually written to disk.
    // The default implementation does nothing.
    virtual void flush (bool fsync);

    // Writelock and unlock the database tables.
    // The user does not need to lock/unlock, but it can increase performance
    // if many small accesses have to be done.
    // The default implementation does nothing.
    // <group>
    virtual void lock (bool lockForWrite);
    virtual void unlock();
    // </group>

    // Get the domain range (freq,time) of the given parameters in the table.
    // This is the minimum and maximum value of these axes for all parameters.
    // An empty name pattern is the same as * (all parms).
    // <group>
    virtual Box getRange (const std::string& parmNamePattern) const = 0;
    virtual Box getRange (const std::vector<std::string>& parmNames) const = 0;
    // </group>

    // Get the default step values for the axes.
    const std::vector<double>& getDefaultSteps() const
      { return itsDefSteps; }

    // Set the default step values.
    virtual void setDefaultSteps (const std::vector<double>&) = 0;

    // Get the parameter values for the given parameters and domain.
    // Only * and ? should be used in the pattern (no [] and {}).
    // The default implementation uses the following getValues function.
    virtual void getValuesPattern (ParmMap& result,
                                   const std::string& parmNamePattern,
                                   const Box& domain);

    // Get the parameter values for the given parameters and domain.
    // The parmids form the indices in the result vector.
    virtual void getValues (std::vector<ParmValueSet>& values,
                            const std::vector<uint>& nameIds,
                            const std::vector<ParmId>& parmIds,
                            const Box& domain) = 0;

    // Put the values for the given parameter name and id.
    // If it is a new value, the new rowid will be stored in the ParmValueSet.
    // If it is a new name, the nameId will be filled in.
    virtual void putValues (const string& parmName, int& nameId,
                            ParmValueSet& values) = 0;

    // Put the value for the given parameters and domain.
    // It only writes the parameters that have the same DBSeqNr as this ParmDB.
    // overwriteMask=true indicates that the solvableMask might be changed
    // and should be overwritten in an existing record.
    //virtual void putValues (ParmMap& parmSet) = 0;

    // Delete the records for the given parameters and domain.
    virtual void deleteValues (const std::string& parmNamePattern,
                               const Box& domain) = 0;

    // Get the default value for the given parameter.
    // If no default value is defined in the ParmDB, the given default value
    // will be used.
    ParmValueSet getDefValue (const std::string& parmName,
                              const ParmValue& defaultValue);

    // Get the default value for the given parameters.
    // Only * and ? should be used in the pattern (no [] and {}).
    virtual void getDefValues (ParmMap& result,
                               const std::string& parmNamePattern) = 0;

    // Put the default value.
    virtual void putDefValue (const string& parmName,
                              const ParmValueSet& value,
                              bool check=true) = 0;

    // Delete the default value records for the given parameters.
    virtual void deleteDefValues (const std::string& parmNamePattern) = 0;

    // Get the names of all parms matching the given (filename like) pattern.
    virtual std::vector<std::string> getNames (const std::string& pattern) = 0;

    // Get the id of a parameter.
    // If not found in the Names table, it returns -1.
    virtual int getNameId (const std::string& parmName) = 0;

    // Clear database or table
    virtual void clearTables() = 0;

    // Set or get the name and type.
    // <group>
    void setParmDBMeta (const ParmDBMeta& ptm)
      { itsPTM = ptm; }
    const ParmDBMeta& getParmDBMeta() const
      { return itsPTM; }
    // </group>

    // Set or get ParmDB sequence nr.
    // <group>
    void setParmDBSeqNr (int seqnr)
      { itsSeqNr = seqnr; }
    int getParmDBSeqNr() const
      { return itsSeqNr; }
    // </group>

    // Set the default value map to being not filled.
    // This is needed after a delete, etc.
    void clearDefFilled()
      { itsDefFilled = false; }

  protected:
    // Set the i-th default step value (i<2) in order x,y.
    void setDefStep (uint i, double value)
      { itsDefSteps[i] = value; }

  private:
    // Fill the map with default values.
    virtual void fillDefMap (ParmMap& defMap) = 0;

    int        itsCount;
    ParmDBMeta itsPTM;
    int        itsSeqNr;
    bool       itsDefFilled;
    ParmMap    itsDefValues;
    std::vector<double> itsDefSteps;
  };


  // @brief Envelope class for a table holding parameters
  class ParmDB
  {
  public:
    // Create the ParmDB object for the given database type.
    // It gets added to the map of open parmDBs.
    explicit ParmDB (const ParmDBMeta& ptm, bool forceNew=false);

    // Copy contructor has reference semantics.
    ParmDB (const ParmDB&);

    // Delete underlying object if no more references to it.
    ~ParmDB()
      { decrCount(); }

    // Assignment has reference semantics.
    ParmDB& operator= (const ParmDB&);

    // Flush possible changes to disk.
    // <br>If <src>fsync=True</src> the file contents are fsync-ed to disk,
    // to ensure that the system buffers are actually written to disk.
    void flush (bool fsync=false)
      { itsRep->flush(fsync); }

    // Lock and unlock the database tables.
    // The user does not need to lock/unlock, but it can increase performance
    // if many small accesses have to be done.
    // <group>
    void lock (bool lockForWrite = true)
      { itsRep->lock (lockForWrite); }
    void unlock()
      { itsRep->unlock(); }

    // Get the domain range (freq,time) of the given parameters in the table.
    // This is the minimum and maximum value of these axes for all parameters.
    // An empty name pattern is the same as * (all parms).
    // <group>
    Box getRange (const std::string& parmNamePattern = "") const
      { return itsRep->getRange (parmNamePattern); }
    Box getRange (const std::vector<std::string>& parmNames) const
      { return itsRep->getRange (parmNames); }
    // </group>

    // Get the default step values for the axes.
    const std::vector<double>& getDefaultSteps() const
      { return itsRep->getDefaultSteps(); }

    // Set the default step values.
    void setDefaultSteps (const std::vector<double>& steps)
      { itsRep->setDefaultSteps (steps); }

    // Get the parameter values for the given parameters and domain.
    // Only * and ? should be used in the pattern (no [] and {}).
    void getValues (ParmMap& result,
                    const std::string& parmNamePattern,
                    const Box& domain) const
      { itsRep->getValuesPattern (result, parmNamePattern, domain); }

    // Get the parameter values for the given parameters and domain.
    // The parmids form the indices in the result vector.
    void getValues (std::vector<ParmValueSet>& values,
                    const std::vector<uint>& nameIds,
                    const std::vector<ParmId>& parmIds,
                    const Box& domain)
      { itsRep->getValues (values, nameIds, parmIds, domain); }

    // Put the values of a parameter.
    // If it is a new value, the new rowid will be stored in the ParmValueSet.
    // If it is a new name, the nameId will be filled in.
    void putValues (const string& name, int& nameId,
                    ParmValueSet& values)
      { itsRep->putValues (name, nameId, values); }

    // Delete the records for the given parameters and domain.
    void deleteValues (const std::string& parmNamePattern,
                       const Box& domain)
      { itsRep->deleteValues (parmNamePattern, domain); }

    // Get the initial value for the given parameter.
    ParmValueSet getDefValue (const std::string& parmName,
                              const ParmValue& defaultValue = ParmValue()) const
      { return itsRep->getDefValue (parmName, defaultValue); }

    // Get the default value for the given parameters.
    // Only * and ? should be used in the pattern (no [] and {}).
    void getDefValues (ParmMap& result,
                       const std::string& parmNamePattern) const
      { itsRep->getDefValues (result, parmNamePattern); }

    // Put the default value for the given parameter.
    void putDefValue (const string& parmName, const ParmValueSet& value,
                      bool check=true)
      { itsRep->putDefValue (parmName, value, check); }

    // Delete the default value records for the given parameters.
    void deleteDefValues (const std::string& parmNamePattern)
      { itsRep->deleteDefValues (parmNamePattern); }

    // Get the names matching the pattern in the table.
    std::vector<std::string> getNames (const std::string& pattern) const
      { return itsRep->getNames (pattern); }

    // Get the id of a parameter.
    // If not found in the Names table, it returns -1.
    int getNameId (const std::string& parmName)
      { return itsRep->getNameId (parmName); }

    // Clear database tables (i.e. remove all rows from all tables).
    void clearTables()
      { itsRep->clearTables(); }

    // Get the name and type of the ParmDB.
    const ParmDBMeta& getParmDBMeta() const
      { return itsRep->getParmDBMeta(); }

    // Get ParmDB sequence nr.
    int getParmDBSeqNr() const
      { return itsRep->getParmDBSeqNr(); }

    // Get the ParmDB object of the opened database for the given index.
    // An exception is thrown if not found.
    static ParmDB getParmDB (uint index);

  private:
    // Create a ParmDB object for an existing ParmDBRep.
    ParmDB (ParmDBRep*);

    // Decrement the refcount and delete if zero.
    void decrCount();

    ParmDBRep* itsRep;

    //# Keep a list of all open ParmDBs.
    static std::map<std::string,int> theirDBNames;
    static std::vector<ParmDBRep*>   theirParmDBs;
  };

  // @}

} // namespace BBS
} // namespace LOFAR

#endif
