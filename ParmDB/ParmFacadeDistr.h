//# ParmFacadeDistr.h: Data access a distributed parameter database
//#
//# Copyright (C) 2009
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
//# $Id: ParmFacadeDistr.h 27747 2013-12-12 11:09:58Z diepen $

#ifndef LOFAR_PARMDB_PARMFACADEDISTR_H
#define LOFAR_PARMDB_PARMFACADEDISTR_H

// \file
// Data access the a distributed parameter database.
/*
//# Includes
#include "ParmFacadeRep.h"

#include "SocketConnectionSet.h"
#include "MWBlobIO.h"

//# Forward Declaration.
namespace casacore {
  class Record;
  class String;
}

namespace LOFAR { namespace BBS {


  // \ingroup ParmDB
  // @{

  // ParmFacadeDistr is the high level interface to a distributed Parameter
  // Data Base.
  // It starts the remote processes and connects to them. At the end it
  // sends them a quit message.
  //
  // The remote processes are started via ssh using startdistproc (in LMWCommon).
  // It starts the script 'parmdbremote-scr' in the background piping its
  // stdout and stderr output to the log file parmdbremote_$USER_$$. Because it
  // is started in the background, no ssh-connection is kept open.
  // In its turn parmdbremote-scr starts the program parmdbremote which connects
  // to ParmFacadeDistr and sends an initial message containing the ParmDB name,
  // the parm names it contains and all default values.
  // Thereafter it waits for requests which can be one of the following:
  // <ul>
  //  <li> Quit: end the remote client.
  //  <li> GetRange: return a vector of 4 elements giving the boundary box
  //       of the domains of the parameters.
  //  <li> GetValues: return the values of the parms calculated on a given grid
  //  <li> GetValuesVec: as above with grid given as vectors
  //  <li> GetValuesGrid: return parm values in box using default freq/time step
  //  <li> GetCoeff: return the coefficients of polynomial parms
  //  <li> ClearTables: clear all tables (remove all values and defaults)
  //  <li> Flush: flush the ParmDB
  //  <li> Lock: lock the ParmDB
  //  <li> Unlock: unlock and flush the ParmDB
  //  <li> SetDefaultSteps: set the default freq and time step
  //  <li> AddDefValues: add one or more default parm values
  //  <li> DeleteDefValues: delete one or more default parm values
  //  <li> DeleteValues: delete one or more parm values
  // </ul>
  //
  // The parameter names can be given as a pattern. This is the same as a
  // file name pattern that can be given in the UNIX shells (e.g. RA:*).
  // Thus it is not a full regular expression.

  class ParmFacadeDistr : public ParmFacadeRep
  {
  public:
    // Define the possible commands.
    enum Command {Quit,
                  GetRange,
                  GetValues,
                  GetValuesVec,
                  GetValuesGrid,
                  GetCoeff,
                  GetDefaultSteps,
                  Flush,
                  Lock,
                  Unlock,
                  ClearTables,
                  SetDefaultSteps,
                  AddDefValues,
                  DeleteDefValues,
                  AddValues,
                  DeleteValues,
    };

    // Make a connection to the given distributed ParmTable.
    // It starts the remote processes which connect to this object.
    ParmFacadeDistr (const string& tableName);

    // The destructor disconnects and sends the remote processes an
    // end message.
    virtual ~ParmFacadeDistr();

    // Get the domain range (as startx,endx,starty,endy) of the given
    // parameters in the table.
    // This is the minimum start value and maximum end value for all parameters.
    // An empty name pattern is the same as * (all parm names).
    virtual vector<double> getRange (const string& parmNamePattern) const;

    // Get parameter names in the table matching the pattern.
    // An empty name pattern is the same as * (all parm names).
    virtual vector<string> getNames (const string& parmNamePattern,
                                     bool includeDefaults) const;

    // Get default parameter names matching the pattern.
    // An empty name pattern is the same as * (all parm names).
    virtual vector<string> getDefNames (const string& parmNamePattern) const;

    // Get the default values of parameters matching the pattern.
    virtual casacore::Record getDefValues (const string& parmNamePattern) const;

    // Add one or more default values to all underlying ParmTables.
    virtual void addDefValues (const casacore::Record&, bool check);

    // Delete the default value records for the given parameters.
    virtual void deleteDefValues (const string& parmNamePattern);

    // Get the values of the given parameters on the given regular grid
    // where v1/v2 represents center/width or start/end.
    // The Record contains a map of parameter name to Array<double>.
    virtual casacore::Record getValues (const string& parmNamePattern,
                                    double freqv1, double freqv2,
                                    double freqStep,
                                    double timev1, double timev2,
                                    double timeStep,
                                    bool asStartEnd,
                                    bool includeDefaults);

    // Get the values of the given parameters on the given grid where v1/v2
    // represents center/width or start/end.
    // The Record contains a map of parameter name to Array<double>.
    virtual casacore::Record getValues (const string& parmNamePattern,
                                    const vector<double>& freqv1,
                                    const vector<double>& freqv2,
                                    const vector<double>& timev1,
                                    const vector<double>& timev2,
                                    bool asStartEnd,
                                    bool includeDefaults);

    // Get the values of the given parameters for the given domain.
    // The Record contains a map of parameter name to Array<value>.
    // Furthermore it contains a subrecord "_grid" containing the grid axes
    // used for each parameters. Their names have the form <parmname>/xx
    // where xx is freqs, freqwidths, times, and timewidths. Their values
    // are the center and width of each cell.
    virtual casacore::Record getValuesGrid (const string& parmNamePattern,
                                        double freqv1, double freqv2,
                                        double timev1, double timev2,
                                        bool asStartEnd);

    // Get coefficients, errors, and domains they belong to.
    virtual casacore::Record getCoeff (const string& parmNamePattern,
                                   double freqv1, double freqv2,
                                   double timev1, double timev2,
                                   bool asStartEnd);

    // Flush possible changes to disk.
    virtual void flush (bool fsync);

    // Writelock and unlock the database tables.
    // The user does not need to lock/unlock, but it can increase performance
    // if many small accesses have to be done.
    // <group>
    virtual void lock (bool lockForWrite);
    virtual void unlock();
    // </group>

    // Clear the tables, thus remove all parameter values and default values.
    virtual void clearTables();

    // Set the default step values.
    virtual void setDefaultSteps (const vector<double>&);

    // Delete the records for the given parameters and domain.
    virtual void deleteValues (const string& parmNamePattern,
                               double freqv1, double freqv2,
                               double timev1, double timev2,
                               bool asStartEnd);


    // The following functions are not implemented for a distributed ParmDB
    // and throw an exception.

    // Get the default step values for the axes.
    virtual vector<double> getDefaultSteps() const;

    // Add the values for the given parameter names and domain.
    virtual void addValues (const casacore::Record& rec);

  private:
    // Send all workers a quit message.
    void quit();

    // Get and free a port.
    // <group>
    string getPort();
    void freePort();
    // </group>

    // Check the return status of a client.
    bool checkStatus (DP3CEP::MWBlobIn& bbi, int i) const;

    // Check the return status of all clients.
    void checkStatusAll() const;

    // Read a Record from the BlobStream.
    void getRecord (BlobIStream& bis, casacore::Record& rec);

    // Write a Record into the BlobStream.
    void putRecord (BlobOStream& bis, const casacore::Record& rec);

    // Check if the names of remote client inx are equal to the first one.
    void checkNames (const vector<string>& firstNames,
                     const vector<string>& names, uint inx) const;

    // Combine the result records from the remote sites.
    casacore::Record combineRemote (const vector<casacore::Record>& recs) const;

    // Find all parm names in the records and add them to the set.
    void findParmNames (const vector<casacore::Record>& recs,
                        set<casacore::String>& names) const;

    // Combine the info for the given parm from all records.
    // The info can be the same in some records meaning that a fully or partial
    // global solve is done and distributed to all parmdbs.
    void combineInfo (const casacore::String& name,
                      const vector<casacore::Record>& recs,
                      casacore::Record& result) const;

    //# Data members
    string                itsPort;      //# declare this before itsConn!!
    mutable DP3CEP::SocketConnectionSet itsConn;
    vector<string>        itsPartNames;
    vector<string>        itsParmNames;
    casacore::Record          itsDefValues;
    static int            theirNextPort;
    static vector<string> theirFreePorts;
  };

  // @}
}

} // namespaces
*/

#endif
