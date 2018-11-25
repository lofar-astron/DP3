//# ParmFacadeRep.h: Data access the parameter database
//#
//# Copyright (C) 2006
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
//# $Id: ParmFacadeRep.h 27639 2013-12-04 08:02:12Z diepen $

#ifndef LOFAR_PARMDB_PARMFACADEREP_H
#define LOFAR_PARMDB_PARMFACADEREP_H

// \file
// Data access the parameter database.

//# Includes
#include "ParmDB.h"

#include <casacore/casa/Containers/Record.h>

namespace DP3 { namespace BBS {


  // \ingroup ParmDB
  // @{

  // ParmFacadeRep is the high level interface to the Parameter Data Base.
  // The current version assumes it is an AIPS++ table; with a few extra
  // constructor arguments it can easily be changed to other types of
  // databases.
  //
  // The class provides a few functions:
  // <ul>
  // <li> getNames returns a vector of the parameter names in the table
  // <li> getRange returns a vector of 4 elements giving the boundary box
  //    of the domains of the parameters.
  // <li> getValues returns the values of the parameters a calculated on
  //      a grid given by the caller.
  // </ul>
  //
  // The parameter names can be given as a pattern. This is the same as a
  // file name pattern that can be given in the UNIX shells (e.g. RA:*).
  // Thus it is not a full regular expression.

  class ParmFacadeRep
  {
  public:
    // Define a shared_ptr for this class.
    typedef std::shared_ptr<ParmFacadeRep> ShPtr;

    // The destructor disconnects.
    virtual ~ParmFacadeRep();

    // Get the domain range (as startfreq,endfreq,starttime,endtime)
    // of the given parameters in the table.
    // This is the minimum start value and maximum end value for all parameters.
    // An empty name pattern is the same as * (all parm names).
    virtual std::vector<double> getRange (const string& parmNamePattern) const = 0;

    // Get parameter names in the table matching the pattern.
    // An empty name pattern is the same as * (all parm names).
    virtual std::vector<string> getNames (const string& parmNamePattern,
                                     bool includeDefaults) const = 0;

    // Get default parameter names matching the pattern.
    // An empty name pattern is the same as * (all parm names).
    virtual std::vector<string> getDefNames (const string& parmNamePattern) const = 0;

    // Get the default values of parameters matching the pattern.
    virtual casacore::Record getDefValues (const string& parmNamePattern) const = 0;

    // Add one or more default values.
    // The name of each field in the record is the parameter name.
    // The values are subrecords containing the parameter values, etc.
    // <br>By default it checks if the name does not exist.
    virtual void addDefValues (const casacore::Record&, bool check=true) = 0;

    // Delete the default value records for the given parameters.
    virtual void deleteDefValues (const string& parmNamePattern) = 0;

    // Get the values of the given parameters on the given regular grid
    // where v1/v2 represents center/width or start/end.
    // The Record contains a map of parameter name to Array<double>.
    virtual casacore::Record getValues (const string& parmNamePattern,
                                    double freqv1, double freqv2,
                                    double freqStep,
                                    double timev1, double timev2,
                                    double timeStep,
                                    bool asStartEnd,
                                    bool includeDefaults) = 0;

    // Get the values of the given parameters on the given grid where v1/v2
    // represents center/width or start/end.
    // The Record contains a map of parameter name to Array<double>.
    virtual casacore::Record getValues (const string& parmNamePattern,
                                    const std::vector<double>& freqv1,
                                    const std::vector<double>& freqv2,
                                    const std::vector<double>& timev1,
                                    const std::vector<double>& timev2,
                                    bool asStartEnd,
                                    bool includeDefaults) = 0;

    // Get the values of the given parameters for the given domain.
    // The Record contains a map of parameter name to Array<value>.
    // Furthermore it contains a subrecord "_grid" containing the grid axes
    // used for each parameters. Their names have the form <parmname>/xx
    // where xx is freqs, freqwidths, times, and timewidths. Their values
    // are the center and width of each cell.
    virtual casacore::Record getValuesGrid (const string& parmNamePattern,
                                        double freqv1, double freqv2,
                                        double timev1, double timev2,
                                        bool asStartEnd) = 0;

    // Get coefficients, errors, and domains they belong to.
    virtual casacore::Record getCoeff (const string& parmNamePattern,
                                   double freqv1, double freqv2,
                                   double timev1, double timev2,
                                   bool asStartEnd) = 0;

    // Flush possible changes to disk.
    virtual void flush (bool fsync) = 0;

    // Writelock and unlock the database tables.
    // The user does not need to lock/unlock, but it can increase performance
    // if many small accesses have to be done.
    // <group>
    virtual void lock (bool lockForWrite) = 0;
    virtual void unlock() = 0;
    // </group>

    // Clear the tables, thus remove all parameter values and default values.
    virtual void clearTables() = 0;

    // Set the default step values.
    virtual void setDefaultSteps (const std::vector<double>&) = 0;

    // Delete the records for the given parameters and domain.
    virtual void deleteValues (const string& parmNamePattern,
                               double freqv1, double freqv2,
                               double timev1, double timev2,
                               bool asStartEnd) = 0;

    // The following functions are only implemented for a local ParmDB.
    // The ParmFacadeDistr functions throw an exception.

    // Get the default step values for the axes.
    virtual std::vector<double> getDefaultSteps() const = 0;

    // Add the values for the given parameter names and domain.
    // The name of each field in the record is the parameter name.
    // The values are subrecords containing the domains, parameter values, etc.
    // <br>It checks if no values exist for the parameters and domains yet.
    virtual void addValues (const casacore::Record&) = 0;
  };

  // @}

}} // namespaces

#endif
