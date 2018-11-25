//# ParmFacadeLocal.h: Data access the parameter database
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
//# $Id: ParmFacadeLocal.h 27639 2013-12-04 08:02:12Z diepen $

#ifndef LOFAR_PARMDB_PARMFACADELOCAL_H
#define LOFAR_PARMDB_PARMFACADELOCAL_H

// \file
// Data access the a local parameter database.

//# Includes
#include "ParmFacadeRep.h"
#include "ParmDB.h"

#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/Vector.h>


namespace DP3 { namespace BBS {


  // \ingroup ParmDB
  // @{

  // ParmFacadeLocal is the high level interface to a local Parameter Data Base.
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

  class ParmFacadeLocal : public ParmFacadeRep
  {
  public:
    // Make a connection to a new or existing ParmTable.
    ParmFacadeLocal (const string& tableName, bool create=false);

    // The destructor disconnects.
    virtual ~ParmFacadeLocal();

    // Get the domain range (as startx,endx,starty,endy) of the given
    // parameters in the table.
    // This is the minimum start value and maximum end value for all parameters.
    // An empty name pattern is the same as * (all parm names).
    virtual std::vector<double> getRange (const string& parmNamePattern) const;

    // Get parameter names in the table matching the pattern.
    // An empty name pattern is the same as * (all parm names).
    virtual std::vector<string> getNames (const string& parmNamePattern,
                                     bool includeDefaults) const;

    // Get default parameter names matching the pattern.
    // An empty name pattern is the same as * (all parm names).
    virtual std::vector<string> getDefNames (const string& parmNamePattern) const;

    // Get the default values of parameters matching the pattern.
    virtual casacore::Record getDefValues (const string& parmNamePattern) const;

    // Add one or more default values.
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
                                    const std::vector<double>& freqv1,
                                    const std::vector<double>& freqv2,
                                    const std::vector<double>& timev1,
                                    const std::vector<double>& timev2,
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

    // Clear the tables, thus remove all parameter values and default values.
    virtual void clearTables();

    // Flush the possible changes to disk.
    virtual void flush (bool fsync);

    // Writelock and unlock the database tables.
    // The user does not need to lock/unlock, but it can increase performance
    // if many small accesses have to be done.
    // <group>
    virtual void lock (bool lockForWrite);
    virtual void unlock();
    // </group>

    // Get the default step values for the axes.
    virtual std::vector<double> getDefaultSteps() const;

    // Set the default step values.
    virtual void setDefaultSteps (const std::vector<double>&);

    // Add the values for the given parameter names and domain.
    virtual void addValues (const casacore::Record& rec);

    // Delete the records for the given parameters and domain.
    virtual void deleteValues (const string& parmNamePattern,
                               double freqv1, double freqv2,
                               double timev1, double timev2,
                               bool asStartEnd);

  private:
    // Get the values for the given predict grid
    casacore::Record doGetValues (const string& parmNamePattern,
                              const Grid& predictGrid,
                              bool includeDefaults);

    // Get the detailed grid from this value set.
    // Limit it to the given domain.
    Grid getGrid (const ParmValueSet& valueSet, const Box& domain);

    // Collect funklet coeff and errors in the record.
    casacore::Record getFunkletCoeff (const ParmValueSet& pvset);

    // Add the default value for a single parm.
    void addDefValue (const string& parmName, const casacore::Record& value,
                      bool check);

    // Add the value for a single parm.
    void addValue (const string& parmName, const casacore::Record& value);

    // Construct the Grid object from the record.
    Grid record2Grid (const casacore::Record& rec) const;

    // Make a RegularAxis or OrderedAxis.
    Axis::ShPtr makeAxis (const casacore::Vector<double>& centers,
                          const casacore::Vector<double>& widths, uint n) const;

    // Convert the string to a funklet type.
    int getType (const string& str) const;


    //# Data members
    ParmDB itsPDB;
  };

  // @}

}} // namespaces

#endif
