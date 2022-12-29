// ParmFacadeLocal.h: Data access the parameter database
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_PARMDB_PARMFACADELOCAL_H
#define LOFAR_PARMDB_PARMFACADELOCAL_H

#include "ParmFacadeRep.h"
#include "ParmDB.h"

#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/Vector.h>

namespace dp3 {
namespace parmdb {

/// \ingroup ParmDB
/// \brief Data access the a local parameter database.

/// @{

/// ParmFacadeLocal is the high level interface to a local Parameter Data Base.
/// The current version assumes it is an AIPS++ table; with a few extra
/// constructor arguments it can easily be changed to other types of
/// databases.
///
/// The class provides a few functions:
/// <ul>
/// <li> getNames returns a vector of the parameter names in the table
/// <li> getRange returns a vector of 4 elements giving the boundary box
///    of the domains of the parameters.
/// <li> getValues returns the values of the parameters a calculated on
///      a grid given by the caller.
/// </ul>
///
/// The parameter names can be given as a pattern. This is the same as a
/// file name pattern that can be given in the UNIX shells (e.g. RA:*).
/// Thus it is not a full regular expression.

class ParmFacadeLocal : public ParmFacadeRep {
 public:
  /// Make a connection to a new or existing ParmTable.
  ParmFacadeLocal(const string& tableName, bool create = false);

  /// The destructor disconnects.
  ~ParmFacadeLocal() override;

  /// Get the domain range (as startx,endx,starty,endy) of the given
  /// parameters in the table.
  /// This is the minimum start value and maximum end value for all parameters.
  /// An empty name pattern is the same as * (all parm names).
  std::vector<double> getRange(const string& parmNamePattern) const override;

  /// Get parameter names in the table matching the pattern.
  /// An empty name pattern is the same as * (all parm names).
  std::vector<string> getNames(const string& parmNamePattern,
                               bool includeDefaults) const override;

  /// Get default parameter names matching the pattern.
  /// An empty name pattern is the same as * (all parm names).
  std::vector<string> getDefNames(const string& parmNamePattern) const override;

  /// Get the default values of parameters matching the pattern.
  casacore::Record getDefValues(const string& parmNamePattern) const override;

  /// Add one or more default values.
  void addDefValues(const casacore::Record&, bool check) override;

  /// Delete the default value records for the given parameters.
  void deleteDefValues(const string& parmNamePattern) override;

  /// Get the values of the given parameters on the given regular grid
  /// where v1/v2 represents center/width or start/end.
  /// The Record contains a map of parameter name to Array<double>.
  casacore::Record getValues(const string& parmNamePattern, double freqv1,
                             double freqv2, double freqStep, double timev1,
                             double timev2, double timeStep, bool asStartEnd,
                             bool includeDefaults) override;

  /// Get the values of the given parameters on the given grid where v1/v2
  /// represents center/width or start/end.
  /// The Record contains a map of parameter name to Array<double>.
  casacore::Record getValues(const string& parmNamePattern,
                             const std::vector<double>& freqv1,
                             const std::vector<double>& freqv2,
                             const std::vector<double>& timev1,
                             const std::vector<double>& timev2, bool asStartEnd,
                             bool includeDefaults) override;

  /// Get the values of the given parameters for the given domain.
  /// The Record contains a map of parameter name to Array<value>.
  /// Furthermore it contains a subrecord "_grid" containing the grid axes
  /// used for each parameters. Their names have the form parmname/xx
  /// where xx is freqs, freqwidths, times, and timewidths. Their values
  /// are the center and width of each cell.
  casacore::Record getValuesGrid(const string& parmNamePattern, double freqv1,
                                 double freqv2, double timev1, double timev2,
                                 bool asStartEnd) override;

  /// Get coefficients, errors, and domains they belong to.
  casacore::Record getCoeff(const string& parmNamePattern, double freqv1,
                            double freqv2, double timev1, double timev2,
                            bool asStartEnd) override;

  /// Clear the tables, thus remove all parameter values and default values.
  void clearTables() override;

  /// Flush the possible changes to disk.
  void flush(bool fsync) override;

  /// Writelock and unlock the database tables.
  /// The user does not need to lock/unlock, but it can increase performance
  /// if many small accesses have to be done.
  ///@{
  void lock(bool lockForWrite) override;
  void unlock() override;
  ///@}

  /// Get the default step values for the axes.
  std::vector<double> getDefaultSteps() const override;

  /// Set the default step values.
  void setDefaultSteps(const std::vector<double>&) override;

  /// Add the values for the given parameter names and domain.
  void addValues(const casacore::Record& rec) override;

  /// Delete the records for the given parameters and domain.
  void deleteValues(const string& parmNamePattern, double freqv1, double freqv2,
                    double timev1, double timev2, bool asStartEnd) override;

 private:
  /// Get the values for the given predict grid
  casacore::Record doGetValues(const string& parmNamePattern,
                               const Grid& predictGrid, bool includeDefaults);

  /// Get the detailed grid from this value set.
  /// Limit it to the given domain.
  Grid getGrid(const ParmValueSet& valueSet, const Box& domain);

  /// Collect funklet coeff and errors in the record.
  casacore::Record getFunkletCoeff(const ParmValueSet& pvset);

  /// Add the default value for a single parm.
  void addDefValue(const string& parmName, const casacore::Record& value,
                   bool check);

  /// Add the value for a single parm.
  void addValue(const string& parmName, const casacore::Record& value);

  /// Construct the Grid object from the record.
  Grid record2Grid(const casacore::Record& rec) const;

  /// Make a RegularAxis or OrderedAxis.
  Axis::ShPtr makeAxis(const casacore::Vector<double>& centers,
                       const casacore::Vector<double>& widths,
                       unsigned int n) const;

  /// Convert the string to a funklet type.
  int getType(const string& str) const;

  ParmDB itsPDB;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
