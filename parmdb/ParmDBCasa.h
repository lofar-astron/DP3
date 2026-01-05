// ParmDBCasa.h: Class to hold parameters in a Casa table
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Class to hold parameters in a Casa table
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMDBCASA_H
#define LOFAR_PARMDB_PARMDBCASA_H

#include "ParmDB.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include "common/Types.h"

namespace dp3 {
namespace parmdb {

/// @ingroup ParmDB
/// @{

/// @brief Class to hold parameters in a Casa table
class ParmDBCasa : public ParmDBRep {
 public:
  explicit ParmDBCasa(const std::string& tableName, bool forceNew = false);

  ~ParmDBCasa() override;

  /// Flush possible changes to disk.
  void flush(bool fsync) override;

  /// Writelock and unlock the table.
  /// It is not necessary to do this, but it can be useful if many
  /// small accesses have to be done.
  ///@{
  void lock(bool lockForWrite) override;
  void unlock() override;
  ///@}

  /// Get the domain range (time,freq) of the given parameters in the table.
  /// This is the minimum and maximum value of these axes for all parameters.
  /// An empty name pattern is the same as * (all parms).
  ///@{
  Box getRange(const std::string& parmNamePattern) const override;
  Box getRange(const std::vector<std::string>& parmNames) const override;
  ///@}

  /// Set the default step values.
  void setDefaultSteps(const std::vector<double>&) override;

  /// Get the parameter values for the given parameters and domain.
  /// The parmids form the indices in the result vector.
  void getValues(std::vector<ParmValueSet>& values,
                 const std::vector<unsigned int>& nameIds,
                 const std::vector<ParmId>& parmIds,
                 const Box& domain) override;

  /// Put the values for the given parameter name and id.
  /// If it is a new value, the new rowid will be stored in the ParmValueSet.
  /// If it is a new name, the nameId will be filled in.
  void putValues(const std::string& parmName, int& nameId,
                 ParmValueSet& values) override;

  /// Delete the value records for the given parameters and domain.
  void deleteValues(const std::string& parmNamePattern,
                    const Box& domain) override;

  /// Get the default value for the given parameters.
  /// Only * and ? should be used in the pattern (no [] and {}).
  void getDefValues(ParmMap& result,
                    const std::string& parmNamePattern) override;

  /// Put the default value.
  void putDefValue(const std::string& name, const ParmValueSet& value,
                   bool check = true) override;

  /// Delete the default value records for the given parameters.
  void deleteDefValues(const std::string& parmNamePattern) override;

  /// Get the names of all parms matching the given (filename like) pattern.
  std::vector<std::string> getNames(const std::string& pattern) override;

  /// Get the id of a parameter.
  /// If not found in the Names table, it returns -1.
  int getNameId(const std::string& parmName) override;

  /// Clear database or table
  void clearTables() override;

 private:
  /// Fill the map with default values.
  void fillDefMap(ParmMap& defMap) override;

  /// Create a parmtable with the given name.
  void createTables(const std::string& tableName);

  /// If not empty, put the domain into the table.
  /// The DOMAIN column is added to the table if it does not exist.
  void putDefDomain(const Box& domain, casacore::Table& tab,
                    unsigned int rownr);

  /// If defined in the table, set the scale domain in the ParmValue.
  Box getDefDomain(const casacore::Table& tab, unsigned int row);

  /// Get a selection from the NAME table.
  ///@{
  casacore::Table getNameSel(const std::string& parmNamePattern) const;
  casacore::Vector<common::rownr_t> getNameIds(
      const std::string& parmNamePattern) const;
  casacore::Vector<common::rownr_t> getNameIds(
      const std::vector<std::string>& parmNames) const;
  ///@}

  /// Find the minmax range in the table.
  Box findRange(const casacore::Table& table) const;

  /// Extract the parm values from a table selection with a single parm name.
  ///@{
  // void extractValues (ParmMap& result, const casacore::Table& table);
  std::pair<string, ParmValueSet> extractDefValue(const casacore::Table& sel,
                                                  int row);
  ///@}

  /// Do the actual put of a value.
  void doPutValue(const std::string& parmName, int& nameId,
                  ParmValueSet& parmSet);

  /// Put the value for an existing parameter/domain.
  void putOldValue(const ParmValue& parmValue, ParmValue::FunkletType type);

  /// Put the value for a new parameter/domain.
  void putNewValue(const std::string& name, int& nameId, ParmValueSet& parmSet,
                   ParmValue& parmValue, const Box& domain);

  /// Put an entry into the NAME table.
  int putName(const std::string& name, const ParmValueSet& pset);

  /// Put the value for a new default parameter.
  void putNewDefValue(const std::string& parmName, const ParmValueSet& value);

  /// Put the begin/end of an irregular axis.
  void putInterval(const Axis& axis, casacore::ArrayColumn<double>& col,
                   unsigned int rownr);

  /// Form an axis from the interval array in the given row.
  /// If no interval array, return a regular axis made from (st,end,n).
  Axis::ShPtr getInterval(casacore::ROArrayColumn<double>& col,
                          unsigned int rownr, double st, double end,
                          unsigned int n);

  /// Find the table subset containing the parameter values for the
  /// requested domain.
  casacore::Table find(const std::string& parmName, const Box& domain);

  /// Create a select expression node on domain.
  casacore::TableExprNode makeExpr(const casacore::Table& table,
                                   const Box& domain) const;

  /// And two table select expressions, where the first one can be null.
  void andExpr(casacore::TableExprNode& expr,
               const casacore::TableExprNode& right) const;

  casacore::Table itsTables[3];  ///< normal,names,default
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
