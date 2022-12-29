// ParmDBBlob.h: Dummy class to hold parameters in a Blob
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Dummy class to hold parameters in a Blob
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMDBBLOB_H
#define LOFAR_PARMDB_PARMDBBLOB_H

#include "ParmDB.h"

namespace dp3 {
namespace parmdb {

/// @ingroup ParmDB
/// @{

/// @brief Dummy class to hold parameters in a Blob

/// This class is only meant to ensure that Blobs can be used for SOurceDB.
/// It throws an exception as soon as it is really used.
class ParmDBBlob : public ParmDBRep {
 public:
  explicit ParmDBBlob(const std::string& tableName, bool forceNew = false);

  ~ParmDBBlob() override;

  /// Flush possible changes to disk.
  /// It does not do anything.
  void flush(bool fsync) override;

  /// Writelock and unlock the table.
  /// They do not do anything.
  ///@{
  void lock(bool lockForWrite) override;
  void unlock() override;
  ///@}

  /// Get the domain range (time,freq) of the given parameters in the table.
  /// They throw a "not implemented" exception.
  ///@{
  Box getRange(const std::string& parmNamePattern) const override;
  Box getRange(const std::vector<std::string>& parmNames) const override;
  ///@}

  /// Set the default step values.
  /// It throws a "not implemented" exception.
  void setDefaultSteps(const std::vector<double>&) override;

  /// Get the parameter values for the given parameters and domain.
  /// It throws a "not implemented" exception.
  void getValues(std::vector<ParmValueSet>& values,
                 const std::vector<unsigned int>& nameIds,
                 const std::vector<ParmId>& parmIds,
                 const Box& domain) override;

  /// Put the values for the given parameter name and id.
  /// It throws a "not implemented" exception.
  void putValues(const std::string& parmName, int& nameId,
                 ParmValueSet& values) override;

  /// Delete the value records for the given parameters and domain.
  /// It throws a "not implemented" exception.
  void deleteValues(const std::string& parmNamePattern,
                    const Box& domain) override;

  /// Get the default value for the given parameters.
  /// It throws a "not implemented" exception.
  void getDefValues(ParmMap& result,
                    const std::string& parmNamePattern) override;

  /// Put the default value.
  /// It throws a "not implemented" exception.
  void putDefValue(const std::string& name, const ParmValueSet& value,
                   bool check = true) override;

  /// Delete the default value records for the given parameters.
  /// It throws a "not implemented" exception.
  void deleteDefValues(const std::string& parmNamePattern) override;

  /// Get the names of all parms matching the given (filename like) pattern.
  /// It throws a "not implemented" exception.
  std::vector<std::string> getNames(const std::string& pattern) override;

  /// Get the id of a parameter.
  /// It throws a "not implemented" exception.
  int getNameId(const std::string& parmName) override;

  /// Clear database or table.
  /// It does not do anything.
  void clearTables() override;

  /// Fill the map with default values.
  void fillDefMap(ParmMap& defMap) override;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
