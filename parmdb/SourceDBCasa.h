// SourceDBCasa.h: Class for a Casa table holding sources and their parameters
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Base class for a table holding sources and their parameters
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_SOURCEDBCASA_H
#define LOFAR_PARMDB_SOURCEDBCASA_H

#include "SourceDB.h"
#include "PatchInfo.h"

#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/version.h>

#include "../common/Types.h"

#include <set>

namespace dp3 {
namespace parmdb {

/// @ingroup ParmDB
/// @{

/// @brief Class for a Casa table holding source parameters.
class SourceDBCasa : public SourceDBRep {
 public:
  SourceDBCasa(const ParmDBMeta& pdm, bool forceNew);

  virtual ~SourceDBCasa();

  /// Writelock and unlock the database tables.
  /// The user does not need to lock/unlock, but it can increase performance
  /// if many small accesses have to be done.
  /// The default implementation does nothing.
  ///@{
  virtual void lock(bool lockForWrite);
  virtual void unlock();
  ///@}

  /// Check for duplicate patches or sources.
  /// An exception is thrown if that is the case.
  virtual void checkDuplicates();

  /// Find non-unique patch names.
  virtual std::vector<string> findDuplicatePatches();

  /// Find non-unique source names.
  virtual std::vector<string> findDuplicateSources();

  /// Test if the patch already exists.
  virtual bool patchExists(const string& patchName);

  /// Test if the source already exists.
  virtual bool sourceExists(const string& sourceName);

  /// Add a patch and return its patchId.
  /// Nomally ra and dec should be filled in, but for moving patches
  /// (e.g. sun) this is not needed.
  /// <br>Optionally it is checked if the patch already exists.
  virtual unsigned int addPatch(const string& patchName, int catType,
                                double apparentBrightness, double ra,
                                double dec, bool check);

  /// Update the ra/dec and apparent brightness of a patch.
  virtual void updatePatch(unsigned int patchId, double apparentBrightness,
                           double ra, double dec);

  /// Add a source to a patch.
  /// Its ra and dec and default parameters will be stored as default
  /// values in the associated ParmDB tables. The names of the parameters
  /// will be preceeded by the source name and a colon.
  /// The map should contain the parameters belonging to the source type.
  /// Missing parameters will default to 0.
  /// <br>Optionally it is checked if the patch already exists.
  ///@{
  virtual void addSource(const SourceInfo& sourceInfo, const string& patchName,
                         const ParmMap& defaultParameters, double ra,
                         double dec, bool check);
  virtual void addSource(const SourceData& source, bool check);
  ///@}

  /// Add a source which forms a patch in itself (with the same name).
  /// <br>Optionally it is checked if the patch or source already exists.
  virtual void addSource(const SourceInfo& sourceInfo, const string& patchName,
                         int catType, double apparentBrightness,
                         const ParmMap& defaultParameters, double ra,
                         double dec, bool check);

  /// Get patch names in order of category and decreasing apparent flux.
  /// category < 0 means all categories.
  /// A brightness < 0 means no test on brightness.
  virtual std::vector<string> getPatches(int category, const string& pattern,
                                         double minBrightness,
                                         double maxBrightness);

  /// Get the info of all patches (name, ra, dec).
  virtual std::vector<PatchInfo> getPatchInfo(int category,
                                              const string& pattern,
                                              double minBrightness,
                                              double maxBrightness);

  /// Get the sources belonging to the given patch.
  virtual std::vector<SourceInfo> getPatchSources(const string& patchName);

  /// Get all data of the sources belonging to the given patch.
  virtual std::vector<SourceData> getPatchSourceData(const string& patchName);

  /// Get the source info of the given source.
  virtual SourceInfo getSource(const string& sourceName);

  /// Get the info of all sources matching the given (filename like) pattern.
  virtual std::vector<SourceInfo> getSources(const string& pattern);

  /// Delete the sources records matching the given (filename like) pattern.
  virtual void deleteSources(const std::string& sourceNamePattern);

  /// Clear database or table
  virtual void clearTables();

  /// Get the next source from the table.
  /// An exception is thrown if there are no more sources.
  virtual void getNextSource(SourceData& src);

  /// Tell if we are the end of the file.
  virtual bool atEnd();

  /// Reset to the beginning of the file.
  virtual void rewind();

 private:
  /// Create the source and patch table.
  void createTables(const string& tableName);

  /// Add a source for the given patch.
  void addSrc(const SourceInfo& sourceInfo, unsigned int patchId,
              const ParmMap& defaultParameters, double ra, double dec);

  /// Write the patch info in the given row.
  void writePatch(double apparentBrightness, double ra, double dec,
                  unsigned int rownr);

  /// Find the duplicate patches or sources.
  std::vector<string> findDuplicates(casacore::Table& table,
                                     const string& columnName);

  /// Fill the patch and source set object from the tables.
  /// They serve as a cache to find out if a patch or source name exists.
  void fillSets();

  /// Read all sources from the table and return them as a vector.
  std::vector<SourceInfo> readSources(const casacore::Table& table);

  /// Create the patches subset matching the given arguments.
  casacore::Table selectPatches(int category, const string& pattern,
                                double minBrightness,
                                double maxBrightness) const;

  /// Read a default parameter.
  double getDefaultParmValue(const string& name);

  casacore::Table itsPatchTable;
  casacore::Table itsSourceTable;
  std::set<std::string> itsPatchSet;
  std::set<std::string> itsSourceSet;
  bool itsSetsFilled;
  casacore::Vector<common::rownr_t> itsRowNr;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
