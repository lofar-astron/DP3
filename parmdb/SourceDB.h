// SourceDB.h: Base class for a table holding sources and their parameters
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Base class for a table holding sources and their parameters
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_SOURCEDB_H
#define LOFAR_PARMDB_SOURCEDB_H

#include "SourceData.h"
#include "PatchInfo.h"
#include "ParmDBMeta.h"
#include "ParmDB.h"

#include <boost/filesystem.hpp>

namespace dp3 {
namespace parmdb {

class ParmMap;

/// @ingroup ParmDB
/// @{

/// @brief Abstract base class for a table holding source parameters.
class SourceDBRep {
 public:
  /// This creates the underlying ParmDB object.
  SourceDBRep(const ParmDBMeta& ptm, bool forceNew);

  virtual ~SourceDBRep();

  /// Link to the DBRep by incrementing the count.
  void link() { itsCount++; }

  /// Unlink by decrementing the count.
  int unlink() { return --itsCount; }

  /// Writelock and unlock the database tables.
  /// The user does not need to lock/unlock, but it can increase performance
  /// if many small accesses have to be done.
  /// The default implementation does nothing.
  ///@{
  virtual void lock(bool lockForWrite);
  virtual void unlock();
  ///@}

  /// Get the associated ParmDB.
  ParmDB& getParmDB() { return itsParmDB; }

  /// Check for duplicate patches or sources.
  /// An exception is thrown if that is the case.
  virtual void checkDuplicates() = 0;

  /// Find non-unique patch names.
  virtual std::vector<string> findDuplicatePatches() = 0;

  /// Find non-unique source names.
  virtual std::vector<string> findDuplicateSources() = 0;

  /// Test if the patch already exists.
  virtual bool patchExists(const string& patchName) = 0;

  /// Test if the source already exists.
  virtual bool sourceExists(const string& sourceName) = 0;

  /// Add a patch and return its patchId.
  /// Nomally ra and dec should be filled in, but for moving patches
  /// (e.g. sun) this is not needed.
  /// <br>Optionally it is checked if the patch already exists.
  virtual unsigned int addPatch(const string& patchName, int catType,
                                double apparentBrightness, double ra,
                                double dec, bool check) = 0;

  /// Update the ra/dec and apparent brightness of a patch.
  virtual void updatePatch(unsigned int patchId, double apparentBrightness,
                           double ra, double dec) = 0;

  /// Add a source to a patch.
  /// Its ra and dec and default parameters will be stored as default
  /// values in the associated ParmDB tables. The names of the parameters
  /// will be succeeded by a colon and the source name.
  /// The map should contain the parameters belonging to the source type.
  /// Missing parameters will default to 0.
  /// <br>Optionally it is checked if the source already exists.
  ///@{
  virtual void addSource(const SourceInfo& sourceInfo, const string& patchName,
                         const ParmMap& defaultParameters, double ra,
                         double dec, bool check) = 0;
  virtual void addSource(const SourceData& source, bool check) = 0;
  ///@}

  /// Add a source which forms a patch in itself (with the same name).
  /// <br>Optionally it is checked if the patch or source already exists.
  virtual void addSource(const SourceInfo& sourceInfo, const string& patchName,
                         int catType, double apparentBrightness,
                         const ParmMap& defaultParameters, double ra,
                         double dec, bool check) = 0;

  /// Get patch names in order of category and decreasing apparent flux.
  /// category < 0 means all categories.
  /// A brightness < 0 means no test on brightness.
  virtual std::vector<string> getPatches(int category, const string& pattern,
                                         double minBrightness,
                                         double maxBrightness) = 0;

  /// Get the info of selected patches (default all patches).
  virtual std::vector<PatchInfo> getPatchInfo(int category,
                                              const string& pattern,
                                              double minBrightness,
                                              double maxBrightness) = 0;

  /// Get the sources belonging to the given patch.
  virtual std::vector<SourceInfo> getPatchSources(const string& patchName) = 0;

  /// Get all data of the sources belonging to the given patch.
  virtual std::vector<SourceData> getPatchSourceData(
      const string& patchName) = 0;

  /// Get the source type of the given source.
  virtual SourceInfo getSource(const string& sourceName) = 0;

  /// Get the info of all sources matching the given (filename like) pattern.
  virtual std::vector<SourceInfo> getSources(const string& pattern) = 0;

  /// Delete the sources records matching the given (filename like) pattern.
  virtual void deleteSources(const std::string& sourceNamePattern) = 0;

  /// Clear database or table
  virtual void clearTables() = 0;

  /// Get the name and type of the SourceDB.
  const ParmDBMeta& getParmDBMeta() const { return itsParmDB.getParmDBMeta(); }

  /// Get the next source from the table.
  /// An exception is thrown if there are no more sources.
  virtual void getNextSource(SourceData& src) = 0;

  /// Tell if we are the end of the file.
  virtual bool atEnd() = 0;

  /// Reset to the beginning of the file.
  virtual void rewind() = 0;

 private:
  int itsCount;
  ParmDB itsParmDB;
};

class SourceDBBase {
 public:
  virtual ~SourceDBBase() = default;

  /// Add a source to a patch.
  /// Its ra and dec and default parameters will be stored as default
  /// values in the associated ParmDB tables. The names of the parameters
  /// will be succeeded by a colon and the source name.
  /// The map should contain the parameters belonging to the source type.
  /// Not all parameters need to be present. The ParmDB classes will
  /// use a default of 0 for missing ones.
  ///@{
  virtual void addSource(const SourceInfo& source_info,
                         const string& patch_name, int cat_type,
                         double apparent_brightness,
                         const ParmMap& default_parameters, double ra,
                         double dec, bool check) = 0;

  virtual void addSource(const SourceInfo& source_info,
                         const string& patch_name,
                         const ParmMap& default_parameters, double ra,
                         double dec, bool check) = 0;
  ///@}

  /// Add a patch and return its patch_id.
  ///
  /// Optionally it is checked if the patch already exists.
  virtual unsigned addPatch(const string& patch_name, int cat_type,
                            double apparent_brightness, double ra, double dec,
                            bool check) = 0;

  virtual void updatePatch(unsigned patch_id, double apparent_brightness,
                           double ra, double dec) = 0;
};

/// @brief Envelope class for a table holding source parameters
class SourceDB final : public SourceDBBase {
 public:
  /// Create the SourceDB object for the given database type.
  /// It gets added to the map of open sourceDBs.
  /// @param mustExist if true and the file does not exist, an exception will be
  /// thrown
  /// @param forceNew if true, an empty file will be created even if the file
  /// exists.
  explicit SourceDB(const ParmDBMeta& ptm, bool mustExist, bool forceNew);

  /// Copy contructor has reference semantics.
  SourceDB(const SourceDB&);

  /// Delete underlying object if no more references to it.
  ~SourceDB();

  /// Assignment has reference semantics.
  SourceDB& operator=(const SourceDB&);

  /// Lock and unlock the database tables.
  /// The user does not need to lock/unlock, but it can increase performance
  /// if many small accesses have to be done.
  ///@{
  void lock(bool lockForWrite = true) { itsRep->lock(lockForWrite); }
  void unlock() { itsRep->unlock(); }
  ///@}

  /// Get the associated ParmDB.
  ParmDB& getParmDB() { return itsRep->getParmDB(); }

  /// Check for duplicate patches or sources.
  /// An exception is thrown if that is the case.
  void checkDuplicates() const { itsRep->checkDuplicates(); }

  /// Find non-unique patch names.
  std::vector<string> findDuplicatePatches() const {
    return itsRep->findDuplicatePatches();
  }

  /// Find non-unique source names.
  std::vector<string> findDuplicateSources() const {
    return itsRep->findDuplicateSources();
  }

  /// Test if the patch already exists.
  bool patchExists(const string& patchName) const {
    return itsRep->patchExists(patchName);
  }

  /// Test if the source already exists.
  bool sourceExists(const string& sourceName) const {
    return itsRep->sourceExists(sourceName);
  }

  /// Add a patch and return its patchId.
  /// Normally ra and dec should be filled in, but for moving patches
  /// (e.g. sun) this is not needed.
  /// <br>Optionally it is checked if the patch already exists.
  unsigned int addPatch(const string& patchName, int catType,
                        double apparentBrightness, double ra, double dec,
                        bool check) override {
    return itsRep->addPatch(patchName, catType, apparentBrightness, ra, dec,
                            check);
  }

  /// Update the ra/dec and apparent brightness of a patch.
  void updatePatch(unsigned int patchId, double apparentBrightness, double ra,
                   double dec) override {
    itsRep->updatePatch(patchId, apparentBrightness, ra, dec);
  }

  /// Add a source to a patch.
  /// Its ra and dec and default parameters will be stored as default
  /// values in the associated ParmDB tables. The names of the parameters
  /// will be succeeded by a colon and the source name.
  /// The map should contain the parameters belonging to the source type.
  /// Not all parameters need to be present. The ParmDB classes will
  /// use a default of 0 for missing ones.
  ///@{
  void addSource(const SourceInfo& sourceInfo, const string& patchName,
                 const ParmMap& defaultParameters, double ra, double dec,
                 bool check) override {
    itsRep->addSource(sourceInfo, patchName, defaultParameters, ra, dec, check);
  }
  void addSource(const SourceData& source, bool check = true) {
    itsRep->addSource(source, check);
  }
  ///@}

  /// Add a source which forms a patch in itself (with the same name).
  void addSource(const SourceInfo& sourceInfo, const string& patchName,
                 int catType, double apparentBrightness,
                 const ParmMap& defaultParameters, double ra, double dec,
                 bool check) override {
    itsRep->addSource(sourceInfo, patchName, catType, apparentBrightness,
                      defaultParameters, ra, dec, check);
  }

  /// Get patch names in order of category and decreasing apparent flux.
  /// category < 0 means all categories.
  /// A brightness < 0 means no test on brightness.
  std::vector<string> getPatches(int category = -1,
                                 const string& pattern = string(),
                                 double minBrightness = -1,
                                 double maxBrightness = -1) const {
    return itsRep->getPatches(category, pattern, minBrightness, maxBrightness);
  }

  /// Get the info of all patches (name, ra, dec).
  std::vector<PatchInfo> getPatchInfo(int category = -1,
                                      const string& pattern = string(),
                                      double minBrightness = -1,
                                      double maxBrightness = -1) const {
    return itsRep->getPatchInfo(category, pattern, minBrightness,
                                maxBrightness);
  }

  /// Get the info of the sources belonging to the given patch.
  std::vector<SourceInfo> getPatchSources(const string& patchName) const {
    return itsRep->getPatchSources(patchName);
  }

  /// Get all data of the sources belonging to the given patch.
  std::vector<SourceData> getPatchSourceData(const string& patchName) const {
    return itsRep->getPatchSourceData(patchName);
  }

  /// Get the source info of the given source.
  SourceInfo getSource(const string& sourceName) const {
    return itsRep->getSource(sourceName);
  }

  /// Get the info of all sources matching the given (filename like) pattern.
  std::vector<SourceInfo> getSources(const string& sourceNamePattern) const {
    return itsRep->getSources(sourceNamePattern);
  }

  /// Delete the sources records matching the given (filename like) pattern.
  void deleteSources(const std::string& sourceNamePattern) {
    itsRep->deleteSources(sourceNamePattern);
  }

  /// Clear database tables (i.e. remove all rows from all tables).
  void clearTables() { itsRep->clearTables(); }

  /// Get the name and type of the SourceDB.
  const ParmDBMeta& getParmDBMeta() const { return itsRep->getParmDBMeta(); }

  /// Get the next source from the table.
  /// An exception is thrown if there are no more sources.
  void getNextSource(SourceData& src) { itsRep->getNextSource(src); }

  /// Tell if we are the end of the file.
  bool atEnd() { return itsRep->atEnd(); }

  /// Reset to the beginning of the file.
  void rewind() { itsRep->rewind(); }

 private:
  /// Create a SourceDB object for an existing SourceDBRep.
  SourceDB(SourceDBRep*);

  /// Decrement the refcount and delete if zero.
  void decrCount();

  SourceDBRep* itsRep;

  /// Path to SourceDB
  boost::filesystem::path itsSourceDBPath;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
