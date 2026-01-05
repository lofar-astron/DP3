// SourceDBBlob.h: Class for a Blob file holding sources and their parameters
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Base class for a table holding sources and their parameters
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_SOURCEDBBLOB_H
#define LOFAR_PARMDB_SOURCEDBBLOB_H

#include "SourceDB.h"

#include <fstream>
#include "blob/BlobOBufStream.h"
#include "blob/BlobIBufStream.h"
#include "blob/BlobOStream.h"
#include "blob/BlobIStream.h"

namespace dp3 {
namespace parmdb {

/// @ingroup ParmDB
/// @{

/// @brief Class for a Blob file holding source parameters.
class SourceDBBlob : public SourceDBRep {
 public:
  SourceDBBlob(const ParmDBMeta& pdm, bool forceNew);

  ~SourceDBBlob() override;

  /// Writelock and unlock the file.
  /// It does not do anything.
  ///@{
  void lock(bool lockForWrite) override;
  void unlock() override;
  ///@}

  /// Check for duplicate patches or sources.
  /// An exception is thrown if that is the case.
  void checkDuplicates() override;

  /// Find non-unique patch names.
  std::vector<std::string> findDuplicatePatches() override;

  /// Find non-unique source names.
  std::vector<std::string> findDuplicateSources() override;

  /// Test if the patch already exists.
  bool patchExists(const std::string& patchName) override;

  /// Test if the source already exists.
  bool sourceExists(const std::string& sourceName) override;

  /// Add a patch and return its patchId.
  /// Nomally ra and dec should be filled in, but for moving patches
  /// (e.g. sun) this is not needed.
  /// <br>Optionally it is checked if the patch already exists.
  unsigned int addPatch(const std::string& patchName, int catType,
                        double apparentBrightness, double ra, double dec,
                        bool check) override;

  /// Update the ra/dec and apparent brightness of a patch.
  void updatePatch(unsigned int patchId, double apparentBrightness, double ra,
                   double dec) override;

  /// Add a source to a patch.
  /// Its ra and dec and default parameters will be stored as default
  /// values in the associated ParmDB tables. The names of the parameters
  /// will be preceeded by the source name and a colon.
  /// The map should contain the parameters belonging to the source type.
  /// Missing parameters will default to 0.
  /// <br>Optionally it is checked if the source already exists.
  ///@{
  void addSource(const SourceInfo& sourceInfo, const std::string& patchName,
                 const ParmMap& defaultParameters, double ra, double dec,
                 bool check) override;
  void addSource(const SourceData& source, bool check) override;
  ///@}

  /// Add a source which forms a patch in itself (with the same name).
  /// <br>Optionally it is checked if the patch or source already exists.
  void addSource(const SourceInfo& sourceInfo, const std::string& patchName,
                 int catType, double apparentBrightness,
                 const ParmMap& defaultParameters, double ra, double dec,
                 bool check) override;

  /// Get patch names in order of category and decreasing apparent flux.
  /// category < 0 means all categories.
  /// A brightness < 0 means no test on brightness.
  std::vector<std::string> getPatches(int category, const std::string& pattern,
                                      double minBrightness,
                                      double maxBrightness) override;

  /// Get the info of all patches (name, ra, dec).
  std::vector<PatchInfo> getPatchInfo(int category, const std::string& pattern,
                                      double minBrightness,
                                      double maxBrightness) override;

  /// Get the sources belonging to the given patch.
  std::vector<SourceInfo> getPatchSources(
      const std::string& patchName) override;

  /// Get all data of the sources belonging to the given patch.
  std::vector<SourceData> getPatchSourceData(
      const std::string& patchName) override;

  /// Get the source info of the given source.
  SourceInfo getSource(const std::string& sourceName) override;

  /// Get the info of all sources matching the given (filename like) pattern.
  std::vector<SourceInfo> getSources(const std::string& pattern) override;

  /// Delete the sources records matching the given (filename like) pattern.
  /// This is not possible yet.
  void deleteSources(const std::string& sourceNamePattern) override;

  /// Clear file (i.e. remove everything).
  void clearTables() override;

  /// Get the next source from the table.
  /// An exception is thrown if there are no more sources.
  void getNextSource(SourceData& src) override;

  /// Tell if we are the end of the file.
  bool atEnd() override;

  /// Reset to the beginning of the file.
  void rewind() override;

 private:
  /// Read all patches and sources filling the maps.
  void readAll();

  std::fstream itsFile;
  std::shared_ptr<blob::BlobIBufStream> itsBufIn;
  std::shared_ptr<blob::BlobOBufStream> itsBufOut;
  std::shared_ptr<blob::BlobIStream> itsBlobIn;
  std::shared_ptr<blob::BlobOStream> itsBlobOut;
  bool itsCanWrite;
  int64_t itsEndPos;
  std::map<std::string, PatchInfo> itsPatches;
  std::map<std::string, std::vector<SourceData>>
      itsSources;  ///< sources per patch
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
