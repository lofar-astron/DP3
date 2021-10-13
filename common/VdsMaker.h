// VdsMaker.h: Class to create the description of an MS
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// Class to create the description of an MS
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_MS_VDSMAKER_H
#define LOFAR_MS_VDSMAKER_H

#include "ClusterDesc.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/casa/Arrays/Vector.h>

namespace dp3 {
namespace common {

/// @ingroup MS
/// @brief Class to create the description of an MS
/// @{

/// Class to create the description of a MeasurementSet (part).
/// The output is a vds file written in parset format
/// (vds = visibility data set).
/// All vds files of an entire MeasurementSet can be combined into
/// a gvds file which is used by the pipeline software to determine
/// on which nodes processses can be run the process the various
/// parts of an MS.
//
/// The vds file can contain extra info. For example, the program
/// extra information is added for the RFI detection software of
/// Peter Fridman.

class VdsMaker {
 public:
  /// Create the description for the given MS and put it in a file
  /// with the given name on the given host.
  /// If the host name is empty, gethostname() will be used.
  /// The given ClusterDesc file will be used to find the file system
  /// on which the MS part os located. If the clusterDescName is empty,
  /// the file system will be set to unknown.
  /// It can be specified if the vectors holding the start and end time
  /// of each time stamp should be made part of the VDS file.
  static void create(const string& msName, const string& outName,
                     const string& clusterDescName,
                     const string& hostName = string(), bool fillTimes = true);

  /// Combine the given VDS file into a global VDS file.
  static void combine(const string& gdsName,
                      const std::vector<string>& vdsNames);

 private:
  /// Get the frequency info for each spectral window in the MS.
  /// The vectors get the start and end frequency of each channel.
  static void getFreqInfo(casacore::MS& ms, std::vector<int>& nrchan,
                          std::vector<casacore::Vector<double>>& startFreq,
                          std::vector<casacore::Vector<double>>& endFreq);

  /// Get the directions of the fields.
  static void getFields(casacore::MS& ms, std::vector<double>& ra,
                        std::vector<double>& dec, std::vector<string>& refType);

  /// Get the names of the antennae (stations).
  static void getAntNames(casacore::MS& ms, std::vector<string>& antNames);

  /// Get the names of the correlations (polarisations).
  static void getCorrInfo(casacore::MS& ms, std::vector<string>& corrTypes);

  /// Find out which file contains the DATA column.
  /// Determine if the DATA are stored in a TSM file of itself.
  /// Determine the cube and tile shape.
  static void getDataFileInfo(casacore::MS& ms, string& name, bool& regular,
                              std::vector<int>& tileShape,
                              std::vector<int>& cubeShape);

  /// Find the file system on which the given file is located.
  /// If the host name is empty, gethostname() will be used.
  static string findFileSys(const string& fileName,
                            const common::ClusterDesc& cdesc,
                            const string& hostName);
};

/// @}

}  // namespace common
}  // namespace dp3

#endif
