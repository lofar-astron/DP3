// ParmSet.h: Set of parameters to be used
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Set of parameters to be used
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMSET_H
#define LOFAR_PARMDB_PARMSET_H

#include <casacore/casa/Arrays/Matrix.h>

#include <map>

namespace dp3 {
namespace parmdb {

class ParmDB;
class ParmValueSet;
class Box;

/// Define the type of a parmId.
typedef unsigned int ParmId;

/// @ingroup ParmDB
/// @{

/// @brief Set of parameters to be used
class ParmSet {
 public:
  /// Create an empty ParmSet.
  /// It can be filled using function addParm.
  ParmSet();

  /// Add a parm for the given ParmDB. If not existing in the parmDB,
  /// it will be added to the ParmDB when its values are written.
  /// <br>It returns a unique parmId.
  /// <br>If the parm has been added before to the ParmSet, it returns
  /// the already known parmid.
  ParmId addParm(ParmDB&, const std::string& name);

  /// Does the parm already exist in the ParmDB?
  bool isInParmDB(ParmId parmid) const {
    return itsParms[parmid].getNameId() >= 0;
  }

  /// Add zero or more existing parms.
  /// A vector of unique parmIds is returned.
  ///@{
  ///    vector<ParmId> addParms (ParmDB&, const vector<string>& names);
  ///    vector<ParmId> addParms (ParmDB&, const string& namePattern);
  ///@}

  /// Get the nr of parameters.
  size_t size() const { return itsParms.size(); }

  /// Find the parmid of a previously added parm.
  /// An exception is thrown if not existing in the ParmSet.
  ParmId find(const std::string& name) const;

  /// Get the ParmDBs used in the ParmSet.
  const std::vector<ParmDB*> getDBs() const { return itsDBs; }

  /// Get the values for the given work domain for all parameters that are
  /// not part of the given vector.
  /// It means that the vector gets extended for all parameters at the end
  /// of itsParms.
  void getValues(std::vector<ParmValueSet>&, const Box& workDomain) const;

  /// Write the parm values for the given parmid.
  /// For parms with a new name, ParmKey::itsNameId will get filled in.
  /// For new values, the rowId in the ParmValueSet will be filled in.
  void write(unsigned int parmId, ParmValueSet&);

  /// Clear the ParmSet.
  void clear();

 private:
  /// Rescale a polynomial if needed.
  void rescale(ParmValueSet& pset, const Box& newDomain) const;

  /// The nested class ParmKey holds the basic info for a parameter.
  class ParmKey {
   public:
    /// Create a parameter key.
    ParmKey(ParmDB* parmdb, const std::string& name, int nameId, ParmId parmId)
        : itsParmDB(parmdb),
          itsName(name),
          itsNameId(nameId),
          itsParmId(parmId) {}

    /// Get the name.
    const std::string& getName() const { return itsName; }

    /// Get the ParmDB used for the parameter.
    ///@{
    ParmDB* getParmDBPtr() { return itsParmDB; }
    const ParmDB* getParmDBPtr() const { return itsParmDB; }
    ///@}

    /// Get the name id.
    /// >= 0 means name is stored in name table;
    /// <0   means not stored.
    ///@{
    int& getNameId() { return itsNameId; }
    int getNameId() const { return itsNameId; }
    ///@}

    /// Get he unique parm id.
    ParmId getParmId() const { return itsParmId; }

   private:
    ParmDB* itsParmDB;
    std::string itsName;
    int itsNameId;
    ParmId itsParmId;
  };

  std::vector<ParmDB*> itsDBs;
  std::vector<ParmKey> itsParms;
  std::map<std::string, int> itsNames;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
