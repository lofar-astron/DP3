// ParameterHandler.h: Handle a LOFAR .parset file
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Handle a LOFAR .parset file
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_LMWCOMMON_PARAMETERHANDLER_H
#define LOFAR_LMWCOMMON_PARAMETERHANDLER_H

#include "ParameterSet.h"

#include "../blob/BlobIStream.h"
#include "../blob/BlobOStream.h"

namespace dp3 {
namespace common {

/// Put ParameterSet into DP3CEP namespace for ease of use.

/// @ingroup LMWCommon
/// @brief Handle a LOFAR .parset file

/// This class handles the processing of a LOFAR .parset file
/// It augments the LOFAR ParameterSet class with functions that can deal
/// with undefined parameters. There is a set of functions that return
/// a default value if undefined and a set of functions that leave the
/// value untouched if undefined.

class ParameterHandler {
 public:
  ParameterHandler(const ParameterSet&);

  /// Get a parameter value.
  /// An exception is thrown if it does not exist.
  /// @{
  std::string getString(const std::string& parm) const;
  double getDouble(const std::string& parm) const;
  unsigned getUint(const std::string& parm) const;
  bool getBool(const std::string& parm) const;
  std::vector<std::string> getStringVector(const std::string& parm) const;
  /// @}

  /// Get a parameter value.
  /// If it does not exist, the default value is used instead.
  /// @{
  std::string getString(const std::string& parm,
                        const std::string& defVal) const;
  double getDouble(const std::string& parm, double defVal) const;
  unsigned getUint(const std::string& parm, unsigned defVal) const;
  bool getBool(const std::string& parm, bool defVal) const;
  std::vector<std::string> getStringVector(
      const std::string& parm, const std::vector<std::string>& defVal) const;
  /// @}

  /// Get a parameter value and fill \a value with it.
  /// If it does not exist, nothing is done.
  /// @{
  void fillString(const std::string& parm, std::string& value) const;
  void fillDouble(const std::string& parm, double& value) const;
  void fillUint(const std::string& parm, unsigned& value) const;
  void fillBool(const std::string& parm, bool& value) const;
  void fillStringVector(const std::string& parm,
                        std::vector<std::string>& value) const;
  /// @}

  /// Convert automatically to a ParameterSet.
  operator const ParameterSet&() const { return itsParms; }

 protected:
  ParameterSet itsParms;
};

/// Write/read a ParameterSet into/from a blob.
/// @{
blob::BlobOStream& operator<<(blob::BlobOStream&, const ParameterSet&);
blob::BlobIStream& operator>>(blob::BlobIStream&, ParameterSet&);
/// @}

inline std::string ParameterHandler::getString(const std::string& parm) const {
  return itsParms.getString(parm);
}
inline double ParameterHandler::getDouble(const std::string& parm) const {
  return itsParms.getDouble(parm);
}
inline unsigned ParameterHandler::getUint(const std::string& parm) const {
  return itsParms.getUint32(parm);
}
inline bool ParameterHandler::getBool(const std::string& parm) const {
  return itsParms.getBool(parm);
}
inline std::vector<std::string> ParameterHandler::getStringVector(
    const std::string& parm) const {
  return itsParms.getStringVector(parm);
}

}  // namespace common
}  // namespace dp3

#endif
