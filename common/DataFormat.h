// DataFormat.h: Get the data format (endian type)
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_DATAFORMAT_H
#define LOFAR_COMMON_DATAFORMAT_H

// Never #include <config.h> or #include <lofar_config.h> in a header file!

namespace dp3 {
namespace common {

/// \brief This file defines an enum for the possible machine data formats.

/// Get the data format (endian type).
/// This file defines an enum for the possible machine data formats.
/// Currently only little and big endian is possible with floating point
/// numbers as IEEE and characters in the ASCII representation.
/// It is used in the Blob classes and the DataConvert functions.
///
/// Furthermore it contains a function giving the data format in use on
/// the machine in use.
enum DataFormat { LittleEndian = 0, BigEndian = 1 };

/// Get the endian type on this machine.
inline DataFormat dataFormat()
#if defined(WORDS_BIGENDIAN)
{
  return BigEndian;
}
#else
{
  return LittleEndian;
}
#endif

}  // namespace common
}  // namespace dp3

#endif
