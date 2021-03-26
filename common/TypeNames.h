// TypeNames.h: Return a string giving the type name to be stored in blobs
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_TYPENAMES_H
#define LOFAR_COMMON_TYPENAMES_H

#include <complex>
#include <string>

namespace dp3 {
namespace common {

/// \ingroup TypeNames
/// \brief Return a string giving the type name to be stored in blobs.

///
/// These global functions return the name of the basic types.
/// They are meant to get the full id of a templated class when such an
/// object is stored in a blob.
/// As much as possible std::complex and builtin complex types get the same
/// name, so they can be read back from a blob in both ways.
/// @{

const std::string& typeName(const void*);
const std::string& typeName(const bool*);
const std::string& typeName(const char*);
const std::string& typeName(const int8_t*);
const std::string& typeName(const uint8_t*);
const std::string& typeName(const int16_t*);
const std::string& typeName(const uint16_t*);
const std::string& typeName(const int32_t*);
const std::string& typeName(const uint32_t*);
const std::string& typeName(const int64_t*);
const std::string& typeName(const uint64_t*);
const std::string& typeName(const float*);
const std::string& typeName(const double*);
const std::string& typeName(const std::complex<float>*);
const std::string& typeName(const std::complex<double>*);
#ifdef LOFAR_BUILTIN_COMPLEXINT
const std::string& typeName(const std::complex<int16>*);
const std::string& typeName(const std::complex<uint16>*);
#endif
template <typename T>
const std::string& typeName(T const* const*);

/// @}

}  // namespace common
}  // namespace dp3

// Include templated implementations.
#include "TypeNames.tcc"

#endif
