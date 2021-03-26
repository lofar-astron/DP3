// Types.h: Define common types.
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Define common types.
/// @author Maik Nijhuis

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <casacore/casa/version.h>

namespace dp3 {
namespace common {

#if CASACORE_MAJOR_VERSION < 3 || \
    (CASACORE_MAJOR_VERSION == 3 && CASACORE_MINOR_VERSION < 4)
typedef unsigned int rownr_t;
#else
typedef casacore::rownr_t rownr_t;
#endif

}  // namespace common
}  // namespace dp3

#endif
