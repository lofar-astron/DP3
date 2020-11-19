// MWError.h: Basic exception for master/worker related errors
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Basic exception for master/worker related errors.
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_LMWCOMMON_MWERROR_H
#define LOFAR_LMWCOMMON_MWERROR_H

namespace DP3 {
namespace CEP {

/// @ingroup LMWCommon
/// @brief Basic exception for master/worker related errors.

/// This class defines the basic MW exception.
/// Only this basic exception is defined so far. In the future, some more
/// fine-grained exceptions might be derived from it.
typedef std::runtime_error MWError;

}  // namespace CEP
}  // namespace DP3

#endif
