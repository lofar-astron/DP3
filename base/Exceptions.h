// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_EXCEPTIONS_H
#define DP3_EXCEPTIONS_H

#include <stdexcept>

// Note: This Exception should be phased out in the future.
// Please use std::runtime_error in new code instead.
namespace dp3 {
using Exception = std::runtime_error;
}

#endif
