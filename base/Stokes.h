// Stokes.h: Complex Stokes vector.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_STOKES_H
#define DPPP_STOKES_H

namespace dp3 {
namespace base {

/// \brief Complex Stokes vector.

/// @{

class Stokes {
 public:
  Stokes();

  double I, Q, U, V;
};

/// @}

}  // namespace base
}  // namespace dp3

#endif
