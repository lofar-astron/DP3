// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_GAIN_SOLVERS_KERNELS_COMMON_H_
#define DP3_DDECAL_GAIN_SOLVERS_KERNELS_COMMON_H_

#include <cudawrappers/cu.hpp>

/// This helper function is needed because the Launch*Kernel functions receive
/// cu::DeviceMemory references, while the GPU kernels require the actual
/// pointer type instead.
template <typename T>
T* Cast(cu::DeviceMemory& m) {
  return reinterpret_cast<T*>(static_cast<CUdeviceptr>(m));
}

#endif  // DP3_DDECAL_GAIN_SOLVERS_KERNELS_COMMON_H_