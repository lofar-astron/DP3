// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <aocommon/py/uniqueptr.h>

#include <dp3/base/DPBuffer.h>

namespace dp3 {
namespace pythondp3 {

/// Wrapper class around DPBuffer, which allows passing a
/// std::unique_ptr<DPBuffer> to Step::process().
using PyDpBuffer = aocommon::py::PyUniquePointer<dp3::base::DPBuffer>;

}  // namespace pythondp3
}  // namespace dp3
