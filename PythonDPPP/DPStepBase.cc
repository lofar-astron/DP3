// DPStepBase.cc: Python base class for a DPStep in python
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <lofar_config.h>
#include <PythonDPPP/DPStepBase.h>

// The C++ PythonStep must be able to call functions in Python.
// But the Python functions must be able to call the C++ functions
// (e.g., to get data, flags, etc.)
// All communication goes through DPStepBase
// Let a python step be created by a static function creating a
// PythonWorker doing the actual work. Its pointer is kept in a static
// and thereafter used to create a PythonStep.

namespace DP3 {
namespace DPPP {

// Define the static.
DPStepBase* DPStepBase::theirPtr = 0;

}  // namespace DPPP
}  // namespace DP3
