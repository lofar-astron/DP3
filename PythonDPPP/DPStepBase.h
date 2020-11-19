// DPStepBase.cc: Python base class for a DPStep in python
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_DPSTEPBASE_H
#define DPPP_DPSTEPBASE_H

#include "PythonStep.h"

#include "../DPPP/DPInfo.h"

namespace DP3 {
namespace DPPP {
/// @brief Python base class for a DPStep in python

/// DPStepBase is the interface for callbacks from the DPPP PythonStep
/// class (and its derivations) to C++. Note it is not used for calls
/// from C++ to Python, which are done in PythonStep.cc using the
/// Boost::Python functionality.
/// DPStepBase is bound to the Python interface by PythonDPPP.cc.
/// See class PythonStep and __init.py__ for more information.
class DPStepBase {
 public:
  /// Keep the this pointer in a static. It is used by PythonStep.cc
  /// to call setStep to have a pointer from the python interface
  /// to the C++ code.
  DPStepBase() : itsStep(0) { theirPtr = this; }

  /// Create the link for callbacks from Python to C++.
  /// It is used by the PythonStep constructor.
  void setStep(PythonStep* step) { itsStep = step; }

  /// Get the visibility data into the ValueHolder array.
  /// The array must have the correct type and shape.
  void _getData(const casacore::ValueHolder& vh) { itsStep->getData(vh); }

  /// Get the flags into the ValueHolder array.
  /// The array must have the correct type and shape.
  void _getFlags(const casacore::ValueHolder& vh) { itsStep->getFlags(vh); }

  /// Get the weights into the ValueHolder array.
  /// The array must have the correct type and shape.
  void _getWeights(const casacore::ValueHolder& vh) { itsStep->getWeights(vh); }

  /// Get the UVW coordinates data into the ValueHolder array.
  /// The array must have the correct type and shape.
  void _getUVW(const casacore::ValueHolder& vh) { itsStep->getUVW(vh); }

  /// Get the model data into the ValueHolder array.
  /// The array must have the correct type and shape.
  void _getModelData(const casacore::ValueHolder& vh) {
    itsStep->getModelData(vh);
  }

  /// Call the process function in the next DPPP step.
  /// The record must contain the changed data fields which will be
  /// stored in the output buffer before calling the next step.
  bool _processNext(const casacore::Record& rec) {
    return itsStep->processNext(rec);
  }

  /// Keep the pointer to this object.
  /// It is filled by the constructor and used shortly thereafter
  /// in the PythonStep constructor.
  static DPStepBase* theirPtr;

 private:
  PythonStep* itsStep;
};

}  // namespace DPPP
}  // namespace DP3

#endif
