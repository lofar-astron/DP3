// ModelComponentVisitor.h: Base class for visitors that visit model component
// hierarchies.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_MODELCOMPONENTVISITOR_H
#define DPPP_MODELCOMPONENTVISITOR_H

namespace DP3 {
namespace DPPP {

class PointSource;
class GaussianSource;

/// \brief Base class for visitors that visit model component hierarchies.

/// @{

class ModelComponentVisitor {
 public:
  virtual ~ModelComponentVisitor();

  virtual void visit(const PointSource&) = 0;
  virtual void visit(const GaussianSource&) = 0;
};

/// @}

}  // namespace DPPP
}  // namespace DP3

#endif
