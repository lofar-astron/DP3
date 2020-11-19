// ModelComponent.h: Base class for model components.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_MODELCOMPONENT_H
#define DPPP_MODELCOMPONENT_H

#include <memory>

namespace DP3 {
namespace DPPP {

class ModelComponentVisitor;
class Position;

/// \brief Base class for model components.

/// @{

class ModelComponent {
 public:
  typedef std::shared_ptr<ModelComponent> Ptr;
  typedef std::shared_ptr<const ModelComponent> ConstPtr;

  virtual ~ModelComponent();
  virtual const Position &position() const = 0;
  virtual void accept(ModelComponentVisitor &) const = 0;
};

/// @}

}  // namespace DPPP
}  // namespace DP3

#endif
