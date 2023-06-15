// ComponentInfo.h: Class for getting information from model component
// hierarchies.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_COMPONENTINFO_H
#define DP3_BASE_COMPONENTINFO_H

#include "ModelComponent.h"
#include "ModelComponentVisitor.h"
#include <dp3/base/Direction.h>
#include "Stokes.h"
#include <vector>

namespace dp3 {
namespace base {

/// \brief Class for visitors that visit model component to extract information.

/// @{

class ComponentInfo : public ModelComponentVisitor {
 public:
  enum SourceType { kPoint, kGaussian };
  ComponentInfo();

  void inspect(const std::shared_ptr<const ModelComponent>& component);

  double ra_, dec_;
  double sI_, sQ_, sU_, sV_;
  SourceType source_type_;
  std::vector<double> spectrum_ = {0.0, 0.0, 0.0};
  double f0_;

  // Only used in Gaussians
  double g_pa_{0.0}, g_major_{1.0}, g_minor_{1.0};

 private:
  void visit(const PointSource&) override;
  void visit(const GaussianSource&) override;

  void update(const PointSource& component);
};

/// @}

}  // namespace base
}  // namespace dp3

#endif
