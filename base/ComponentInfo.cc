// ComponentInfo.h: Class for getting information from model component
// hierarchies.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <iostream>
#include <algorithm>

#include "ComponentInfo.h"
#include "GaussianSource.h"
#include "PointSource.h"

namespace dp3 {
namespace base {

ComponentInfo::ComponentInfo() {}

void ComponentInfo::inspect(
    const std::shared_ptr<const ModelComponent>& component) {
  component->accept(*this);
}

void ComponentInfo::update(const PointSource& component) {
  const Direction& direction = component.direction();
  ra_ = direction.ra;
  dec_ = direction.dec;
  const Stokes& stokes = component.stokes();
  sI_ = stokes.I;
  sQ_ = stokes.Q;
  sU_ = stokes.U;
  sV_ = stokes.V;
  const std::vector<double>& spec = component.spectrum();
  // Note: we can only handle spectra with 3 terms
  if (spec.size() > 3) {
    throw std::runtime_error(
        "Cannot handle spectra with more than three components");
  }
  const size_t spec_size = (spec.size() > 3 ? 3 : spec.size());
  for (size_t index = 0; index < spec_size; index++) {
    spectrum_[index] = spec[index];
  }
  f0_ = component.referenceFreq();
}

void ComponentInfo::visit(const PointSource& component) {
  update(component);
  source_type_ = kPoint;
}

void ComponentInfo::visit(const GaussianSource& component) {
  update(component);
  source_type_ = kGaussian;
  g_pa_ = component.getPositionAngle();
  g_major_ = component.getMajorAxis();
  g_minor_ = component.getMinorAxis();
}
}  // namespace base
}  // namespace dp3
