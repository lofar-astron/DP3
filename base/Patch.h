// Patch.h: A set of sources for which direction dependent effects are assumed
// to be equal.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_PATCH_H
#define DPPP_PATCH_H

#include "ModelComponent.h"
#include "Direction.h"

#include <memory>
#include <string>
#include <vector>

namespace dp3 {
namespace base {

/// \brief A set of sources for which direction dependent effects are assumed to
/// be equal.

/// @{

class Patch {
 public:
  typedef std::shared_ptr<Patch> Ptr;
  typedef std::shared_ptr<const Patch> ConstPtr;
  typedef std::vector<ModelComponent::ConstPtr>::const_iterator const_iterator;

  template <typename T>
  Patch(const std::string &name, T first, T last);

  const std::string &name() const;
  const Direction &direction() const;
  double brightness() const;
  void setDirection(const Direction &);
  void setBrightness(double);

  size_t nComponents() const;
  ModelComponent::ConstPtr component(size_t i) const;

  const_iterator begin() const;
  const_iterator end() const;

  /// Compute the direction as the average of the positions of the components.
  void computeDirection();

 private:
  std::string itsName;
  Direction itsDirection;
  double itsBrightness;
  std::vector<ModelComponent::ConstPtr> itsComponents;
};

/// @}

// -------------------------------------------------------------------------- //
// - Implementation: Patch                                                  - //
// -------------------------------------------------------------------------- //

template <typename T>
Patch::Patch(const std::string &name, T first, T last)
    : itsName(name), itsComponents(first, last) {
  computeDirection();
}

inline const std::string &Patch::name() const { return itsName; }

inline const Direction &Patch::direction() const { return itsDirection; }

inline double Patch::brightness() const { return itsBrightness; }

inline void Patch::setDirection(const Direction &pos) { itsDirection = pos; }

inline void Patch::setBrightness(double brightness) {
  itsBrightness = brightness;
}

inline size_t Patch::nComponents() const { return itsComponents.size(); }

inline ModelComponent::ConstPtr Patch::component(size_t i) const {
  return itsComponents[i];
}

inline Patch::const_iterator Patch::begin() const {
  return itsComponents.begin();
}

inline Patch::const_iterator Patch::end() const { return itsComponents.end(); }

}  // namespace base
}  // namespace dp3

#endif
