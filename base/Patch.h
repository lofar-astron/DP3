// Patch.h: A set of sources for which direction dependent effects are assumed
// to be equal.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_PATCH_H
#define DPPP_PATCH_H

#include "ModelComponent.h"
#include "Position.h"

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
  const Position &position() const;
  double brightness() const;
  void setPosition(const Position &);
  void setBrightness(double);

  size_t nComponents() const;
  ModelComponent::ConstPtr component(size_t i) const;

  const_iterator begin() const;
  const_iterator end() const;

  /// Compute the position as the average of the positions of the components.
  void computePosition();

 private:
  std::string itsName;
  Position itsPosition;
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
  computePosition();
}

inline const std::string &Patch::name() const { return itsName; }

inline const Position &Patch::position() const { return itsPosition; }

inline double Patch::brightness() const { return itsBrightness; }

inline void Patch::setPosition(const Position &pos) { itsPosition = pos; }

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
