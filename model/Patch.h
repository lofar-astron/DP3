// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_PATCH_H_
#define DP3_BASE_PATCH_H_

#include <dp3/base/Direction.h>

#include "../base/ModelComponent.h"

#include <memory>
#include <string>
#include <vector>

namespace dp3::model {

/// \brief A set of sources for which direction dependent effects are assumed to
/// be equal.

class Patch {
 public:
  using iterator = std::vector<std::shared_ptr<base::ModelComponent>>::iterator;
  using const_iterator =
      std::vector<std::shared_ptr<const base::ModelComponent>>::const_iterator;

  template <typename T>
  Patch(const std::string &name, T first, T last)
      : name_(name), components_(first, last) {
    ComputeDirection();
  }

  const std::string &Name() const { return name_; }
  const base::Direction &Direction() const { return direction_; }
  double Brightness() const { return brightness_; }
  /**
   * Index is used by OnePredict to be able to know what beam
   * values belong to what patch. It is the index of this
   * patch in the full list of patches being predicted.
   */
  size_t Index() const { return index_; }

  void SetDirection(const base::Direction &pos) { direction_ = pos; }
  void SetBrightness(double brightness) { brightness_ = brightness; }
  void SetIndex(size_t index) { index_ = index; }

  size_t NComponents() const { return components_.size(); }
  std::shared_ptr<base::ModelComponent> component(size_t i) {
    return components_[i];
  }
  std::shared_ptr<const base::ModelComponent> component(size_t i) const {
    return components_[i];
  }

  iterator begin() { return components_.begin(); }
  iterator end() { return components_.end(); }

  /// Compute the direction as the average of the positions of the components.
  void ComputeDirection();

 private:
  size_t index_ = 0;
  std::string name_;
  base::Direction direction_;
  double brightness_ = 0.0;
  std::vector<std::shared_ptr<base::ModelComponent>> components_;
};

}  // namespace dp3::model

#endif
