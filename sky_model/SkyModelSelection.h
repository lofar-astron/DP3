#ifndef DP3_SKY_MODEL_SKY_MODEL_SELECTION_H_
#define DP3_SKY_MODEL_SKY_MODEL_SELECTION_H_

#include <optional>
#include <string>
#include <vector>

#include "Patch.h"
#include "SkyModel.h"

namespace dp3::sky_model {

class SkyModelSelection {
 public:
  explicit SkyModelSelection(sky_model::SkyModel sky_model)
      : sky_model_(sky_model) {}

  void SelectAllPatches();

  void SelectPatchList(const std::vector<std::string> &patch_names);

  void SelectMatchingPatches(const std::vector<std::string> &filter);

  std::vector<std::shared_ptr<Patch>> MakePatchList();

  bool CheckPolarized();

  bool CheckAnyOrientationIsAbsolute();

 private:
  std::vector<std::string> patch_names_;
  sky_model::SkyModel sky_model_;
};

}  // namespace dp3::sky_model

#endif
