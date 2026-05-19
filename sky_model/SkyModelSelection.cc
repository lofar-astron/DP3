#include "ReadSkyModel.h"
#include "SkyModelFunctions.h"
#include "SkyModelSelection.h"

namespace dp3::sky_model {

std::vector<std::shared_ptr<Patch>> SkyModelSelection::MakePatchList() {
  return MakePatches(sky_model_, patch_names_);
}

bool SkyModelSelection::CheckPolarized() {
  return sky_model::CheckPolarized(sky_model_, patch_names_);
}

bool SkyModelSelection::CheckAnyOrientationIsAbsolute() {
  return sky_model::CheckAnyOrientationIsAbsolute(sky_model_, patch_names_);
}

void SkyModelSelection::SelectAllPatches() {
  patch_names_ = sky_model_.GetPatchNames();
}

void SkyModelSelection::SelectMatchingPatches(
    const std::vector<std::string>& filter) {
  if (std::find(filter.cbegin(), filter.cend(), "") != filter.end()) {
    throw std::runtime_error("Empty source pattern not allowed");
  }
  patch_names_ = sky_model::MakePatchList(sky_model_, filter);
}

void SkyModelSelection::SelectPatchList(
    const std::vector<std::string>& patch_names) {
  if (std::find(patch_names.cbegin(), patch_names.cend(), "") !=
      patch_names.end()) {
    throw std::runtime_error("Empty patch name not allowed");
  }
  patch_names_ = patch_names;
}

}  // namespace dp3::sky_model
