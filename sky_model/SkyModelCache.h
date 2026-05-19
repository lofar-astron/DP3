#ifndef SKY_MODEL_CACHE_H_
#define SKY_MODEL_CACHE_H_

#include <map>
#include <string>

#include "ReadSkyModel.h"
#include "SkyModel.h"
#include "SkyModelSelection.h"

namespace dp3::sky_model {

class SkyModelCache {
 public:
  static SkyModelCache& GetInstance() {
    static SkyModelCache instance;
    return instance;
  }

  SkyModel GetSkyModel(const std::string& filename) {
    auto iterator = cache_.find(filename);
    if (iterator == cache_.end()) {
      return cache_.insert_or_assign(filename, ReadSkyModel(filename))
          .first->second;
    } else {
      return iterator->second;
    }
  }

  SkyModelSelection GetMatchedSkyModel(
      const std::string& filename, const std::vector<std::string>& patterns) {
    SkyModelSelection selection(GetSkyModel(filename));
    selection.SelectMatchingPatches(patterns);
    return selection;
  }

  SkyModelSelection GetSkyModelPatches(
      const std::string& filename,
      const std::vector<std::string>& patch_names) {
    SkyModelSelection selection(GetSkyModel(filename));
    selection.SelectPatchList(patch_names);
    return selection;
  }

  void Clear() { cache_.clear(); }

 private:
  std::map<std::string, SkyModel> cache_;
};

}  // namespace dp3::sky_model

#endif
