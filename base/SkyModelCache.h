#ifndef SKY_MODEL_CACHE_H_
#define SKY_MODEL_CACHE_H_

#include <map>
#include <string>

#include "SourceDBUtil.h"

namespace dp3::base {

class SkyModelCache {
 public:
  static SkyModelCache& GetInstance() {
    static SkyModelCache instance;
    return instance;
  }

  SourceDBWrapper GetSkyModel(const std::string& filename) {
    if (!HasSkymodelExtension(filename)) {
      // SourceDBs can not be value-copied (they use reference copy semantics).
      // Therefore, don't use the cache for them
      return SourceDBWrapper(filename);
    } else {
      auto iterator = cache_.find(filename);
      if (iterator == cache_.end()) {
        return cache_.insert_or_assign(filename, SourceDBWrapper(filename))
            .first->second;
      } else {
        return iterator->second;
      }
    }
  }

  void Clear() { cache_.clear(); }

 private:
  std::map<std::string, base::SourceDBWrapper> cache_;
};

}  // namespace dp3::base

#endif
