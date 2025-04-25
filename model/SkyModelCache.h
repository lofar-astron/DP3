#ifndef SKY_MODEL_CACHE_H_
#define SKY_MODEL_CACHE_H_

#include <map>
#include <stdexcept>
#include <string>

#include "SourceDBUtil.h"

namespace dp3::model {

class SkyModelCache {
 public:
  static SkyModelCache& GetInstance() {
    static SkyModelCache instance;
    return instance;
  }

  SourceDBWrapper GetSkyModel(const std::string& filename) {
    if (!HasSkymodelExtension(filename)) {
      throw std::runtime_error(
          "Sky models in the 'source db' format are no longer supported by "
          "DP3. Models need to be provided in the textual format and do not "
          "require conversion using the 'makesourcedb' tool. If you do not "
          "have the original text model file, the 'showsourcedb' tool can be "
          "used to convert a sourcedb back to text format.");
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
  std::map<std::string, SourceDBWrapper> cache_;
};

}  // namespace dp3::model

#endif
