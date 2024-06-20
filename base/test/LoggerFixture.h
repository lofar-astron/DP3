#ifndef DP3_BASE_LOGGER_FIXTURE_H_
#define DP3_BASE_LOGGER_FIXTURE_H_

#include <aocommon/logger.h>

namespace dp3::base::test {

class LoggerFixture {
 public:
  LoggerFixture(aocommon::LogVerbosityLevel test_verbosity =
                    aocommon::LogVerbosityLevel::kQuiet) {
    aocommon::Logger::SetVerbosity(test_verbosity);
  }
  ~LoggerFixture() {
    aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kNormal);
  }
};

}  // namespace dp3::base::test

#endif
