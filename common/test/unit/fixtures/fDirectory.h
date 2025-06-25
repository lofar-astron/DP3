// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_COMMON_TEST_UNIT_FIXTURES_FDIRECTORY_H
#define DP3_COMMON_TEST_UNIT_FIXTURES_FDIRECTORY_H

#include <cstdlib>
#include <string>

#include <boost/filesystem.hpp>

#include "test_config.h"

namespace dp3 {
namespace common {
namespace test {

/// Fixture for Boost unit tests.
/// Creates a temporary directory and set it as working directory.
/// When the test ends the directory is removed.
/// The working directory will be something like <path>/build/abc123/
/// More info:
/// https://www.boost.org/doc/libs/1_60_0/libs/test/doc/html/boost_test/tests_organization/fixtures/case.html
class FixtureDirectory {
 public:
  /// Create the temporary directory and set it as working directory
  FixtureDirectory() {
    boost::filesystem::create_directories(kPath);
    boost::filesystem::current_path(kPath);
  }

  FixtureDirectory(const FixtureDirectory&) = delete;
  FixtureDirectory& operator=(const FixtureDirectory&) = delete;

  /// Remove the temporary diectory
  /// Will always run
  ~FixtureDirectory() {
    boost::filesystem::current_path(kWorkDir);
    boost::filesystem::remove_all(kPath);
  }

  /// Extracts a test resource tgz file into the current directory.
  static void ExtractResource(const std::string& tgz) {
    const std::string command = "tar xfz " DP3_RESOURCE_DIR "/" + tgz;
    const int status = std::system(command.c_str());
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
      throw std::runtime_error("Error while extracting " + tgz);
    }
  }

 private:
  const boost::filesystem::path kPath = boost::filesystem::unique_path();
  const boost::filesystem::path kWorkDir = boost::filesystem::current_path();
};

}  // namespace test
}  // namespace common
}  // namespace dp3

#endif
