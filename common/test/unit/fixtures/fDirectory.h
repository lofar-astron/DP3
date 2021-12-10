// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_COMMON_TEST_UNIT_FIXTURES_FDIRECTORY_H
#define DP3_COMMON_TEST_UNIT_FIXTURES_FDIRECTORY_H

#include <boost/filesystem.hpp>

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

 private:
  const boost::filesystem::path kPath = boost::filesystem::unique_path();
  const boost::filesystem::path kWorkDir = boost::filesystem::current_path();
};

#endif
