// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)

#include <boost/filesystem.hpp>

using boost::filesystem::create_directories;
using boost::filesystem::current_path;
using boost::filesystem::path;
using boost::filesystem::remove_all;

using namespace boost::filesystem;

/// Fixture for Boost unit tests.
/// Creates a temporary directory and set it as working directory.
/// When the test ends the directory is removed.
/// The working directory will be something like <path>/build/abc123/
/// More info:
/// https://www.boost.org/doc/libs/1_60_0/libs/test/doc/html/boost_test/tests_organization/fixtures/case.html
struct FixtureDirectory {
  const path kPath = unique_path();
  const path kWorkDir = current_path();

  /// Create the temporary directory and set it as working directory
  FixtureDirectory() {
    create_directories(kPath);
    current_path(kPath);
  }

  /// Remove the temporary diectory
  /// Will always run
  ~FixtureDirectory() {
    current_path(kWorkDir);
    remove_all(kPath);
  }
};