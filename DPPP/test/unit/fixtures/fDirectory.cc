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
struct FixtureDirectory {
  const path ph = unique_path();
  const path workDir = current_path();

  int i;

  FixtureDirectory() : i(0) {
    create_directories(ph);
    current_path(ph);
  }

  /// Remove the temporary diectory
  /// Will run if a test fails or exits due to errors
  ~FixtureDirectory() {
    current_path(workDir);
    remove_all(ph);
  }
};