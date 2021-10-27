// Copyright (C) 2021
// ASTRON (Netherlands Institute for Radio Astronomy)

#ifndef DP3_COMMON_TEST_UNIT_FIXTURES_FSKYMODEL_H
#define DP3_COMMON_TEST_UNIT_FIXTURES_FSKYMODEL_H

#include "../../../base/Patch.h"

#include <array>
#include <string>

/// A fixture to write a skymodel information in the current working directory.
///
/// The fixture always writes a skymodel file and optionally a sourceDB. It's
/// recommened to used the @ref FixtureDirectory to get a temporary directory.
///
/// @note The created files will not be cleaned up automatically.
class FixtureSkymodel {
 public:
  /// The constructor arguments.
  ///
  /// A fixture can have only one argument. Using this helper struct works
  /// around the limitation.
  struct Arguments {
#if __cplusplus < 201402L
    // Prior to C++14 the class isn't considered an aggregate. This work-around
    // can be removed once C++14 or newer is required to build DP3.
    explicit Arguments(const std::string& skymodel_name)
        : skymodel_name(skymodel_name) {}

    explicit Arguments(const std::string& skymodel_name,
                       const std::string& source_db_name)
        : skymodel_name(skymodel_name), source_db_name(source_db_name) {}

    explicit Arguments(const std::string& skymodel_name,
                       const std::string& source_db_name,
                       const std::string& skymodel_contents)
        : skymodel_name(skymodel_name),
          source_db_name(source_db_name),
          skymodel_contents(skymodel_contents) {}
#endif

    /// The filename of the skymodel file.
    ///
    /// This argument must contain an non-empty string.
    std::string skymodel_name;

    /// The name of the sourcedb file.
    ///
    /// When this name is not-empty the fixture will create a sourceDB based
    /// on the contents of the skymodel file.
    std::string source_db_name;

    /// The contents of the skymodel file to be written to the disc.
    std::string skymodel_contents =
        R"(FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'
center, POINT, 16:38:28.205000, +63.44.34.314000, 1, , , , ,
ra_off, POINT, 16:58:28.205000, +63.44.34.314000, 0.5, , , , ,
radec_off, POINT, 16:38:28.205000, +65.44.34.314000, 0.25, , , , ,
)";
  };

  explicit FixtureSkymodel(const Arguments& arguments);
};

/// Helper function to be used with @ref FixtureSkymodel.
namespace test_source_db {
struct Patch {
  std::string name;
  double ra;
  double dec;
  double brightness;
  int n_components;
};

/// Contains a set of expected patches.
///
/// This matched the result of using the default value of
/// @ref FixtureSkymodel::Arguments::skymodel_contents.
extern const std::array<Patch, 3> Expected;

/// Helper to validate patch.
///
/// Uses the BOOST_CHECK_x macros to validate whether the arguments are
/// considered equal.
void CheckEqual(const dp3::base::Patch& lhs, const Patch& rhs);

}  // namespace test_source_db
#endif  // DP3_COMMON_TEST_UNIT_FIXTURES_FSKYMODEL_H
