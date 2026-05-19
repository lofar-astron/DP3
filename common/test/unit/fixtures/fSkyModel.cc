// Copyright (C) 2021
// ASTRON (Netherlands Institute for Radio Astronomy)

#include "fSkyModel.h"

#include "sky_model/ReadSkyModel.h"

#include <boost/test/unit_test.hpp>

#include <cassert>
#include <fstream>

FixtureSkyModel::FixtureSkyModel(const FixtureSkyModel::Arguments& arguments) {
  assert(!arguments.sky_model_name.empty() &&
         "The sky_model_name should contain a name.");

  {  // Scoped so the file will be closed before it's convert to a source db.
    std::ofstream file(arguments.sky_model_name);
    file << arguments.sky_model_contents;
  }
}

namespace test_source_db {
const std::array<Patch, 3> Expected{
    {{"center", 4.356648, 1.112523, 1.0, 1},
     {"ra_off", 4.356648, 1.112523, 0.5, 1},
     {"radec_off", 4.356648, 1.112523, 0.25, 1}}};

void CheckEqual(const dp3::sky_model::Patch& lhs,
                const test_source_db::Patch& rhs) {
  BOOST_CHECK_EQUAL(lhs.Name(), rhs.name);
  BOOST_CHECK_CLOSE(lhs.Direction().ra, rhs.ra, 1e6);
  BOOST_CHECK_CLOSE(lhs.Direction().ra, rhs.ra, 1e6);
  BOOST_CHECK_EQUAL(lhs.Brightness(), rhs.brightness);
  BOOST_CHECK_EQUAL(lhs.NComponents(), rhs.n_components);
}
}  // namespace test_source_db
