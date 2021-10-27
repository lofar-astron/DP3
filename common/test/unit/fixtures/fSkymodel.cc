// Copyright (C) 2021
// ASTRON (Netherlands Institute for Radio Astronomy)

#include "fSkymodel.h"

#include "../../../parmdb/SkymodelToSourceDB.h"

#include <boost/test/unit_test.hpp>

#include <cassert>
#include <fstream>

FixtureSkymodel::FixtureSkymodel(const FixtureSkymodel::Arguments& arguments) {
  assert(!arguments.skymodel_name.empty() &&
         "The skymodel_name should contain a name.");

  {  // Scoped so the file will be closed before it's convert to a source db.
    std::ofstream file(arguments.skymodel_name);
    file << arguments.skymodel_contents;
  }

  if (!arguments.source_db_name.empty()) {
    dp3::parmdb::skymodel_to_source_db::MakeSourceDb(
        arguments.skymodel_name, arguments.source_db_name, std::string(),
        dp3::parmdb::skymodel_to_source_db::ReadFormat("",
                                                       arguments.skymodel_name),
        "", "", false, true, false,
        dp3::parmdb::skymodel_to_source_db::GetSearchInfo("", "", ""));
  }
}

namespace test_source_db {
const std::array<Patch, 3> Expected{
    {{"center", 4.356648, 1.112523, 1.0, 1},
     {"ra_off", 4.356648, 1.112523, 0.5, 1},
     {"radec_off", 4.356648, 1.112523, 0.25, 1}}};

void CheckEqual(const dp3::base::Patch& lhs, const Patch& rhs) {
  BOOST_CHECK_EQUAL(lhs.name(), rhs.name);
  BOOST_CHECK_CLOSE(lhs.direction().ra, rhs.ra, 1e6);
  BOOST_CHECK_CLOSE(lhs.direction().ra, rhs.ra, 1e6);
  BOOST_CHECK_EQUAL(lhs.brightness(), rhs.brightness);
  BOOST_CHECK_EQUAL(lhs.nComponents(), rhs.n_components);
}
}  // namespace test_source_db
