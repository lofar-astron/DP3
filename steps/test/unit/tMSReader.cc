// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../MSReader.h"
#include "../../../common/ParameterSet.h"

#include <EveryBeam/load.h>
#include <EveryBeam/station.h>

#include <casacore/casa/Arrays/Vector.h>
#include <vector>

using dp3::steps::MSReader;

BOOST_AUTO_TEST_SUITE(msreader)

// Test reading a LOFAR measurement set
BOOST_AUTO_TEST_CASE(read_lofar) {
  std::string msname = "tNDPPP-generic.MS";
  dp3::common::ParameterSet parset;
  MSReader msreader(msname, parset, "");

  std::vector<std::shared_ptr<everybeam::Station>> stations;
  std::vector<std::string> ant_vec = {"CS001HBA0", "CS002HBA0"};
  casacore::Vector<casacore::String> ant_names(ant_vec.size());

  for (size_t i = 0; i < ant_vec.size(); ++i) {
    ant_names[i] = ant_vec[i];
  }

  msreader.fillBeamInfo(stations, ant_names,
                        everybeam::ElementResponseModel::kHamaker);

  for (size_t i = 0; i < stations.size(); ++i) {
    BOOST_CHECK_EQUAL(stations[i]->GetName(), ant_vec[i]);
  }
}

// Reading an OSKAR measurement set
BOOST_AUTO_TEST_CASE(read_oskar) {
  std::string msname = "tOSKAR.in_MS";
  dp3::common::ParameterSet parset;
  MSReader msreader(msname, parset, "");

  std::vector<std::shared_ptr<everybeam::Station>> stations;
  std::vector<std::string> ant_vec = {"s0012", "s0013", "s0015"};
  casacore::Vector<casacore::String> ant_names(ant_vec.size());

  for (size_t i = 0; i < ant_vec.size(); ++i) {
    ant_names[i] = ant_vec[i];
  }

  msreader.fillBeamInfo(stations, ant_names,
                        everybeam::ElementResponseModel::kOSKARSphericalWave);

  for (size_t i = 0; i < stations.size(); ++i) {
    BOOST_CHECK_EQUAL(stations[i]->GetName(), ant_vec[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
