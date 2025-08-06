// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Counter.h"

#include <boost/test/unit_test.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "../../../common/ParameterSet.h"
#include "mock/MockStep.h"

#include <iostream>
#include <filesystem>
#include <fstream>

using dp3::steps::Counter;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(counter, *boost::unit_test::tolerance(0.001) *
                                   boost::unit_test::tolerance(0.001f))

BOOST_AUTO_TEST_CASE(fields) {
  dp3::common::ParameterSet parset;

  const Counter counter(parset, "");
  BOOST_TEST(counter.getRequiredFields() == Step::kFlagsField);
  BOOST_TEST(counter.getProvidedFields() == dp3::common::Fields());
}

BOOST_AUTO_TEST_CASE(save_ratios_to_json) {
  unsigned int n_corr = 4;
  unsigned int n_chan = 4;
  unsigned int n_baselines = 3;

  const std::vector<std::string> ant_names{"ant0", "ant1", "ant2"};

  const std::vector<casacore::MPosition> ant_pos{
      casacore::MVPosition{0, 0, 0}, casacore::MVPosition{300, 0, 0},
      casacore::MVPosition{200, 0, 0}};
  const std::vector<double> ant_diam(3, 1.0);
  const std::vector<int> ant1{0, 1, 2};
  const std::vector<int> ant2{1, 2, 0};

  dp3::base::DPInfo info(n_corr, n_chan);
  info.setTimes(0.0, 0.0, 1.0);
  info.setAntennas(ant_names, ant_diam, ant_pos, ant1, ant2);

  dp3::common::ParameterSet parset;
  std::string test_filename = "flag_counter_test.json";
  parset.add("savetojson", "true");
  parset.add("jsonfilename", test_filename);

  const std::vector<std::size_t> channel_counts(n_chan, 1);
  auto buffer = std::make_unique<dp3::base::DPBuffer>();
  buffer->GetFlags().resize({n_baselines, n_chan, n_corr});
  buffer->GetFlags().fill(false);

  std::vector<std::vector<std::vector<bool>>> flag_values;
  // Define flags for baseline 1 (antennas: 0, 1)
  flag_values.push_back({{true, true, true, true},
                         {true, true, true, true},
                         {true, true, true, true},
                         {true, true, true, true}});
  // Define flags for baseline 2 (antennas: 1, 2)
  flag_values.push_back({{false, false, false, false},
                         {false, false, false, false},
                         {false, false, false, false},
                         {false, false, false, false}});
  // Define flags for baseline 3 (antennas: 2, 0)
  flag_values.push_back({{true, true, true, true},
                         {true, true, true, true},
                         {false, false, false, false},
                         {false, false, false, false}});

  // Expected ratios per antenna (based on flags above)
  std::vector<std::string> expected_ratio_per_antenna{"0.75", "0.5", "0.25"};

  bool current_val = true;
  for (unsigned int bl = 0; bl < n_baselines; ++bl) {
    current_val = !current_val;
    for (unsigned int chan = 0; chan < n_chan; ++chan) {
      for (unsigned int corr = 0; corr < n_corr; ++corr) {
        buffer->GetFlags()(bl, chan, corr) = flag_values[bl][chan][corr];
      }
    }
  }

  auto mock_step = std::make_shared<dp3::steps::MockStep>();

  Counter counter(parset, "");
  counter.updateInfo(info);

  counter.setNextStep(mock_step);

  counter.process(std::move(buffer));

  std::ostringstream os;
  counter.showCounts(os);

  // assert that file exists
  BOOST_REQUIRE(std::filesystem::exists(test_filename));

  // assert that values are correct
  boost::property_tree::ptree root;
  // Load the json file in this ptree
  boost::property_tree::read_json(test_filename, root);
  boost::property_tree::ptree& content =
      root.get_child("flagged_fraction_dict");

  int i = 0;
  for (auto it = content.begin(); it != content.end(); ++it) {
    BOOST_CHECK_EQUAL(it->first, ant_names[i]);
    BOOST_CHECK_EQUAL(it->second.get_value<std::string>(it->first),
                      expected_ratio_per_antenna[i]);
    i++;
  }

  // delete file
  std::filesystem::remove(test_filename);
}

BOOST_AUTO_TEST_SUITE_END()
