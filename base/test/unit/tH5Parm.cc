// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <sstream>
#include <stdio.h>
#include <vector>
#include <casacore/casa/BasicMath/Math.h>

#include <boost/test/unit_test.hpp>

#include "../../H5Parm.h"

using DP3::H5Parm;
using std::vector;

BOOST_AUTO_TEST_SUITE(h5parm)

void checkAxes(H5Parm::SolTab& soltab, size_t ntimes) {
  BOOST_CHECK_EQUAL(soltab.nAxes(), size_t{3});
  BOOST_CHECK(soltab.hasAxis("ant"));
  BOOST_CHECK(soltab.hasAxis("time"));
  BOOST_CHECK(soltab.hasAxis("bla"));
  BOOST_CHECK_EQUAL(soltab.getAxis(0).name, "ant");
  BOOST_CHECK_EQUAL(soltab.getAxis(1).name, "time");
  BOOST_CHECK_EQUAL(soltab.getAxis(2).name, "bla");
  BOOST_CHECK_EQUAL(soltab.getAxis(0).size, size_t{3});
  BOOST_CHECK_EQUAL(soltab.getAxis(1).size, ntimes);
  BOOST_CHECK_EQUAL(soltab.getAxis(2).size, size_t{1});
}

BOOST_AUTO_TEST_CASE(test_gridinterpolate) {
  {
    size_t ntimes = 7;
    {
      // Create a new H5Parm
      H5Parm h5parm("tH5Parm_tmp.h5", true);

      // Check that something is created
      BOOST_CHECK_EQUAL(((H5::H5File&)(h5parm)).getNumObjs(), size_t{1});

      // Check the name of the new solset "sol000"
      BOOST_CHECK_EQUAL(h5parm.getSolSetName(), "sol000");

      // Add some metadata
      vector<string> antNames;
      vector<double> oneAntPos(3);
      vector<vector<double>> antPositions;
      for (unsigned int i = 0; i < 5; ++i) {
        std::stringstream antNameStr;
        antNameStr << "Antenna" << i;
        antNames.push_back(antNameStr.str());
        antPositions.push_back(oneAntPos);
      }
      h5parm.addAntennas(antNames, antPositions);

      vector<H5Parm::AxisInfo> axes;
      axes.push_back(H5Parm::AxisInfo("ant", 3));
      axes.push_back(H5Parm::AxisInfo("time", ntimes));
      axes.push_back(H5Parm::AxisInfo("bla", 1));

      H5Parm::SolTab a = h5parm.createSolTab("mysol", "mytype", axes);

      // Check that the soltab exists
      BOOST_CHECK_EQUAL(h5parm.nSolTabs(), size_t{1});
      BOOST_CHECK(h5parm.hasSolTab("mysol"));

      // Check the axes
      H5Parm::SolTab soltab = h5parm.getSolTab("mysol");
      BOOST_CHECK_EQUAL(soltab.getType(), "mytype");
      checkAxes(soltab, ntimes);

      // Add some data
      vector<double> vals(3 * ntimes);
      vector<double> weights(3 * ntimes);
      for (size_t ant = 0; ant < 3; ++ant) {
        for (size_t time = 0; time < ntimes; ++time) {
          vals[ant * ntimes + time] = 10 * ant + time;
          weights[ant * ntimes + time] = 0.4;
        }
      }

      soltab.setValues(vals, weights, "CREATE with DPPP");

      // Add metadata for stations
      vector<string> someAntNames;
      someAntNames.push_back("Antenna1");
      someAntNames.push_back("Antenna2");
      someAntNames.push_back("Antenna3");
      soltab.setAntennas(someAntNames);

      // Add metadata for times
      vector<double> times;
      for (size_t time = 0; time < ntimes; ++time) {
        times.push_back(57878.5 + 2.0 * time);
      }
      soltab.setTimes(times);

      // Add metadata for freqs;
      vector<double> freqs{130e6, 131e6, 135e6, 137e6};
      soltab.setFreqs(freqs);
    }

    // TODO tests start failing from here
    {
      H5Parm h5parm("tH5Parm_tmp.h5", false, true);
      BOOST_CHECK_EQUAL(h5parm.getSolSetName(), "sol001");
    }

    {
      H5Parm h5parm("tH5Parm_tmp.h5", false, true, "harry");
      BOOST_CHECK_EQUAL(h5parm.getSolSetName(), "harry");
    }

    {
      H5Parm h5parm("tH5Parm_tmp.h5", false, false, "sol000");
      BOOST_CHECK_EQUAL(h5parm.getSolSetName(), "sol000");
      BOOST_CHECK_EQUAL(h5parm.nSolTabs(), size_t{1});
      BOOST_CHECK(h5parm.hasSolTab("mysol"));
      BOOST_CHECK(!h5parm.hasSolTab("nonexistingsol"));

      // Check the axes
      H5Parm::SolTab soltab = h5parm.getSolTab("mysol");
      BOOST_CHECK_EQUAL(soltab.getType(), "mytype");
      checkAxes(soltab, ntimes);

      double starttime = 57878.49999;
      hsize_t starttimeindex = soltab.getTimeIndex(starttime);
      vector<double> val = soltab.getValues("Antenna2", starttimeindex, ntimes,
                                            2, 0, 4, 0, 4, 0);
      BOOST_CHECK(casacore::near(val[0], 10.));
      BOOST_CHECK(casacore::near(val[1], 11.));
      BOOST_CHECK(casacore::near(val[2], 12.));
      BOOST_CHECK(casacore::near(val[3], 13.));
      starttime = 57880.5;
      starttimeindex = soltab.getTimeIndex(starttime);
      BOOST_CHECK_EQUAL(starttimeindex, hsize_t{1});
      vector<double> val2 =
          soltab.getValues("Antenna3", starttimeindex, 2, 2, 0, 4, 0, 4, 0);
      BOOST_CHECK(casacore::near(val2[0], 21.));
      BOOST_CHECK(casacore::near(val2[1], 23.));
      BOOST_CHECK(casacore::near(soltab.getTimeInterval(), 2.));
      vector<string> antennas = soltab.getStringAxis("ant");
      BOOST_CHECK_EQUAL(antennas.size(), size_t{3});
      BOOST_CHECK_EQUAL(antennas[0], "Antenna1");
      BOOST_CHECK_EQUAL(antennas[1], "Antenna2");
      BOOST_CHECK_EQUAL(antennas[2], "Antenna3");
      BOOST_CHECK(casacore::near(soltab.getFreqInterval(0), 1e6));
      BOOST_CHECK(casacore::near(soltab.getFreqInterval(1), 4e6));
      BOOST_CHECK(casacore::near(soltab.getFreqInterval(2), 2e6));

      vector<double> freqs{130e6, 131e6};

      vector<double> times;
      for (size_t time = 0; time < ntimes; ++time) {
        times.push_back(57878.5 + 2.0 * time);
      }

      vector<double> newgridvals = soltab.getValuesOrWeights(
          "val", "Antenna1", times, freqs, 0, 0, true);
      BOOST_CHECK_EQUAL(newgridvals.size(), times.size() * freqs.size());
      size_t idx = 0;
      for (size_t time = 0; time < times.size(); ++time) {
        for (size_t freq = 0; freq < freqs.size(); ++freq) {
          BOOST_CHECK(casacore::near(newgridvals[idx++], double(time)));
        }
      }

      times.clear();
      for (size_t time = 0; time < 3 * ntimes + 2; ++time) {
        times.push_back(57878.5 + 2.0 * time / 3.);
      }
      newgridvals = soltab.getValuesOrWeights("val", "Antenna1", times, freqs,
                                              0, 0, true);
      BOOST_CHECK_EQUAL(newgridvals.size(), times.size() * freqs.size());
      idx = 0;
      for (int time = 0; time < int(times.size()); ++time) {
        for (size_t freq = 0; freq < freqs.size(); ++freq) {
          BOOST_CHECK(
              casacore::near(newgridvals[idx++],
                             min(double((time + 1) / 3), double(ntimes - 1))));
        }
      }
    }
    // Remove the file
    remove("tH5Parm_tmp.h5");
  }
}

BOOST_AUTO_TEST_SUITE_END()
