// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../Simulator.h"
#include "../../Stokes.h"
#include "../../Direction.h"
#include "../../PointSource.h"

#include <sstream>
#include <complex>
#include <vector>

namespace dp3 {
namespace base {
namespace test {
BOOST_AUTO_TEST_SUITE(tsimulator)

const Direction kReference(0.5, 0.1);
// Offset 1/2 deg in radians in RA and DEC
const Direction kOffsetSource(kReference.ra + 0.02, kReference.dec + 0.02);
const size_t kNStations = 4;
const size_t kNChan = 2;

Simulator MakeSimulator(bool correct_freq_smearing, bool stokes_i_only,
                        casacore::Cube<std::complex<double>>& buffer) {
  std::vector<Baseline> baselines;
  for (size_t st1 = 0; st1 < kNStations - 1; ++st1) {
    for (size_t st2 = st1 + 1; st2 < kNStations; ++st2) {
      baselines.emplace_back(Baseline(st1, st2));
    }
  }

  std::vector<double> chan_freqs;
  for (size_t chan = 0; chan < kNChan; ++chan) {
    chan_freqs.push_back(130.0e6 + chan * 1.0e6);
  }
  std::vector<double> chan_widths(kNChan, 1.0e6);

  casacore::Matrix<double> uvw(3, kNStations);

  for (size_t st = 0; st < kNStations; ++st) {
    uvw(0, st) = st * 5000;
    uvw(1, st) = st * 1000;
    uvw(2, st) = 0;
  }

  return Simulator(kReference, kNStations, baselines, chan_freqs, chan_widths,
                   uvw, buffer, correct_freq_smearing, stokes_i_only);
}

BOOST_AUTO_TEST_CASE(test_pointsource_onlyI) {
  Stokes unit;
  unit.I = 1.0;
  unit.Q = 0.0;
  unit.U = 0.0;
  unit.V = 0.0;

  auto pointsource = std::make_shared<PointSource>(kOffsetSource, unit);

  const size_t nbaselines = kNStations * (kNStations - 1) / 2;
  casacore::Cube<std::complex<double>> buffer(1, kNChan, nbaselines);

  Simulator sim = MakeSimulator(false, true, buffer);

  sim.simulate(pointsource);

  BOOST_CHECK_CLOSE(abs(buffer(0, 0, 0)), 1.0, 1.0e-3);  // Channel 0
  BOOST_CHECK_CLOSE(abs(buffer(1, 0, 0)), 1.0, 1.0e-3);  // Channel 1
}

BOOST_AUTO_TEST_CASE(test_pointsource_fullstokes) {
  Stokes unit;
  unit.I = 1.0;
  unit.Q = 0.0;
  unit.U = 0.0;
  unit.V = 0.0;

  auto phasecenter_point = std::make_shared<PointSource>(kOffsetSource, unit);

  const size_t nbaselines = kNStations * (kNStations - 1) / 2;
  casacore::Cube<std::complex<double>> buffer(4, kNChan, nbaselines);

  Simulator sim = MakeSimulator(false, false, buffer);

  sim.simulate(phasecenter_point);

  BOOST_CHECK_CLOSE(abs(buffer(0, 0, 0)), 1.0, 1.0e-3);  // Channel 0, XX
  BOOST_CHECK_CLOSE(abs(buffer(1, 0, 0)), 0.0, 1.0e-3);  // Channel 0, XY
  BOOST_CHECK_CLOSE(abs(buffer(2, 0, 0)), 0.0, 1.0e-3);  // Channel 0, YX
  BOOST_CHECK_CLOSE(abs(buffer(3, 0, 0)), 1.0, 1.0e-3);  // Channel 0, YY
  BOOST_CHECK_CLOSE(abs(buffer(0, 1, 0)), 1.0, 1.0e-3);  // Channel 1, XX
  BOOST_CHECK_CLOSE(abs(buffer(1, 1, 0)), 0.0, 1.0e-3);  // Channel 1, XY
  BOOST_CHECK_CLOSE(abs(buffer(2, 1, 0)), 0.0, 1.0e-3);  // Channel 1, YX
  BOOST_CHECK_CLOSE(abs(buffer(3, 1, 0)), 1.0, 1.0e-3);  // Channel 1, YY
}

BOOST_AUTO_TEST_CASE(test_pointsource_onlyI_freqsmear) {
  Stokes unit;
  unit.I = 1.0;
  unit.Q = 0.0;
  unit.U = 0.0;
  unit.V = 0.0;

  auto phasecenter_offset = std::make_shared<PointSource>(kOffsetSource, unit);

  const size_t nbaselines = kNStations * (kNStations - 1) / 2;

  casacore::Cube<std::complex<double>> buffer(1, kNChan, nbaselines);

  Simulator sim = MakeSimulator(true, true, buffer);

  sim.simulate(phasecenter_offset);

  BOOST_CHECK_CLOSE(std::abs(buffer(0, 0, kNStations - 1)), 0.759154, 1.0e-3);
  BOOST_CHECK_CLOSE(std::abs(buffer(1, 0, kNStations - 1)), 0.759154, 1.0e-3);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace test
}  // namespace base
}  // namespace dp3
