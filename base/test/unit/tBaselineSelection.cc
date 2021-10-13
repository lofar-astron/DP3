// tBaselineSelection.cc: Test program for class BaselneSelection
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>

#include <boost/test/unit_test.hpp>

#include "../../BaselineSelection.h"
#include "../../DPInfo.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StreamUtil.h"

using dp3::base::BaselineSelection;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using std::vector;

BOOST_AUTO_TEST_SUITE(baselineselection)

// Define the info object containing the antenna/baseline info.
DPInfo makeInfo(int nbl) {
  DPInfo info;
  info.init(4, 0, 16, 1, 0.5, 5., string(), string());
  // Fill the baseline stations; use 4 stations.
  // So they are called 00 01 02 03 10 11 12 13 20, etc.
  vector<int> ant1(nbl);
  vector<int> ant2(nbl);
  int st1 = 0;
  int st2 = 0;
  for (int i = 0; i < nbl; ++i) {
    ant1[i] = st1;
    ant2[i] = st2;
    if (++st2 == 4) {
      st2 = 0;
      if (++st1 == 4) {
        st1 = 0;
      }
    }
  }
  vector<string> antNames{"rs01.s01", "rs02.s01", "cs01.s01", "cs01.s02"};
  // Define their positions (more or less WSRT RT0-3).
  // Baseline lengths are:
  //  0-1  144.01
  //  0-2  288.021
  //  0-3  431.914
  //  1-2  144.01
  //  1-3  287.904
  //  2-3  143.896
  vector<casacore::MPosition> antPos(4);
  vector<double> vals(3);
  vals[0] = 3828763;
  vals[1] = 442449;
  vals[2] = 5064923;
  antPos[0] = casacore::MPosition(
      casacore::Quantum<casacore::Vector<double>>(vals, "m"),
      casacore::MPosition::ITRF);
  vals[0] = 3828746;
  vals[1] = 442592;
  vals[2] = 5064924;
  antPos[1] = casacore::MPosition(
      casacore::Quantum<casacore::Vector<double>>(vals, "m"),
      casacore::MPosition::ITRF);
  vals[0] = 3828729;
  vals[1] = 442735;
  vals[2] = 5064925;
  antPos[2] = casacore::MPosition(
      casacore::Quantum<casacore::Vector<double>>(vals, "m"),
      casacore::MPosition::ITRF);
  vals[0] = 3828713;
  vals[1] = 442878;
  vals[2] = 5064926;
  antPos[3] = casacore::MPosition(
      casacore::Quantum<casacore::Vector<double>>(vals, "m"),
      casacore::MPosition::ITRF);
  vector<double> antDiam(4, 70.);
  info.set(antNames, antDiam, antPos, ant1, ant2);
  return info;
}

void test1(const DPInfo& info) {
  ParameterSet ps;
  ps.add("baseline", "[]");
  BaselineSelection selBL(ps, "");
  casacore::Matrix<bool> res = selBL.apply(info);
  BOOST_CHECK(allEQ(res, true));
}

void test2(const DPInfo& info) {
  ParameterSet ps;
  ps.add("baseline", "[]");
  ps.add("corrtype", "AUTO");
  BaselineSelection selBL(ps, "");
  casacore::Matrix<bool> res = selBL.apply(info);
  BOOST_CHECK(allEQ(res.diagonal(), true));
  res.diagonal() = false;
  BOOST_CHECK(allEQ(res, false));
}

void test3(const DPInfo& info) {
  ParameterSet ps;
  ps.add("corrtype", "CROSS");
  BaselineSelection selBL(ps, "");
  casacore::Matrix<bool> res = selBL.apply(info);
  BOOST_CHECK(allEQ(res.diagonal(), false));
  res.diagonal() = true;
  BOOST_CHECK(allEQ(res, true));
}

void test4(const DPInfo& info) {
  casacore::Matrix<double> blength(4, 4, 0.);
  blength(0, 1) = blength(1, 0) = 144.01;
  blength(0, 2) = blength(2, 0) = 288.021;
  blength(0, 3) = blength(3, 0) = 431.914;
  blength(1, 2) = blength(2, 1) = 144.01;
  blength(1, 3) = blength(3, 1) = 287.904;
  blength(2, 3) = blength(3, 2) = 143.896;
  ParameterSet ps;
  ps.add("blmin", "145");
  BaselineSelection selBL(ps, "", true);
  BOOST_CHECK(allEQ(selBL.apply(info), blength <= 145.));
  ps.add("blmax", "288");
  selBL = BaselineSelection(ps, "", true);
  BOOST_CHECK(allEQ(selBL.apply(info), blength <= 145. || blength >= 288.));
  ps.add("corrtype", "cross");
  selBL = BaselineSelection(ps, "", true);
  BOOST_CHECK(allEQ(selBL.apply(info),
                    (blength > 0. && (blength <= 145. || blength >= 288.))));
  ps.add("blrange", "[0,144,288,430]");
  selBL = BaselineSelection(ps, "", true);
  BOOST_CHECK(allEQ(selBL.apply(info),
                    (blength > 0. && (blength <= 145. || blength >= 288.))));
  selBL = BaselineSelection(ps, "", false);
  BOOST_CHECK(allEQ(selBL.apply(info), ((blength > 0. && blength <= 144.) ||
                                        (blength >= 288. && blength <= 430.))));
  ps.replace("corrtype", "");
  selBL = BaselineSelection(ps, "", false);
  BOOST_CHECK(allEQ(selBL.apply(info), ((blength >= 0. && blength <= 144.) ||
                                        (blength >= 288. && blength <= 430.))));
}

void test5(const DPInfo& info) {
  casacore::Matrix<int> bl(4, 4);
  bl(0, 0) = 00;
  bl(1, 1) = 11;
  bl(2, 2) = 22;
  bl(3, 3) = 33;
  bl(0, 1) = bl(1, 0) = 01;
  bl(0, 2) = bl(2, 0) = 02;
  bl(0, 3) = bl(3, 0) = 03;
  bl(1, 2) = bl(2, 1) = 12;
  bl(1, 3) = bl(3, 1) = 13;
  bl(2, 3) = bl(3, 2) = 23;
  ParameterSet ps;
  ps.add("baseline", "[[rs01.s01]]");
  BaselineSelection selBL(ps, "");
  BOOST_CHECK(
      allEQ(selBL.apply(info), bl == 00 || bl == 01 || bl == 02 || bl == 03));
  ps.replace("baseline", "[[rs01.s01,rs*]]");
  selBL = BaselineSelection(ps, "");
  BOOST_CHECK(allEQ(selBL.apply(info), bl == 00 || bl == 01));
  ps.replace("baseline", "[[rs01.s01],[rs02*]]");
  selBL = BaselineSelection(ps, "");
  BOOST_CHECK(allEQ(selBL.apply(info), !(bl == 22 || bl == 23 || bl == 33)));
  ps.replace("baseline", "[rs01.s01,'rs02*']");
  selBL = BaselineSelection(ps, "");
  BOOST_CHECK(allEQ(selBL.apply(info), !(bl == 22 || bl == 23 || bl == 33)));
  ps.replace("corrtype", "cross");
  selBL = BaselineSelection(ps, "");
  BOOST_CHECK(allEQ(selBL.apply(info), !(bl == 00 || bl == 11 || bl == 22 ||
                                         bl == 23 || bl == 33)));
  // Note that the MSSelection syntax is not tested, because it requires
  // a MeasurementSet.
}

BOOST_AUTO_TEST_CASE(test_1) {
  DPInfo info(makeInfo(16));
  test1(info);
}

BOOST_AUTO_TEST_CASE(test_2) {
  DPInfo info(makeInfo(16));
  test2(info);
}

BOOST_AUTO_TEST_CASE(test_3) {
  DPInfo info(makeInfo(16));
  test3(info);
}

BOOST_AUTO_TEST_CASE(test_4) {
  DPInfo info(makeInfo(16));
  test4(info);
}

BOOST_AUTO_TEST_CASE(test_5) {
  DPInfo info(makeInfo(16));
  test5(info);
}

BOOST_AUTO_TEST_SUITE_END()
