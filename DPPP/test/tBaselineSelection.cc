//# tBaselineSelection.cc: Test program for class BaselneSelection
//# Copyright (C) 2012
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/BaselineSelection.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <Common/LofarLogger.h>
#include <Common/StreamUtil.h>
#include <casa/Arrays/ArrayLogical.h>
#include <iostream>

using namespace LOFAR;
using namespace LOFAR::DPPP;
using namespace casa;
using namespace std;


// Define the info object containing the antenna/baseline info.
DPInfo makeInfo (int nbl)
{
  DPInfo info;
  info.init (4, 16, 1, 0.5, 5., string());
  // Fill the baseline stations; use 4 stations.
  // So they are called 00 01 02 03 10 11 12 13 20, etc.
  Vector<Int> ant1(nbl);
  Vector<Int> ant2(nbl);
  int st1 = 0;
  int st2 = 0;
  for (int i=0; i<nbl; ++i) {
    ant1[i] = st1;
    ant2[i] = st2;
    if (++st2 == 4) {
      st2 = 0;
      if (++st1 == 4) {
        st1 = 0;
      }
    }
  }
  Vector<String> antNames(4);
  antNames[0] = "rs01.s01";
  antNames[1] = "rs02.s01";
  antNames[2] = "cs01.s01";
  antNames[3] = "cs01.s02";
  // Define their positions (more or less WSRT RT0-3).
  // Baseline lengths are:
  //  0-1  144.01
  //  0-2  288.021
  //  0-3  431.914
  //  1-2  144.01
  //  1-3  287.904
  //  2-3  143.896
  vector<MPosition> antPos(4);
  Vector<double> vals(3);
  vals[0] = 3828763; vals[1] = 442449; vals[2] = 5064923;
  antPos[0] = MPosition(Quantum<Vector<double> >(vals,"m"),
                        MPosition::ITRF);
  vals[0] = 3828746; vals[1] = 442592; vals[2] = 5064924;
  antPos[1] = MPosition(Quantum<Vector<double> >(vals,"m"),
                        MPosition::ITRF);
  vals[0] = 3828729; vals[1] = 442735; vals[2] = 5064925;
  antPos[2] = MPosition(Quantum<Vector<double> >(vals,"m"),
                        MPosition::ITRF);
  vals[0] = 3828713; vals[1] = 442878; vals[2] = 5064926;
  antPos[3] = MPosition(Quantum<Vector<double> >(vals,"m"),
                        MPosition::ITRF);
  Vector<double> antDiam(4, 70.);
  info.set (antNames, antDiam, antPos, ant1, ant2);
  return info;
}

void test1 (const DPInfo& info)
{
  cout << "test empty" << endl;
  ParameterSet ps;
  ps.add ("baseline", "[]");
  BaselineSelection selBL (ps, "");
  selBL.show (cout);
  Matrix<bool> res = selBL.apply (info);
  ASSERT (allEQ(res, true));
}

void test2 (const DPInfo& info)
{
  cout << "test auto" << endl;
  ParameterSet ps;
  ps.add ("baseline", "[]");
  ps.add ("corrtype", "AUTO");
  BaselineSelection selBL (ps, "");
  selBL.show (cout);
  Matrix<bool> res = selBL.apply (info);
  ASSERT (allEQ(res.diagonal(), true));
  res.diagonal() = false;
  ASSERT (allEQ(res, false));
}

void test3 (const DPInfo& info)
{
  cout << "test cross" << endl;
  ParameterSet ps;
  ps.add ("corrtype", "CROSS");
  BaselineSelection selBL (ps, "");
  selBL.show (cout);
  Matrix<bool> res = selBL.apply (info);
  ASSERT (allEQ(res.diagonal(), false));
  res.diagonal() = true;
  ASSERT (allEQ(res, true));
}

void test4 (const DPInfo& info)
{
  cout << "test length" << endl;
  Matrix<double> blength(4,4, 0.);
  blength(0,1) = blength(1,0) = 144.01;
  blength(0,2) = blength(2,0) = 288.021;
  blength(0,3) = blength(3,0) = 431.914;
  blength(1,2) = blength(2,1) = 144.01;
  blength(1,3) = blength(3,1) = 287.904;
  blength(2,3) = blength(3,2) = 143.896;
  ParameterSet ps;
  ps.add ("blmin", "145");
  BaselineSelection selBL (ps, "", true);
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), blength<=145.));
  ps.add ("blmax", "288");
  selBL = BaselineSelection(ps, "", true);
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), blength<=145. || blength>=288.));
  ps.add ("corrtype", "cross");
  selBL = BaselineSelection(ps, "", true);
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), (blength > 0.  &&
                                    (blength<=145. || blength>=288.))));
  ps.add ("blrange", "[0,144,288,430]");
  selBL = BaselineSelection(ps, "", true);
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), (blength > 0.  &&
                                    (blength<=145. || blength>=288.))));
  selBL = BaselineSelection(ps, "", false);
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), ((blength > 0.  &&  blength<=144.)  ||
                                    (blength>=288. &&  blength<=430.))));
  ps.replace ("corrtype", "");
  selBL = BaselineSelection(ps, "", false);
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), ((blength >=0.  &&  blength<=144.)  ||
                                    (blength>=288. &&  blength<=430.))));
}

void test5 (const DPInfo& info)
{
  cout << "test baseline" << endl;
  //antNames[0] = "rs01.s01";
  //antNames[1] = "rs02.s01";
  //antNames[2] = "cs01.s01";
  //antNames[3] = "cs01.s02";
  Matrix<int> bl(4,4);
  bl(0,0) = 00;
  bl(1,1) = 11;
  bl(2,2) = 22;
  bl(3,3) = 33;
  bl(0,1) = bl(1,0) = 01;
  bl(0,2) = bl(2,0) = 02;
  bl(0,3) = bl(3,0) = 03;
  bl(1,2) = bl(2,1) = 12;
  bl(1,3) = bl(3,1) = 13;
  bl(2,3) = bl(3,2) = 23;
  ParameterSet ps;
  ps.add ("baseline", "[[rs01.s01]]");
  BaselineSelection selBL (ps, "");
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), bl==00 || bl==01 || bl==02 || bl==03));
  ps.replace ("baseline", "[[rs01.s01,rs*]]");
  selBL = BaselineSelection(ps, "");
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), bl==00 || bl==01));
  ps.replace ("baseline", "[[rs01.s01],[rs02*]]");
  selBL = BaselineSelection(ps, "");
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), !(bl==22 || bl==23 || bl==33)));
  ps.replace ("baseline", "[rs01.s01,'rs02*']");
  selBL = BaselineSelection(ps, "");
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), !(bl==22 || bl==23 || bl==33)));
  ps.replace ("corrtype", "cross");
  selBL = BaselineSelection(ps, "");
  selBL.show (cout);
  ASSERT (allEQ(selBL.apply(info), !(bl==00 || bl==11 || bl==22 ||
                                     bl==23 || bl==33)));
  // Note that the MSSelection syntax is not tested, because it requires
  // a MeasurementSet.
}

int main()
{
  try {
    DPInfo info(makeInfo(16));
    test1(info);
    test2(info);
    test3(info);
    test4(info);
    test5(info);
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
