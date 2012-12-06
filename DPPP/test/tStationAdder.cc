//# tStationAdder.cc: Test program for class StationAdder
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
#include <DPPP/StationAdder.h>
#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <Common/StringUtil.h>
#include <Common/StreamUtil.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayIO.h>
#include <iostream>

using namespace LOFAR;
using namespace LOFAR::DPPP;
using namespace casa;
using namespace std;

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput: public DPInput
{
public:
  TestInput(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
  {
    info().init (ncorr, nchan, ntime, 0., 5., string());
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
    info().set (antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    Vector<double> chanWidth(nchan, 1000000.);
    Vector<double> chanFreqs(nchan);
    indgen (chanFreqs, 10500000., 1000000.);
    info().set (chanFreqs, chanWidth);
  }
private:
  virtual bool process (const DPBuffer&)
  {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(i+itsCount*10,i-10+itsCount*6);
    }
    Cube<Float> weights(itsNCorr, itsNChan, itsNBl);
    indgen (weights, 0.5f, 0.01f);
    Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }
    DPBuffer buf;
    buf.setTime (itsCount*30 + 4472025740.0);
    buf.setData (data);
    buf.setWeights (weights);
    buf.setUVW  (uvw);
    Cube<bool> flags(data.shape());
    flags = false;
    buf.setFlags (flags);
    Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = false;
    buf.setFullResFlags (fullResFlags);
    getNextStep()->process (buf);
    ++itsCount;
    return true;
  }

  virtual void finish() {getNextStep()->finish();}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo&) {}

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Class to check result of TestInput run by test1.
class TestOutput: public DPStep
{
public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
  {}
private:
  void addData (Cube<Complex>& to, const Cube<Complex>& from, int bl)
  {
    to += from(IPosition(3,0,0,bl), IPosition(3,to.nrow()-1,to.ncolumn()-1,bl));
  }
  void addConjData (Cube<Complex>& to, const Cube<Complex>& from, int bl)
  {
    to += conj(from(IPosition(3,0,0,bl), IPosition(3,to.nrow()-1,to.ncolumn()-1,bl)));
  }
  virtual bool process (const DPBuffer& buf)
  {
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(i+itsCount*10,i-10+itsCount*6);
    }
    Cube<Float> weights(itsNCorr, itsNChan, itsNBl);
    indgen (weights, 0.5f, 0.01f);
    Cube<Complex> databl0 (itsNCorr, itsNChan, 1);
    Cube<Complex> databl1 (itsNCorr, itsNChan, 1);
    addData (databl0, data, 0);
    addData (databl0, data, 5);
    addData (databl0, data, 15);
    addData (databl1, data, 8);
    addData (databl1, data, 9);
    addData (databl1, data, 11);
    addConjData (databl1, data, 2);
    addConjData (databl1, data, 6);
    addConjData (databl1, data, 14);
    Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }
    IPosition end(3,itsNCorr-1,itsNChan-1,itsNBl-1);
    ASSERT (allEQ (buf.getData()(IPosition(3,0), end), data));
    ASSERT (buf.getFlags().shape() == IPosition(3,itsNCorr,itsNChan,itsNBl+2));
    ASSERT (allEQ (buf.getFlags(), false));
    ASSERT (allEQ (buf.getWeights()(IPosition(3,0), end), weights));
    ASSERT (allEQ (buf.getUVW()(IPosition(2,0),
                                IPosition(2,2,itsNBl-1)), uvw));
    ASSERT (buf.getFullResFlags().shape() == IPosition(3,itsNChan,1,itsNBl+2));
    ASSERT (allEQ (buf.getFullResFlags(), false));
    // Now check data of new baselines.
    end[2] = itsNBl;
    ASSERT (allNear (buf.getData()(IPosition(3,0,0,itsNBl), end), databl0, 1e-5));
    ASSERT (allNear (buf.getWeights()(IPosition(3,0,0,itsNBl), end), 3.f, 1e-5));
    end[2] = itsNBl+1;
    ASSERT (allNear (buf.getData()(IPosition(3,0,0,itsNBl+1), end), databl1, 1e-5));
    ASSERT (allNear (buf.getWeights()(IPosition(3,0,0,itsNBl+1), end), 6.f, 1e-5));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& infoIn)
  {
    info() = infoIn;
    ASSERT (int(infoIn.origNChan())==itsNChan);
    ASSERT (int(infoIn.nchan())==itsNChan);
    ASSERT (int(infoIn.ntime())==itsNTime);
    ASSERT (infoIn.timeInterval()==5);
    ASSERT (int(infoIn.nchanAvg())==1);
    ASSERT (int(infoIn.ntimeAvg())==1);
    ASSERT (int(infoIn.nbaselines())==itsNBl+2);
    ASSERT (int(infoIn.antennaNames().size())==5);
    ASSERT (int(infoIn.antennaDiam().size())==5);
    ASSERT (int(infoIn.antennaPos().size())==5);
    ASSERT (infoIn.antennaNames()[4]=="ns");
    Vector<Double> pos1 (infoIn.antennaPos()[4].getValue().getValue());
    ASSERT (near(pos1[0], (3828763.+3828746.+3828713.)/3));
    ASSERT (near(pos1[1], ( 442449.+ 442592.+ 442878.)/3));
    ASSERT (near(pos1[2], (5064923.+5064924.+5064926.)/3));
    // Check diam.
    double d1 = sqrt ((pos1[0]-3828763) * (pos1[0]-3828763) +
                      (pos1[1]- 442449) * (pos1[1]- 442449) +
                      (pos1[2]-5064923) * (pos1[2]-5064923));
    double d2 = sqrt ((pos1[0]-3828746) * (pos1[0]-3828746) +
                      (pos1[1]- 442592) * (pos1[1]- 442592) +
                      (pos1[2]-5064924) * (pos1[2]-5064924));
    double d3 = sqrt ((pos1[0]-3828713) * (pos1[0]-3828713) +
                      (pos1[1]- 442878) * (pos1[1]- 442878) +
                      (pos1[2]-5064926) * (pos1[2]-5064926));
    ASSERT (near(infoIn.antennaDiam()[4], 70+2*max(d1,max(d2,d3))));
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsNAvgTime, itsNAvgChan;
};

// Class to check result of flagged, unaveraged TestInput run by test2.
class TestOutput2: public DPStep
{
public:
  TestOutput2(int ntime, int nbl, int nchan, int ncorr)
    : itsCount(0), itsNTime(ntime), itsNBl(nbl), itsNChan(nchan),
      itsNCorr(ncorr)
  {}
private:
  void addData (Cube<Complex>& to, const Cube<Complex>& from,
                Cube<Float>& tow, const Cube<Float>& weights, int bl)
  {
    to += from(IPosition(3,0,0,bl), IPosition(3,to.nrow()-1,to.ncolumn()-1,bl));
    tow += weights(IPosition(3,0,0,bl), IPosition(3,to.nrow()-1,to.ncolumn()-1,bl));
  }
  void addConjData (Cube<Complex>& to, const Cube<Complex>& from,
                    Cube<Float>& tow, const Cube<Float>& weights, int bl)
  {
    to += conj(from(IPosition(3,0,0,bl), IPosition(3,to.nrow()-1,to.ncolumn()-1,bl)));
    tow += weights(IPosition(3,0,0,bl), IPosition(3,to.nrow()-1,to.ncolumn()-1,bl));
  }
  virtual bool process (const DPBuffer& buf)
  {
    Cube<Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i=0; i<int(data.size()); ++i) {
      data.data()[i] = Complex(i+itsCount*10,i-10+itsCount*6);
    }
    Cube<Float> weights(itsNCorr, itsNChan, itsNBl);
    indgen (weights, 0.5f, 0.01f);
    Cube<Complex> databl0 (itsNCorr, itsNChan, 1);
    Cube<Complex> databl1 (itsNCorr, itsNChan, 1);
    Cube<Complex> databl2 (itsNCorr, itsNChan, 1);
    Cube<Complex> databl3 (itsNCorr, itsNChan, 1);
    Cube<Complex> databl4 (itsNCorr, itsNChan, 1);
    Cube<Float> weightbl0 (itsNCorr, itsNChan, 1, 0.);
    Cube<Float> weightbl1 (itsNCorr, itsNChan, 1, 0.);
    Cube<Float> weightbl2 (itsNCorr, itsNChan, 1, 0.);
    Cube<Float> weightbl3 (itsNCorr, itsNChan, 1, 0.);
    Cube<Float> weightbl4 (itsNCorr, itsNChan, 1, 0.);
    addData (databl0, data, weightbl0, weights, 8);
    addData (databl0, data, weightbl0, weights, 9);
    addData (databl1, data, weightbl1, weights, 12);
    addData (databl1, data, weightbl1, weights, 13);
    addData (databl2, data, weightbl2, weights, 2);
    addData (databl2, data, weightbl2, weights, 3);
    addData (databl3, data, weightbl3, weights, 6);
    addData (databl3, data, weightbl3, weights, 7);
    addConjData (databl0, data, weightbl0, weights, 2);
    addConjData (databl0, data, weightbl0, weights, 6);
    addConjData (databl1, data, weightbl1, weights, 3);
    addConjData (databl1, data, weightbl1, weights, 7);
    addConjData (databl2, data, weightbl2, weights, 8);
    addConjData (databl2, data, weightbl2, weights, 12);
    addConjData (databl3, data, weightbl3, weights, 9);
    addConjData (databl3, data, weightbl3, weights, 13);
    addConjData (databl4, databl0, weightbl4, weightbl0, 0);
    addConjData (databl4, databl1, weightbl4, weightbl1, 0);
    Matrix<double> uvw(3, itsNBl);
    for (int i=0; i<itsNBl; ++i) {
      uvw(0,i) = 1 + itsCount + i;
      uvw(1,i) = 2 + itsCount + i;
      uvw(2,i) = 3 + itsCount + i;
    }
    IPosition end(3,itsNCorr-1,itsNChan-1,itsNBl-1);
    ASSERT (allEQ (buf.getData()(IPosition(3,0), end), data));
    ASSERT (buf.getFlags().shape() == IPosition(3,itsNCorr,itsNChan,itsNBl+5));
    ASSERT (allEQ (buf.getFlags(), false));
    ASSERT (allEQ (buf.getWeights()(IPosition(3,0), end), weights));
    ASSERT (allEQ (buf.getUVW()(IPosition(2,0),
                                IPosition(2,2,itsNBl-1)), uvw));
    ASSERT (buf.getFullResFlags().shape() == IPosition(3,itsNChan,1,itsNBl+5));
    ASSERT (allEQ (buf.getFullResFlags(), false));
    // Now check data of new baselines.
    end[2] = itsNBl;
    cout<< buf.getUVW()(IPosition(2,0,itsNBl-1), IPosition(2,2,itsNBl+4));
    ASSERT (allNear (buf.getData()(IPosition(3,0,0,itsNBl), end), databl0, 1e-5));
    ASSERT (allNear (buf.getWeights()(IPosition(3,0,0,itsNBl), end), weightbl0, 1e-5));
    end[2] = itsNBl+1;
    ASSERT (allNear (buf.getData()(IPosition(3,0,0,itsNBl+1), end), databl1, 1e-5));
    ASSERT (allNear (buf.getWeights()(IPosition(3,0,0,itsNBl+1), end), weightbl1, 1e-5));
    end[2] = itsNBl+2;
    ASSERT (allNear (buf.getData()(IPosition(3,0,0,itsNBl+2), end), databl2, 1e-5));
    ASSERT (allNear (buf.getWeights()(IPosition(3,0,0,itsNBl+2), end), weightbl2, 1e-5));
    end[2] = itsNBl+3;
    ASSERT (allNear (buf.getData()(IPosition(3,0,0,itsNBl+3), end), databl3, 1e-5));
    ASSERT (allNear (buf.getWeights()(IPosition(3,0,0,itsNBl+3), end), weightbl3, 1e-5));
    end[2] = itsNBl+4;
    ASSERT (allNear (buf.getData()(IPosition(3,0,0,itsNBl+4), end), databl4, 1e-5));
    ASSERT (allNear (buf.getWeights()(IPosition(3,0,0,itsNBl+4), end), weightbl4, 1e-5));
    itsCount++;
    return true;
    ///cout << buf.getFlags() << endl << result << endl;
    ASSERT (allEQ(buf.getFlags(), false));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show (std::ostream&) const {}
  virtual void updateInfo (const DPInfo& infoIn)
  {
    info() = infoIn;
    ASSERT (int(infoIn.origNChan())==itsNChan);
    ASSERT (int(infoIn.nchan())==itsNChan);
    ASSERT (int(infoIn.ntime())==itsNTime);
    ASSERT (infoIn.timeInterval()==5);
    ASSERT (int(infoIn.nchanAvg())==1);
    ASSERT (int(infoIn.ntimeAvg())==1);
    ASSERT (int(infoIn.nbaselines())==itsNBl+5);
    ASSERT (int(infoIn.antennaNames().size())==6);
    ASSERT (infoIn.antennaNames()[4]=="ns1");
    ASSERT (infoIn.antennaNames()[5]=="ns2");
    Vector<Double> pos1 (infoIn.antennaPos()[4].getValue().getValue());
    ASSERT (near(pos1[0], (3828763.+3828746.)/2));
    ASSERT (near(pos1[1], ( 442449.+ 442592.)/2));
    ASSERT (near(pos1[2], (5064923.+5064924.)/2));
    Vector<Double> pos2 (infoIn.antennaPos()[5].getValue().getValue());
    ASSERT (near(pos2[0], (3828729.+3828713.)/2));
    ASSERT (near(pos2[1], ( 442735.+ 442878.)/2));
    ASSERT (near(pos2[2], (5064925.+5064926.)/2));
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsNAvgTime, itsNAvgChan;
};


// Execute steps.
void execute (const DPStep::ShPtr& step1)
{
  // Set DPInfo.
  step1->setInfo (DPInfo());
  step1->getNextStep()->show (cout);
  // Execute the steps.
  DPBuffer buf;
  while (step1->process(buf));
  step1->finish();
}

// Test adding 3 stations.
void test1(int ntime, int nbl, int nchan, int ncorr)
{
  cout << "test1: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("stations",
              "{ns:[rs01.s01, rs02.s01, cs01.s02]}");
  parset.add ("autocorr", "true");
  parset.add ("useweights", "false");
  DPStep::ShPtr step2(new StationAdder(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput(ntime, nbl, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
}

// Test adding two groups of 2 stations.
void test2(int ntime, int nbl, int nchan, int ncorr)
{
  cout << "test2: ntime=" << ntime << " nrbl=" << nbl << " nchan=" << nchan
       << " ncorr=" << ncorr << endl;
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("stations",
              "{ns1:[rs01.s01, rs02.s01], ns2:[cs01.s02, cs01.s01]}");
  parset.add ("autocorr", "false");
  DPStep::ShPtr step2(new StationAdder(in, parset, ""));
  DPStep::ShPtr step3(new TestOutput2(ntime, nbl, nchan, ncorr));
  step1->setNextStep (step2);
  step2->setNextStep (step3);
  execute (step1);
  step2->showCounts (cout);
}

void test3 (const string& stations)
{
  // Do some erronous attempts.
  TestInput* in = new TestInput(2, 8, 4, 4);
  DPStep::ShPtr step1(in);
  ParameterSet parset;
  parset.add ("stations", stations);
  parset.add ("autocorr", "true");
  DPStep::ShPtr step2(new StationAdder(in, parset, ""));
  step1->setNextStep (step2);
  bool ok = true;
  try {
    execute (step1);
  } catch (std::exception& x) {
    cout << "Expected exception: " << x.what() << endl;
    ok = false;
  }
  ASSERT (!ok);
}

void testPatterns()
{
  Vector<String> antNames(10);
  antNames[0] = "CS001HBA0";   antNames[1] = "CS001HBA1";
  antNames[2] = "CS002HBA0";   antNames[3] = "CS002HBA1";
  antNames[4] = "CS003HBA0";   antNames[5] = "CS003HBA1";
  antNames[6] = "CS004HBA0";   antNames[7] = "CS004HBA1";
  antNames[8] = "CS005HBA0";   antNames[9] = "CS005HBA1";
  vector<string> patterns;
  patterns.push_back ("CS00[0-9]*");
  cout << StationAdder::getMatchingStations (antNames, patterns) << endl;
  patterns[0] = "CS00[0-9]*";
  cout << StationAdder::getMatchingStations (antNames, patterns) << endl;
  patterns.push_back ("!CS00[45]*");
  cout << StationAdder::getMatchingStations (antNames, patterns) << endl;
  patterns.push_back ("CS00[124]HBA0");
  cout << StationAdder::getMatchingStations (antNames, patterns) << endl;
}


int main()
{
  INIT_LOGGER ("tUVWFlagger");
  try {
    // Test the station selection patterns.
    testPatterns();
    // Test must be done with with 16 baselines.
    test1( 10,  16, 32, 4);
    test2( 10,  16, 32, 4);
    // Unknown station.
    test3("{ns1:rs01.s1, ns2:[cs01.s02, cs01.s01]}");
    // New station already used.
    test3("{ns1:[rs01.s01, rs02.s01], cs01.s02:[cs01.s02, cs01.s01]}");
    // Old station doubly used.
    test3("{ns1:[rs01.s01, rs02.s01], ns2:[rs01.s01, cs01.s01]}");
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
