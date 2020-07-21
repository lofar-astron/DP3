// tDPPP.cc: test program for DPPP
// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$
//
// @author Ger van Diepen

#include <casacore/tables/Tables.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayPartMath.h>
#include <iostream>
#include <stdexcept>

#include <boost/test/unit_test.hpp>

#include "../../DPRun.h"

using namespace DP3::DPPP;
using namespace casacore;

BOOST_AUTO_TEST_SUITE(dppp)

// This test program uses the MS in tNDPPP.in_MS.tgz.
// The MS contains 4 corr, 16 freq, 6 baselines, 18 time slots of 30 sec.
// Two time slots are missing between time slot 2 and 3.

void checkCopy (const String& in, const String& out, int nms)
{
  Table tin(in);
  Table tout(out);
  BOOST_CHECK_EQUAL (tout.nrow(), size_t {6*24});
  for (int j=0; j<nms; ++j) {
    // A few dummy time slots were inserted, so ignore those.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() not in [0:6, 4*6:6*6, 21*6:24*6]",
       tout);
    ROArrayColumn<Complex> data(t1, "DATA");
    ROArrayColumn<Bool> flag(t1, "FLAG");
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    IPosition dshape (2, 4, 16);
    IPosition dst (2, 0, j*16);
    Slicer dslicer (dst, dshape);
    BOOST_CHECK_EQUAL (data(0).shape(), IPosition(2,4,16*nms));
    BOOST_CHECK_EQUAL (flag(0).shape(), IPosition(2,4,16*nms));
    BOOST_CHECK_EQUAL (oflag(0).shape(), IPosition(2,16*nms/8,1));
    BOOST_CHECK (allEQ(data.getColumn(dslicer),
                  ROArrayColumn<Complex>(tin,"DATA").getColumn()));
    BOOST_CHECK (allEQ(flag.getColumn(), false));
    BOOST_CHECK (allEQ(oflag.getColumn(), uChar(0)));
    BOOST_CHECK (allEQ(ROArrayColumn<float>(t1,"WEIGHT_SPECTRUM").getColumn(),
                  float(1)));
    //    cout<<ROArrayColumn<double>(t1,"UVW").getColumn()<<
    //                  ROArrayColumn<double>(tin,"UVW").getColumn();
    //    BOOST_CHECK (allEQ(ROArrayColumn<double>(t1,"UVW").getColumn(),
    //                  ROArrayColumn<double>(tin,"UVW").getColumn()));
    BOOST_CHECK (allEQ(ROScalarColumn<double>(t1,"TIME").getColumn(),
                  ROScalarColumn<double>(tin,"TIME").getColumn()));
    BOOST_CHECK (allEQ(ROScalarColumn<double>(t1,"TIME_CENTROID").getColumn(),
                  ROScalarColumn<double>(tin,"TIME_CENTROID").getColumn()));
    BOOST_CHECK (allEQ(ROScalarColumn<double>(t1,"INTERVAL").getColumn(),
                  ROScalarColumn<double>(tin,"INTERVAL").getColumn()));
    BOOST_CHECK (allEQ(ROScalarColumn<double>(t1,"EXPOSURE").getColumn(),
                  ROScalarColumn<double>(tin,"EXPOSURE").getColumn()));
  }
  {
    // Check the inserted time slots.
    // The MS misses a few time slots (3 and 4).
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [0:6, 4*6:6*6, 21*6:24*6]",
       tout);
    ROArrayColumn<Complex> data(t1, "DATA");
    ROArrayColumn<Bool> flag(t1, "FLAG");
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL (data(0).shape(), IPosition(2,4,16*nms));
    BOOST_CHECK_EQUAL (flag(0).shape(), IPosition(2,4,16*nms));
    BOOST_CHECK_EQUAL (oflag(0).shape(), IPosition(2,16*nms/8,1));
    BOOST_CHECK (allEQ(data.getColumn(), Complex()));
    BOOST_CHECK (allEQ(flag.getColumn(), true));
    BOOST_CHECK (allEQ(oflag.getColumn(), uChar(0xff)));
    BOOST_CHECK (allEQ(ROArrayColumn<float>(t1,"WEIGHT_SPECTRUM").getColumn(),
                  float(0)));
    BOOST_CHECK (allEQ(ROScalarColumn<double>(t1,"INTERVAL").getColumn(), 30.));
    BOOST_CHECK (allEQ(ROScalarColumn<double>(t1,"EXPOSURE").getColumn(), 30.));
    double time = ROScalarColumn<double>(tin,"TIME")(0);
    for (unsigned int i=0; i<6; ++i) {
      double timec = time - 30;
      BOOST_CHECK (near(ROScalarColumn<double>(t1,"TIME")(i), timec));
      BOOST_CHECK (near(ROScalarColumn<double>(t1,"TIME_CENTROID")(i), timec));
    }
    time = ROScalarColumn<double>(tin,"TIME")(2*6);
    for (unsigned int i=6; i<18; ++i) {
      double timec = time + (i/6)*30.;
      BOOST_CHECK (near(ROScalarColumn<double>(t1,"TIME")(i), timec));
      BOOST_CHECK (near(ROScalarColumn<double>(t1,"TIME_CENTROID")(i), timec));
    }
    time = ROScalarColumn<double>(tin,"TIME")(17*6);
    for (unsigned int i=18; i<36; ++i) {
      double timec = time + (i/6-2)*30.;
      BOOST_CHECK (near(ROScalarColumn<double>(t1,"TIME")(i), timec));
      BOOST_CHECK (near(ROScalarColumn<double>(t1,"TIME_CENTROID")(i), timec));
    }
  }
  // Now check if the SPECTRAL_WINDOW table is fine.
  Table spwin(tin.keywordSet().asTable("SPECTRAL_WINDOW"));
  Table spwout(tout.keywordSet().asTable("SPECTRAL_WINDOW"));
  for (int j=0; j<nms; ++j) {
    IPosition dshape (1, 16);
    IPosition dst (1, j*16);
    Slicer dslicer (dst, dshape);
    BOOST_CHECK (allEQ (ROArrayColumn<double>(spwin, "CHAN_FREQ").getColumn(),
                   ROArrayColumn<double>(spwout,"CHAN_FREQ").getColumn(dslicer)));
    BOOST_CHECK (allEQ (ROArrayColumn<double>(spwin, "CHAN_WIDTH").getColumn(),
                   ROArrayColumn<double>(spwout,"CHAN_WIDTH").getColumn(dslicer)));
    BOOST_CHECK (allEQ (ROArrayColumn<double>(spwin, "EFFECTIVE_BW").getColumn(),
                   ROArrayColumn<double>(spwout,"EFFECTIVE_BW").getColumn(dslicer)));
    BOOST_CHECK (allEQ (ROArrayColumn<double>(spwin, "RESOLUTION").getColumn(),
                   ROArrayColumn<double>(spwout,"RESOLUTION").getColumn(dslicer)));
  }
  BOOST_CHECK (allEQ (double(nms)*ROScalarColumn<double>(spwin, "TOTAL_BANDWIDTH").getColumn(),
                 ROScalarColumn<double>(spwout,"TOTAL_BANDWIDTH").getColumn()));
  if (nms == 1) {
    BOOST_CHECK (allEQ (ROScalarColumn<double>(spwin, "REF_FREQUENCY").getColumn(),
                   ROScalarColumn<double>(spwout,"REF_FREQUENCY").getColumn()));
  }
  // Check the TIME_RANGE in the OBSERVATION table.
  Table obsout(tout.keywordSet().asTable("OBSERVATION"));
  Vector<double> timeRange
    (ROArrayColumn<double>(obsout, "TIME_RANGE").getColumn());
  BOOST_CHECK (near(timeRange(0), ROScalarColumn<double>(tout,"TIME")(0) - 15));
  BOOST_CHECK (near(timeRange(1), ROScalarColumn<double>(tout,"TIME")(143) + 15));
}

void checkCopyColumn (const String& in)
{
  Table tin(in);
  BOOST_CHECK_EQUAL (tin.nrow(), size_t {6*24});
  ROArrayColumn<Complex> data1(tin, "DATA");
  ROArrayColumn<Complex> data2(tin, "COPY_DATA");
  ROArrayColumn<float> weight1(tin, "NEW_WEIGHT_SPECTRUM");
  ROArrayColumn<float> weight2(tin, "COPY_NEW_WEIGHT_SPECTRUM");
  BOOST_CHECK (allEQ(data1.getColumn(), data2.getColumn()));
  BOOST_CHECK (allEQ(weight1.getColumn(), weight2.getColumn()));
}

void testCopy()
{
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    // Give starttime 35 sec before actual, hence 1 missing timeslot.
    ostr << "msin.starttime=03-Aug-2000/13:21:45" << endl;
    // Give endtime 90 sec after actual, hence 3 missing timeslots.
    ostr << "msin.endtime=03-Aug-2000/13:33:15" << endl;
    ostr << "msout=tNDPPP_tmp.MS1" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  checkCopy ("tNDPPP_tmp.MS", "tNDPPP_tmp.MS1", 1);
}

void testCopyColumn()
{
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS1" << endl;
    ostr << "msout=." << endl;
    ostr << "msout.datacolumn=COPY_DATA" << endl;
    ostr << "msout.weightcolumn=NEW_WEIGHT_SPECTRUM" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");

  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS1" << endl;
    ostr << "msin.datacolumn=COPY_DATA" << endl;
    ostr << "msin.weightcolumn=NEW_WEIGHT_SPECTRUM" << endl;
    ostr << "msout=." << endl;
    ostr << "msout.datacolumn=DATA" << endl;
    ostr << "msout.weightcolumn=COPY_NEW_WEIGHT_SPECTRUM" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");

  checkCopyColumn ("tNDPPP_tmp.MS1");
}

void testMultiIn()
{
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=[tNDPPP_tmp.MS1, tNDPPP_tmp.MS1]" << endl;
    ostr << "msout=tNDPPP_tmp.MS1a" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  checkCopy ("tNDPPP_tmp.MS", "tNDPPP_tmp.MS1a", 2);
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=[tNDPPP_tmp.MS1, tNDPPP_tmp.MS1]" << endl;
    ostr << "msin.datacolumn=CORRECTED_DATA" << endl;
    ostr << "msin.weightcolumn=NEW_WEIGHT_SPECTRUM" << endl;
    ostr << "msin.missingdata=true" << endl;
    ostr << "msin.baseline=0,2&6" << endl;
    ostr << "msout=tNDPPP_tmp.MS1a" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  Table tab("tNDPPP_tmp.MS1a");
  BOOST_CHECK_EQUAL (tab.nrow(), size_t {48});
  BOOST_CHECK (allEQ (ROArrayColumn<Complex>(tab,"DATA").getColumn(), Complex()));
  BOOST_CHECK (tab.tableDesc().isColumn("WEIGHT_SPECTRUM"));
  BOOST_CHECK (allEQ (ROArrayColumn<Bool>(tab,"FLAG").getColumn(), True));
  BOOST_CHECK (allEQ (ROScalarColumn<Int>(tab,"ANTENNA2").getColumn(), 6));
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=[notexist, tNDPPP_tmp.MS1, notexist, notexist]" << endl;
    ostr << "msin.datacolumn=CORRECTED_DATA" << endl;
    ostr << "msin.missingdata=true" << endl;
    ostr << "msin.orderms=false" << endl;
    ostr << "msin.baseline=0,2&6" << endl;
    ostr << "msout=tNDPPP_tmp.MS1b" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  tab = Table("tNDPPP_tmp.MS1b");
  BOOST_CHECK_EQUAL (tab.nrow(), size_t {48});
  BOOST_CHECK (allEQ (ROArrayColumn<Complex>(tab,"DATA").getColumn(), Complex()));
  BOOST_CHECK (allEQ (ROArrayColumn<Bool>(tab,"FLAG").getColumn(), True));
  BOOST_CHECK (allEQ (ROScalarColumn<Int>(tab,"ANTENNA2").getColumn(), 6));
}

void checkAvg (const String& outName)
{
  Table tin("tNDPPP_tmp.MS");
  Table tout(outName);
  ROScalarColumn<double> timeCol(tin, "TIME");
  double time = 0.5 * (timeCol(tin.nrow()-1) + timeCol(0));
  BOOST_CHECK_EQUAL (tout.nrow(), size_t {6});
  Block<String> colNames(2);
  colNames[0] = "ANTENNA1";
  colNames[1] = "ANTENNA2";
  TableIterator iterin(tin, colNames);
  TableIterator iterout(tout, colNames);
  // Iterate over baseline to be able to average the input in an easy way.
  while (! iterin.pastEnd()) {
    Table t2(iterin.table());
    Table t1(iterout.table());
    ROArrayColumn<Complex> data(t1, "DATA");
    ROArrayColumn<Bool> flag(t1, "FLAG");
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL (data(0).shape(), IPosition(2,4,1));
    BOOST_CHECK_EQUAL (flag(0).shape(), IPosition(2,4,1));
    BOOST_CHECK_EQUAL (oflag(0).shape(), IPosition(2,16/8,20));
    // Average the original data over all channels and times.
    Array<Complex> dataAvg = partialMeans
      (ROArrayColumn<Complex>(t2,"DATA").getColumn(), IPosition(2,1,2));
    Array<Complex> dataRes = dataAvg.reform (IPosition(3,4,1,1));
    BOOST_CHECK (allNear(data.getColumn(), dataRes, 1e-5));
    BOOST_CHECK (allNear(ROArrayColumn<float>(t1,"WEIGHT_SPECTRUM").getColumn(),
                    float(18*16), 1e-5));
    BOOST_CHECK (allEQ(flag.getColumn(), false));
    BOOST_CHECK (allNear(ROScalarColumn<double>(t1,"TIME").getColumn(), time, 1e-5));
    BOOST_CHECK (allNear(ROScalarColumn<double>(t1,"TIME_CENTROID").getColumn(),
                    time, 1e-5));
    BOOST_CHECK (allEQ(ROScalarColumn<double>(t1,"INTERVAL").getColumn(), 20*30.));
    BOOST_CHECK (allEQ(ROScalarColumn<double>(t1,"EXPOSURE").getColumn(), 20*30.));
    // Two time entries should be flagged.
    Array<uChar> of = oflag.getColumn();
    BOOST_CHECK (allEQ(of(Slicer(IPosition(3,0,0,0),IPosition(3,2,3,1))),uChar(0)));
    BOOST_CHECK (allEQ(of(Slicer(IPosition(3,0,3,0),IPosition(3,2,2,1))),uChar(0xff)));
    BOOST_CHECK (allEQ(of(Slicer(IPosition(3,0,5,0),IPosition(3,2,15,1))),uChar(0)));
    iterin.next();
    iterout.next();
  }
  // Now check if the SPECTRAL_WINDOW table is fine.
  Table spwin(tin.keywordSet().asTable("SPECTRAL_WINDOW"));
  Table spwout(tout.keywordSet().asTable("SPECTRAL_WINDOW"));
  Matrix<double> cw = ROArrayColumn<double>(spwout, "CHAN_WIDTH").getColumn();
  BOOST_CHECK_EQUAL (cw.size(), size_t {1});
  BOOST_CHECK (near(cw(0,0),
               sum(ROArrayColumn<double>(spwin, "CHAN_WIDTH").getColumn())));
  Matrix<double> cfi = ROArrayColumn<double>(spwin,  "CHAN_FREQ").getColumn();
  Matrix<double> cfo = ROArrayColumn<double>(spwout, "CHAN_FREQ").getColumn();
  BOOST_CHECK (near(cfo(0,0), 0.5*(cfi(0,0) + cfi(15,0))));
  Matrix<double> ce = ROArrayColumn<double>(spwout, "EFFECTIVE_BW").getColumn();
  BOOST_CHECK_EQUAL (ce.size(), size_t {1});
  BOOST_CHECK (near(ce(0,0), cw(0,0)));
  Matrix<double> cr = ROArrayColumn<double>(spwout, "RESOLUTION").getColumn();
  BOOST_CHECK_EQUAL (cr.size(), size_t {1});
  BOOST_CHECK (near(cr(0,0), cw(0,0)));
  BOOST_CHECK (near(ROScalarColumn<double>(spwin, "TOTAL_BANDWIDTH")(0),
               cw(0,0)));
  BOOST_CHECK (allEQ (ROScalarColumn<double>(spwin, "REF_FREQUENCY").getColumn(),
                 ROScalarColumn<double>(spwout,"REF_FREQUENCY").getColumn()));
  // Check the TIME_RANGE in the OBSERVATION table.
  Table obsout(tout.keywordSet().asTable("OBSERVATION"));
  Vector<double> timeRange
    (ROArrayColumn<double>(obsout, "TIME_RANGE").getColumn());
  BOOST_CHECK (near(timeRange(0), ROScalarColumn<double>(tout,"TIME")(0) - 300));
  BOOST_CHECK (near(timeRange(1), ROScalarColumn<double>(tout,"TIME")(0) + 295));
}

void testAvg1()
{
  {
    // Average in a single step.
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin.name=tNDPPP_tmp.MS" << endl;
    // Give start and end time as actual, hence no missing timeslots.
    ostr << "msin.starttime=03-Aug-2000/13:22:20" << endl;
    ostr << "msin.endtime=03-Aug-2000/13:31:45" << endl;
    ostr << "msout.name=tNDPPP_tmp.MS2" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[avg,count]" << endl;
    ostr << "avg.type=average" << endl;
    ostr << "avg.timestep=20" << endl;
    ostr << "avg.freqstep=100" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  checkAvg ("tNDPPP_tmp.MS2");
}

void testAvg2()
{
  // Averaging in multiple steps should be the same as above.
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msout=tNDPPP_tmp.MS3" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[avg1,avg2,avg3,avg4]" << endl;
    ostr << "avg1.type=average" << endl;
    ostr << "avg1.timestep=5" << endl;
    ostr << "avg1.freqstep=2" << endl;
    ostr << "avg2.type=average" << endl;
    ostr << "avg2.timestep=1" << endl;
    ostr << "avg2.freqstep=2" << endl;
    ostr << "avg3.type=average" << endl;
    ostr << "avg3.timestep=2" << endl;
    ostr << "avg3.freqstep=1" << endl;
    ostr << "avg4.type=average" << endl;
    ostr << "avg4.timestep=2" << endl;
    ostr << "avg4.freqstep=4" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  checkAvg ("tNDPPP_tmp.MS3");
}

void testAvg3()
{
  // Averaging in multiple steps with multiple outputs should be the same
  // as above.
  {
    ofstream ostr1("tNDPPP_tmp.parset1");
    ofstream ostr2("tNDPPP_tmp.parset2");
    ofstream ostr3("tNDPPP_tmp.parset3");
    ofstream ostr4("tNDPPP_tmp.parset4");
    ostr1 << "msin=tNDPPP_tmp.MS" << endl;
    ostr1 << "msout=tNDPPP_tmp.MS4a" << endl;
    ostr1 << "msout.overwrite=true" << endl;
    ostr1 << "steps=[avg1]" << endl;
    ostr1 << "avg1.type=average" << endl;
    ostr1 << "avg1.timestep=5" << endl;
    ostr1 << "avg1.freqstep=2" << endl;
    ostr2 << "msin=tNDPPP_tmp.MS4a" << endl;
    ostr2 << "msout=tNDPPP_tmp.MS4b" << endl;
    ostr2 << "msout.overwrite=true" << endl;
    ostr2 << "steps=[avg2]" << endl;
    ostr2 << "avg2.type=average" << endl;
    ostr2 << "avg2.timestep=1" << endl;
    ostr2 << "avg2.freqstep=2" << endl;
    ostr3 << "msin=tNDPPP_tmp.MS4b" << endl;
    ostr3 << "msout=tNDPPP_tmp.MS4c" << endl;
    ostr3 << "msout.overwrite=true" << endl;
    ostr3 << "steps=[avg3]" << endl;
    ostr3 << "avg3.type=average" << endl;
    ostr3 << "avg3.timestep=2" << endl;
    ostr3 << "avg3.freqstep=1" << endl;
    ostr4 << "msin=tNDPPP_tmp.MS4c" << endl;
    ostr4 << "msout=tNDPPP_tmp.MS4d" << endl;
    ostr4 << "msout.overwrite=true" << endl;
    ostr4 << "steps=[avg4]" << endl;
    ostr4 << "avg4.type=average" << endl;
    ostr4 << "avg4.timestep=2" << endl;
    ostr4 << "avg4.freqstep=4" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset1");
  DPRun::execute ("tNDPPP_tmp.parset2");
  DPRun::execute ("tNDPPP_tmp.parset3");
  DPRun::execute ("tNDPPP_tmp.parset4");
  checkAvg ("tNDPPP_tmp.MS4d");
}

// This function tests if the correct start time is used when selecting times.
void testAvg4()
{
  {
    // Average in a single step.
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    // Give start a few seconds after first one, hence skip first time slot.
    ostr << "msin.starttime=03-Aug-2000/13:22:25" << endl;
    ostr << "msout=tNDPPP_tmp.MS5" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[avg]" << endl;
    ostr << "avg.type=average" << endl;
    ostr << "avg.timestep=2" << endl;
    ostr << "avg.freqstep=100" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  {
    // Only check the times;
    Table tab ("tNDPPP_tmp.MS");
    // First time to be used.
    double time = ROScalarColumn<double>(tab, "TIME")(6) + 15;
    Table t2 ("tNDPPP_tmp.MS5");
    BOOST_CHECK_EQUAL (t2.nrow(), size_t {6*10});
    ROScalarColumn<double> timeCol(t2, "TIME");
    for (unsigned int i=0; i<t2.nrow(); ++i) {
      BOOST_CHECK (near(timeCol(i), time));
      if (i%6 == 5) time += 60;
    }
  }
}

void testUpdate1()
{
  // Test if update works fine.
  // In fact, it does not do anything apart from rewriting the current flags.
  // However, it should ignore the inserted time slots.
  Array<bool> flags;
  {
    Table tab("tNDPPP_tmp.MS");
    flags = ROArrayColumn<bool>(tab,"FLAG").getColumn();
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msout=''" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  // Check that the flags did not change.
  {
    Table tab("tNDPPP_tmp.MS");
    BOOST_CHECK (allEQ(ROArrayColumn<bool>(tab,"FLAG").getColumn(), flags));
  }
}

void testUpdate2()
{
  // Test if update all flags works fine.
  {
    Table tab("tNDPPP_tmp.MS");
    tab.deepCopy ("tNDPPP_tmp.MS_copy1", Table::New);
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS_copy1" << endl;
    ostr << "msout=." << endl;
    ostr << "steps=[preflag]" << endl;
    ostr << "preflag.blmin=1e6" << endl;     // should flag all data
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  // Check that all flags are true.
  {
    Table tab("tNDPPP_tmp.MS_copy1");
    BOOST_CHECK (allEQ(ROArrayColumn<bool>(tab,"FLAG").getColumn(), true));
  }
}

void testUpdateScale()
{
  // Test if update data works fine.
  Array<Complex> data;
  Array<bool> flags;
  {
    Table tab("tNDPPP_tmp.MS");
    data = ROArrayColumn<Complex>(tab,"DATA").getColumn();
    flags = ROArrayColumn<bool>(tab,"FLAG").getColumn();
    tab.deepCopy ("tNDPPP_tmp.MS_copy1", Table::New);
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS_copy1" << endl;
    ostr << "msout=tNDPPP_tmp.MS_copy1" << endl;   // same name means update
    ostr << "steps=[scaledata]" << endl;
    ostr << "scaledata.coeffs=2" << endl;
    ostr << "scaledata.stations=*" << endl;
    ostr << "scaledata.scalesize=false" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  // Check that all data is doubled.
  {
    Table tab("tNDPPP_tmp.MS_copy1");
    data *= Complex(2,0);
    BOOST_CHECK (allNear(ROArrayColumn<Complex>(tab,"DATA").getColumn(), data, 1e-5));
    BOOST_CHECK (allEQ(ROArrayColumn<bool>(tab,"FLAG").getColumn(), flags));
  }
}

void checkFlags (const string& outName)
{
  // Only check the FULL_RES_FLAGS and table size.
  // The flags are created in various ways, but should be the same in all cases.
  // 3 time slots are averaged to 1.
  Table tout(outName);
  BOOST_CHECK_EQUAL (tout.nrow(), size_t {6*4});
  // Check the full-res-flags.
  // Channels 0,2,6,7,8 are flagged everywhere.
  {
  // Input time slots 3,4 are inserted, thus flagged.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [0:6]",
       tout);
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL (oflag(0).shape(), IPosition(2,2,6));
    Array<uChar> flags = oflag.getColumn();
    BOOST_CHECK (allEQ(flags(IPosition(3,0,0,0), IPosition(3,0,2,5)), uChar(0xc5)));
    BOOST_CHECK (allEQ(flags(IPosition(3,1,0,0), IPosition(3,1,2,5)), uChar(0x01)));
    BOOST_CHECK (allEQ(flags(IPosition(3,0,3,0), IPosition(3,0,4,5)), uChar(0xff)));
    BOOST_CHECK (allEQ(flags(IPosition(3,1,3,0), IPosition(3,1,4,5)), uChar(0x0f)));
    BOOST_CHECK (allEQ(flags(IPosition(3,0,5,0), IPosition(3,0,5,5)), uChar(0xc5)));
    BOOST_CHECK (allEQ(flags(IPosition(3,1,5,0), IPosition(3,1,5,5)), uChar(0x01)));
  }
  {
  // Input time slots 10,11 are flagged.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [6:12]",
       tout);
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL (oflag(0).shape(), IPosition(2,2,6));
    Array<uChar> flags = oflag.getColumn();
    BOOST_CHECK (allEQ(flags(IPosition(3,0,0,0), IPosition(3,0,3,5)), uChar(0xc5)));
    BOOST_CHECK (allEQ(flags(IPosition(3,1,0,0), IPosition(3,1,3,5)), uChar(0x01)));
    BOOST_CHECK (allEQ(flags(IPosition(3,0,4,0), IPosition(3,0,5,5)), uChar(0xff)));
    BOOST_CHECK (allEQ(flags(IPosition(3,1,4,0), IPosition(3,1,5,5)), uChar(0x0f)));
  }
  {
    // Input time slots 12-17 are not flagged.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [12:18]",
       tout);
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL (oflag(0).shape(), IPosition(2,2,6));
    Array<uChar> flags = oflag.getColumn();
    BOOST_CHECK (allEQ(flags(IPosition(3,0,0,0), IPosition(3,0,5,5)), uChar(0xc5)));
    BOOST_CHECK (allEQ(flags(IPosition(3,1,0,0), IPosition(3,1,5,5)), uChar(0x01)));
  }
  {
    // Input time slot 20-23 did not exist, thus flagged in average.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [18:24]",
       tout);
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL (oflag(0).shape(), IPosition(2,2,6));
    Array<uChar> flags = oflag.getColumn();
    BOOST_CHECK (allEQ(flags(IPosition(3,0,0,0), IPosition(3,0,1,5)), uChar(0xc5)));
    BOOST_CHECK (allEQ(flags(IPosition(3,1,0,0), IPosition(3,1,1,5)), uChar(0x01)));
    BOOST_CHECK (allEQ(flags(IPosition(3,0,2,0), IPosition(3,0,5,5)), uChar(0xff)));
    BOOST_CHECK (allEQ(flags(IPosition(3,1,2,0), IPosition(3,1,5,5)), uChar(0x0f)));
  }
}

void testFlags1()
{
  {
    // Most simple case.
    // Just flag some channels and time stamps.
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msin.startchan=1" << endl;
    ostr << "msin.nchan=12" << endl;
    ostr << "msout=tNDPPP_tmp.MS5" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[preflag,average]" << endl;
    ostr << "preflag.expr='flag1 or flag2'" << endl;
    ostr << "preflag.flag1.timeslot=[10,11]" << endl;
    ostr << "preflag.flag2.chan=[0,2,6..8]" << endl;
    ostr << "average.timestep=6" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  checkFlags ("tNDPPP_tmp.MS5");
}

void testFlags2()
{
  {
    // A more advanced case in two NDPPP steps.
    // Flag the channels, but shifted 2. In the next step the first 2 channels
    // are skipped, thus the channel numbers shift back 2.
    // An averaged time slot is flagged, so the original time slots should
    // come out flagged.
    ofstream ostr1("tNDPPP_tmp.parset1");
    ostr1 << "msin=tNDPPP_tmp.MS" << endl;
    ostr1 << "msin.nchan=15" << endl;
    ostr1 << "msout=tNDPPP_tmp.MS6a" << endl;
    ostr1 << "msout.overwrite=true" << endl;
    ostr1 << "steps=[preflag,average]" << endl;
    ostr1 << "preflag.chan=[2,4,8..10]" << endl;
    ostr1 << "average.timestep=2" << endl;
    ofstream ostr2("tNDPPP_tmp.parset2");
    ostr2 << "msin=tNDPPP_tmp.MS6a" << endl;
    ostr2 << "msin.startchan=2*1" << endl;    // output chan 0,2 are now flagged
    ostr2 << "msin.nchan=nchan-3" << endl;
    ostr2 << "msout=tNDPPP_tmp.MS6b" << endl;
    ostr2 << "msout.overwrite=true" << endl;
    ostr2 << "steps=[preflag,average]" << endl;
    ostr2 << "preflag.timeslot=5" << endl;  // is 10,11 in input
    ostr2 << "average.timestep=3" << endl;
    ostr2 << "average.freqstep=2" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset1");
  DPRun::execute ("tNDPPP_tmp.parset2");
  checkFlags ("tNDPPP_tmp.MS6b");
}

void testFlags3()
{
  {
    // Even a bit more advanced, also in two NDPPP runs.
    // Input channels 6,7,8 are flagged by flagging their averaged channel.
    // This is done at the end of run 1, so the averager of run 2 should pick
    // up those flags.
    ofstream ostr1("tNDPPP_tmp.parset1");
    ostr1 << "msin=tNDPPP_tmp.MS" << endl;
    ostr1 << "msin.nchan=15" << endl;
    ostr1 << "msout=tNDPPP_tmp.MS7a" << endl;
    ostr1 << "msout.overwrite=true" << endl;
    ostr1 << "steps=[preflag,average,pre2]" << endl;
    ostr1 << "preflag.chan=[0,2]" << endl;
    ostr1 << "average.timestep=2" << endl;
    ostr1 << "average.freqstep=3" << endl;
    ostr1 << "pre2.type=preflag" << endl;
    ostr1 << "pre2.chan=2" << endl;         // is input channel 6,7,8
    ofstream ostr2("tNDPPP_tmp.parset2");
    ostr2 << "msin=tNDPPP_tmp.MS7a" << endl;
    ostr2 << "msin.nchan=4" << endl;
    ostr2 << "msout=tNDPPP_tmp.MS7b" << endl;
    ostr2 << "msout.overwrite=true" << endl;
    ostr2 << "steps=[preflag,average]" << endl;
    ostr2 << "preflag.timeslot=5" << endl;  // is 10,11 in input
    ostr2 << "average.timestep=3" << endl;
    ostr2 << "average.freqstep=2" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset1");
  DPRun::execute ("tNDPPP_tmp.parset2");
  checkFlags ("tNDPPP_tmp.MS7b");
}

void testStationAdd()
{
  // Add station RT0, 1 and 2.
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msout=tNDPPP_tmp.MSa" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[stationadd]" << endl;
    ostr << "stationadd.stations={RTnew:[RT0..2]}" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  Table t1("tNDPPP_tmp.MS/ANTENNA");
  Table t2("tNDPPP_tmp.MSa/ANTENNA");
  BOOST_CHECK_EQUAL (t2.nrow(), t1.nrow()+1);     // 1 antenna has been added
  BOOST_CHECK_EQUAL (ROScalarColumn<String>(t2,"NAME")(t2.nrow()-1), "RTnew");
  Int oldNant = t1.nrow();
  t1 = Table("tNDPPP_tmp.MS/FEED");
  t2 = Table("tNDPPP_tmp.MSa/FEED");
  BOOST_CHECK_EQUAL (t2.nrow(), t1.nrow()+1);     // 1 antenna has been added
  BOOST_CHECK_EQUAL (ROScalarColumn<Int>(t2,"ANTENNA_ID")(t2.nrow()-1), oldNant);
  t1 = Table("tNDPPP_tmp.MS");
  t2 = Table("tNDPPP_tmp.MSa");
  BOOST_CHECK_EQUAL (t2.nrow(), t1.nrow()+40+12); // 2 baselines and 2 time slots added
}

void testFilter1()
{
  // Remove all baselines containing station RT1 or 6.
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msout=tNDPPP_tmp.MSa" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[filter]" << endl;
    ostr << "filter.baseline=!RT[16]&&*" << endl;
    ostr << "filter.remove=True" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  Table t1("tNDPPP_tmp.MS/ANTENNA");
  Table t2("tNDPPP_tmp.MSa/ANTENNA");
  // Note: the ANTENNA table also contained RT8, RT9, etc., but they do not
  // have baselines. So these were removed as well meaning only 0,2,7 are left.
  Vector<uInt> rownrs(3);
  rownrs[0]=0; rownrs[1]=2; rownrs[2]=7;
  Table t1s = t1(rownrs);
  BOOST_CHECK_EQUAL (t2.nrow(), t1s.nrow());
  BOOST_CHECK (allEQ (ROScalarColumn<String>(t2,"NAME").getColumn(),
                 ROScalarColumn<String>(t1s,"NAME").getColumn()));
  t1 = Table("tNDPPP_tmp.MS/FEED");
  t2 = Table("tNDPPP_tmp.MSa/FEED");
  t1s = t1(rownrs);
  BOOST_CHECK_EQUAL (t2.nrow(), t1s.nrow());
  // The ANTENNA_IDs in the FEED table must be 0,1,2.
  Vector<Int> ids(t2.nrow());
  indgen (ids);
  BOOST_CHECK (allEQ (ROScalarColumn<Int>(t2,"ANTENNA_ID").getColumn(), ids));
  // Check the main table.
  t1 = Table("tNDPPP_tmp.MS");
  t2 = Table("tNDPPP_tmp.MSa");
  BOOST_CHECK_EQUAL (t2.nrow(), t1.nrow()-72+4); // 4 baselines removed, 2 timeslots added
  t1s = t1((t1.col("ANTENNA1")==0 || t1.col("ANTENNA1")==2) &&
           t1.col("ANTENNA2")==7);
  // A few dummy time slots were inserted, so ignore those.
  Table t2s = t2(t2.nodeRownr() < 6  ||  t2.nodeRownr() >= 10);
  BOOST_CHECK (allEQ (ROArrayColumn<Complex>(t2s,"DATA").getColumn(),
                 ROArrayColumn<Complex>(t1s,"DATA").getColumn()));
  t2s = t2(t2.nodeRownr() % 2 == 0);
  BOOST_CHECK (allEQ (ROScalarColumn<Int>(t2s,"ANTENNA1").getColumn(), 0));
  BOOST_CHECK (allEQ (ROScalarColumn<Int>(t2s,"ANTENNA2").getColumn(), 2));
  t2s = t2(t2.nodeRownr() % 2 == 1);
  BOOST_CHECK (allEQ (ROScalarColumn<Int>(t2s,"ANTENNA1").getColumn(), 1));
  BOOST_CHECK (allEQ (ROScalarColumn<Int>(t2s,"ANTENNA2").getColumn(), 2));
}

void testFilter2()
{
  // Keep all baselines.
  // First by not specifying baseline selection, second by all baselines.
  // Also alter between remove and !remove.
  for (int iter=0; iter<4; ++iter) {
    {
      ofstream ostr("tNDPPP_tmp.parset");
      ostr << "msin=tNDPPP_tmp.MS" << endl;
      ostr << "msout=tNDPPP_tmp.MSa" << endl;
      ostr << "msout.overwrite=true" << endl;
      ostr << "steps=[filter]" << endl;
      if (iter%2 == 1) {
        ostr << "filter.baseline=*&&*" << endl;
      }
      if (iter/2 == 1) {
        ostr << "filter.remove=True" << endl;
      }
    }
    DPRun::execute ("tNDPPP_tmp.parset");
    //cout << "check ANTENNA"<<endl;
    Table t1("tNDPPP_tmp.MS/ANTENNA");
    Table t2("tNDPPP_tmp.MSa/ANTENNA");
    // Note: the ANTENNA table also contained RT8, RT9, etc., but they do not
    // have baselines. So these were removed meaning only 0,1,2,6,7 are left.
    Vector<uInt> rownrs(5);
    rownrs[0]=0; rownrs[1]=1; rownrs[2]=2; rownrs[3]=6; rownrs[4]=7;
    Table t1s(t1);
    if (iter/2 == 1) {
      t1s = t1(rownrs);
    }
    BOOST_CHECK_EQUAL (t2.nrow(), t1s.nrow());
    BOOST_CHECK (allEQ (ROScalarColumn<String>(t2,"NAME").getColumn(),
                   ROScalarColumn<String>(t1s,"NAME").getColumn()));
    //cout << "check FEED"<<endl;
    t1 = Table("tNDPPP_tmp.MS/FEED");
    t2 = Table("tNDPPP_tmp.MSa/FEED");
    t1s = t1;
    if (iter/2 == 1) {
      t1s = t1(rownrs);
    }
    BOOST_CHECK_EQUAL (t2.nrow(), t1s.nrow());
    // The ANTENNA_IDs in the FEED table must be 0,1,2.
    Vector<Int> ids(t2.nrow());
    indgen (ids);
    BOOST_CHECK (allEQ (ROScalarColumn<Int>(t2,"ANTENNA_ID").getColumn(), ids));
    // Check the main table.
    t1 = Table("tNDPPP_tmp.MS");
    t2 = Table("tNDPPP_tmp.MSa");
    BOOST_CHECK_EQUAL (t2.nrow(), t1.nrow()+12); // 2 timeslots added
    // A few dummy time slots were inserted, so ignore those.
    Table t2s = t2(t2.nodeRownr() < 18  ||  t2.nodeRownr() >= 30);
    BOOST_CHECK (allEQ (ROArrayColumn<Complex>(t2s,"DATA").getColumn(),
                   ROArrayColumn<Complex>(t1,"DATA").getColumn()));
    int ant1[] = {0,0,1,1,2,2};
    int ant2[] = {6,7,6,7,6,7};
    int sub = (iter/2 == 0 ? 0:3);    // if remove, ant2 6->3 and 7->4
    for (int i=0; i<6; ++i) {
      t2s = t2(t2.nodeRownr() % 6 == i);
      BOOST_CHECK (allEQ (ROScalarColumn<Int>(t2s,"ANTENNA1").getColumn(),
                     ant1[i]));
      BOOST_CHECK (allEQ (ROScalarColumn<Int>(t2s,"ANTENNA2").getColumn(),
                     ant2[i]-sub));
    }
  }
}


void testFilter3()
{
  // Remove some baselines, update original file with different data column
  // This test justs tests if it runs without throwing exceptions
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msout=." << endl;
    ostr << "msout.datacolumn=DATA_FILTER" << endl;
    ostr << "steps=[filter]" << endl;
    ostr << "filter.baseline=!RT[16]&&*" << endl;
    ostr << "filter.remove=False" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
}


void testClear()
{
  Array<bool> flags;
  // First flag in the same way as testFlags1.
  testFlags1();
  // Get the resulting flags.
  {
    Table tab("tNDPPP_tmp.MS5");
    flags.reference (ROArrayColumn<bool>(tab, "FLAG").getColumn());
  }
  // Flag all data.
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS5" << endl;
    ostr << "msout=tNDPPP_tmp.MS5a" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[preflag]" << endl;
    ostr << "preflag.baseline=[[*]]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  checkFlags ("tNDPPP_tmp.MS5a");
  {
    Table tab("tNDPPP_tmp.MS5a");
    BOOST_CHECK (allEQ(ROArrayColumn<bool>(tab, "FLAG").getColumn(), True));
  }
  // Clear the flags.
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS5a" << endl;
    ostr << "msout=tNDPPP_tmp.MS5b" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[preflag]" << endl;
    ostr << "preflag.mode=clear" << endl;
    ostr << "preflag.baseline=[[*]]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  checkFlags ("tNDPPP_tmp.MS5b");
  {
    Table tab("tNDPPP_tmp.MS5b");
    BOOST_CHECK (allEQ(ROArrayColumn<bool>(tab, "FLAG").getColumn(), false));
  }
}


void testMultiOut()
{
  {
    // First make the reference output MS.
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msout=tNDPPP_tmp.MS_copy" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  // Test if update data works fine with multiple outputs:
  // read from tNDPPP_tmp.MS, write to copy3, update to copy3
  Array<Complex> data;
  Array<bool> flags;
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "steps=[scaledata,out1,scaledata,out2]" << endl;
    ostr << "scaledata.coeffs=2" << endl;
    ostr << "scaledata.stations=*" << endl;
    ostr << "scaledata.scalesize=false" << endl;
    ostr << "out1.type=out" << endl;
    ostr << "out1.name=tNDPPP_tmp.MS_copy3" << endl;
    ostr << "out1.overwrite=true" << endl;
    ostr << "out2.type=out" << endl;
    ostr << "out2.name=." << endl; // Defaults to the previous out, so _copy3
    ostr << "out2.datacolumn=DATA_2" << endl;
    ostr << "msout=tNDPPP_tmp.MS_copy4" << endl;   // same name means update
    ostr << "msout.overwrite=true" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  // Check that tables exist, contain the specified columns
  {
    Table tab1("tNDPPP_tmp.MS_copy");
    Table tab2("tNDPPP_tmp.MS_copy3");
    ///cout<<ROArrayColumn<Complex>(tab2,"DATA_2").getColumn();
    ///cout<<Complex(2,0)*ROArrayColumn<Complex>(tab1,"DATA").getColumn();
    BOOST_CHECK (allNear(ROArrayColumn<Complex>(tab2,"DATA").getColumn(),
                    Complex(2,0)*ROArrayColumn<Complex>(tab1,"DATA").getColumn(),
                    1e-5));
    BOOST_CHECK (allNear(ROArrayColumn<Complex>(tab2,"DATA_2").getColumn(),
                    Complex(4,0)*ROArrayColumn<Complex>(tab1,"DATA").getColumn(),
                    1e-5));
    BOOST_CHECK (allNear(ROArrayColumn<Float>(tab2,"WEIGHT_SPECTRUM").getColumn(),
                    ROArrayColumn<Float>(tab1,"WEIGHT_SPECTRUM").getColumn(),
                    1e-5));
    BOOST_CHECK (allEQ(ROArrayColumn<Bool>(tab2,"FLAG").getColumn(),
                  ROArrayColumn<Bool>(tab1,"FLAG").getColumn()));
    Table tab3("tNDPPP_tmp.MS_copy4");
    BOOST_CHECK (allNear(ROArrayColumn<Complex>(tab3,"DATA").getColumn(),
                    Complex(4,0)*ROArrayColumn<Complex>(tab1,"DATA").getColumn(),
                    1e-5));
    BOOST_CHECK (allNear(ROArrayColumn<Float>(tab3,"WEIGHT_SPECTRUM").getColumn(),
                    ROArrayColumn<Float>(tab1,"WEIGHT_SPECTRUM").getColumn(),
                    1e-5));
    BOOST_CHECK (allEQ(ROArrayColumn<Bool>(tab3,"FLAG").getColumn(),
                  ROArrayColumn<Bool>(tab1,"FLAG").getColumn()));
  }

}

void tryErr (const string& parsetName)
{
  bool err = false;
  try {
    DPRun::execute (parsetName);
  } catch (const std::exception& x) {
    err = true;
  }
  BOOST_CHECK (err);
}

void testErrorOut()
{
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "steps=[filter,out1,average,out2]" << endl;
    ostr << "out1.type=out" << endl;
    ostr << "out1.name=''" << endl;
    ostr << "out2.type=out" << endl;
    ostr << "out2.name=." << endl;       // update not possible when avg
    ostr << "msout=''" << endl;
    tryErr ("tNDPPP_tmp.parset");
  }
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "steps=[average,out1,filter,out2]" << endl;
    ostr << "out1.type=out" << endl;
    ostr << "out1.name=tNDPPP_tmp.MSx" << endl;
    ostr << "out1.overwrite=true" << endl;
    ostr << "filter.remove=true" << endl;
    ostr << "out2.type=out" << endl;
    ostr << "out2.name=./tNDPPP_tmp.MSx" << endl;  // update not possible (filter)
    ostr << "msout=''" << endl;
    tryErr ("tNDPPP_tmp.parset");
  }
}

BOOST_AUTO_TEST_CASE( test_copy ) {
  testCopy();
}

BOOST_AUTO_TEST_CASE( test_copy_column ) {
  testCopyColumn();
}

BOOST_AUTO_TEST_CASE( test_multi_in ) {
  testMultiIn();
}

BOOST_AUTO_TEST_CASE( test_avg1 ) {
  testAvg1();
}

BOOST_AUTO_TEST_CASE( test_avg2 ) {
  testAvg2();
}

BOOST_AUTO_TEST_CASE( test_avg3 ) {
  testAvg3();
}

BOOST_AUTO_TEST_CASE( test_avg4 ) {
  testAvg4();
}

BOOST_AUTO_TEST_CASE( test_update1 ) {
  testUpdate1();
}

BOOST_AUTO_TEST_CASE( test_update2 ) {
  testUpdate2();
}

BOOST_AUTO_TEST_CASE( test_update_scale ) {
  testUpdateScale();
}

BOOST_AUTO_TEST_CASE( test_flags1 ) {
  testFlags1();
}

BOOST_AUTO_TEST_CASE( test_flags2 ) {
  testFlags2();
}

BOOST_AUTO_TEST_CASE( test_flags3 ) {
  testFlags3();
}

BOOST_AUTO_TEST_CASE( test_station_add ) {
  testStationAdd();
}

BOOST_AUTO_TEST_CASE( test_filter1 ) {
  testFilter1();
}

BOOST_AUTO_TEST_CASE( test_filter2 ) {
  testFilter2();
}

BOOST_AUTO_TEST_CASE( test_filter3 ) {
  testFilter3();
}

BOOST_AUTO_TEST_CASE( test_clear ) {
  testClear();
}

BOOST_AUTO_TEST_CASE( test_multi_out ) {
  testMultiOut();
}

BOOST_AUTO_TEST_CASE( test_error_out ) {
  testErrorOut();
}

BOOST_AUTO_TEST_SUITE_END()
