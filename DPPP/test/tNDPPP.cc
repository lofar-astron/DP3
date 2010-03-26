//# tNDPPP.cc: test program for NDPPP
//# Copyright (C) 2010
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
#include <DPPP/DPRun.h>
#include <tables/Tables.h>
#include <tables/Tables/TableIter.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayPartMath.h>
#include <Common/LofarLogger.h>
#include <iostream>
#include <stdexcept>

using namespace LOFAR::DPPP;
using namespace casa;

// This test program uses the MS in tNDPPP.in_MS.tgz.
// The MS contains 4 corr, 16 freq, 6 baselines, 18 time slots of 30 sec.
// Two time slots are missing between time slot 2 and 3.

void testCopy()
{
  cout << endl << "** testCopy **" << endl;
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    // Give starttime 30 sec before actual, hence 1 missing timeslot.
    ostr << "msin.starttime=03-Aug-2000/13:21:50" << endl;
    // Give endtime 90 sec after actual, hence 3 missing timeslots.
    ostr << "msin.endtime=03-Aug-2000/13:33:15" << endl;
    ostr << "msout=tNDPPP_tmp.MS1" << endl;
    ostr << "msout.overwrite=true" << endl;
    ostr << "steps=[]" << endl;
  }
  DPRun::execute ("tNDPPP_tmp.parset");
  Table tin("tNDPPP_tmp.MS");
  Table tout("tNDPPP_tmp.MS1");
  ASSERT (tout.nrow() == 6*24);
  {
    // A few dummy time slots were inserted, so ignore those.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() not in [0:6, 4*6:6*6, 21*6:24*6]",
       tout);
    ROArrayColumn<Complex> data(t1, "DATA");
    ROArrayColumn<Bool> flag(t1, "FLAG");
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    ASSERT (data(0).shape() == IPosition(2,4,16));
    ASSERT (flag(0).shape() == IPosition(2,4,16));
    ASSERT (oflag(0).shape() == IPosition(2,16/8,1));
    ASSERT (allEQ(data.getColumn(),
                  ROArrayColumn<Complex>(tin,"DATA").getColumn()));
    ASSERT (allEQ(flag.getColumn(), false));
    ASSERT (allEQ(oflag.getColumn(), uChar(0)));
    ASSERT (allEQ(ROArrayColumn<float>(t1,"WEIGHT_SPECTRUM").getColumn(),
                  float(1)));
    ASSERT (allEQ(ROArrayColumn<double>(t1,"UVW").getColumn(),
                  ROArrayColumn<double>(tin,"UVW").getColumn()));
    ASSERT (allEQ(ROScalarColumn<double>(t1,"TIME").getColumn(),
                  ROScalarColumn<double>(tin,"TIME").getColumn()));
    ASSERT (allEQ(ROScalarColumn<double>(t1,"TIME_CENTROID").getColumn(),
                  ROScalarColumn<double>(tin,"TIME_CENTROID").getColumn()));
    ASSERT (allEQ(ROScalarColumn<double>(t1,"INTERVAL").getColumn(),
                  ROScalarColumn<double>(tin,"INTERVAL").getColumn()));
    ASSERT (allEQ(ROScalarColumn<double>(t1,"EXPOSURE").getColumn(),
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
    ASSERT (data(0).shape() == IPosition(2,4,16));
    ASSERT (flag(0).shape() == IPosition(2,4,16));
    ASSERT (oflag(0).shape() == IPosition(2,16/8,1));
    ASSERT (allEQ(data.getColumn(), Complex()));
    ASSERT (allEQ(flag.getColumn(), true));
    ASSERT (allEQ(oflag.getColumn(), uChar(0xff)));
    ASSERT (allEQ(ROArrayColumn<float>(t1,"WEIGHT_SPECTRUM").getColumn(),
                  float(0)));
    ASSERT (allEQ(ROScalarColumn<double>(t1,"INTERVAL").getColumn(), 30.));
    ASSERT (allEQ(ROScalarColumn<double>(t1,"EXPOSURE").getColumn(), 30.));
    double time = ROScalarColumn<double>(tin,"TIME")(0);
    for (uint i=0; i<6; ++i) {
      double timec = time - 30;
      ASSERT (near(ROScalarColumn<double>(t1,"TIME")(i), timec));
      ASSERT (near(ROScalarColumn<double>(t1,"TIME_CENTROID")(i), timec));
    }
    time = ROScalarColumn<double>(tin,"TIME")(2*6);
    for (uint i=6; i<18; ++i) {
      double timec = time + (i/6)*30.;
      ASSERT (near(ROScalarColumn<double>(t1,"TIME")(i), timec));
      ASSERT (near(ROScalarColumn<double>(t1,"TIME_CENTROID")(i), timec));
    }
    time = ROScalarColumn<double>(tin,"TIME")(17*6);
    for (uint i=18; i<36; ++i) {
      double timec = time + (i/6-2)*30.;
      ASSERT (near(ROScalarColumn<double>(t1,"TIME")(i), timec));
      ASSERT (near(ROScalarColumn<double>(t1,"TIME_CENTROID")(i), timec));
    }
  }
  // Now check if the SPECTRAL_WINDOW table is fine.
  Table spwin(tin.keywordSet().asTable("SPECTRAL_WINDOW"));
  Table spwout(tout.keywordSet().asTable("SPECTRAL_WINDOW"));
  ASSERT (allEQ (ROArrayColumn<double>(spwin, "CHAN_FREQ").getColumn(),
                 ROArrayColumn<double>(spwout,"CHAN_FREQ").getColumn()));
  ASSERT (allEQ (ROArrayColumn<double>(spwin, "CHAN_WIDTH").getColumn(),
                 ROArrayColumn<double>(spwout,"CHAN_WIDTH").getColumn()));
  ASSERT (allEQ (ROArrayColumn<double>(spwin, "EFFECTIVE_BW").getColumn(),
                 ROArrayColumn<double>(spwout,"EFFECTIVE_BW").getColumn()));
  ASSERT (allEQ (ROArrayColumn<double>(spwin, "RESOLUTION").getColumn(),
                 ROArrayColumn<double>(spwout,"RESOLUTION").getColumn()));
  ASSERT (allEQ (ROScalarColumn<double>(spwin, "TOTAL_BANDWIDTH").getColumn(),
                 ROScalarColumn<double>(spwout,"TOTAL_BANDWIDTH").getColumn()));
  ASSERT (allEQ (ROScalarColumn<double>(spwin, "REF_FREQUENCY").getColumn(),
                 ROScalarColumn<double>(spwout,"REF_FREQUENCY").getColumn()));
  // Check the TIME_RANGE in the OBSERVATION table.
  Table obsout(tout.keywordSet().asTable("OBSERVATION"));
  Vector<double> timeRange
    (ROArrayColumn<double>(obsout, "TIME_RANGE").getColumn());
  ASSERT (near(timeRange(0), ROScalarColumn<double>(tout,"TIME")(0) - 15));
  ASSERT (near(timeRange(1), ROScalarColumn<double>(tout,"TIME")(143) + 15));
}

void checkAvg (const String& outName)
{
  Table tin("tNDPPP_tmp.MS");
  Table tout(outName);
  ROScalarColumn<double> timeCol(tin, "TIME");
  double time = 0.5 * (timeCol(tin.nrow()-1) + timeCol(0));
  ASSERT (tout.nrow() == 6);
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
    ASSERT (data(0).shape() == IPosition(2,4,1));
    ASSERT (flag(0).shape() == IPosition(2,4,1));
    ASSERT (oflag(0).shape() == IPosition(2,16/8,20));
    // Average the original data over all channels and times.
    Array<Complex> dataAvg = partialMeans
      (ROArrayColumn<Complex>(t2,"DATA").getColumn(), IPosition(2,1,2));
    Array<Complex> dataRes = dataAvg.reform (IPosition(3,4,1,1));
    ASSERT (allNear(data.getColumn(), dataRes, 1e-5));
    ASSERT (allNear(ROArrayColumn<float>(t1,"WEIGHT_SPECTRUM").getColumn(),
                    float(18*16), 1e-5));
    ASSERT (allEQ(flag.getColumn(), false));
    ASSERT (allNear(ROScalarColumn<double>(t1,"TIME").getColumn(), time, 1e-5));
    ASSERT (allNear(ROScalarColumn<double>(t1,"TIME_CENTROID").getColumn(),
                    time, 1e-5));
    ASSERT (allEQ(ROScalarColumn<double>(t1,"INTERVAL").getColumn(), 20*30.));
    ASSERT (allEQ(ROScalarColumn<double>(t1,"EXPOSURE").getColumn(), 20*30.));
    // Two time entries should be flagged.
    Array<uChar> of = oflag.getColumn();
    ASSERT (allEQ(of(Slicer(IPosition(3,0,0,0),IPosition(3,2,3,1))),uChar(0)));
    ASSERT (allEQ(of(Slicer(IPosition(3,0,3,0),IPosition(3,2,2,1))),uChar(0xff)));
    ASSERT (allEQ(of(Slicer(IPosition(3,0,5,0),IPosition(3,2,15,1))),uChar(0)));
    iterin.next();
    iterout.next();
  }
  // Now check if the SPECTRAL_WINDOW table is fine.
  Table spwin(tin.keywordSet().asTable("SPECTRAL_WINDOW"));
  Table spwout(tout.keywordSet().asTable("SPECTRAL_WINDOW"));
  Matrix<double> cw = ROArrayColumn<double>(spwout, "CHAN_WIDTH").getColumn();
  ASSERT (cw.size() == 1);
  ASSERT (near(cw(0,0),
               sum(ROArrayColumn<double>(spwin, "CHAN_WIDTH").getColumn())));
  Matrix<double> cfi = ROArrayColumn<double>(spwin,  "CHAN_FREQ").getColumn();
  Matrix<double> cfo = ROArrayColumn<double>(spwout, "CHAN_FREQ").getColumn();
  ASSERT (near(cfo(0,0), 0.5*(cfi(0,0) + cfi(15,0))));
  Matrix<double> ce = ROArrayColumn<double>(spwout, "EFFECTIVE_BW").getColumn();
  ASSERT (ce.size() == 1);
  ASSERT (near(ce(0,0), cw(0,0)));
  Matrix<double> cr = ROArrayColumn<double>(spwout, "RESOLUTION").getColumn();
  ASSERT (cr.size() == 1);
  ASSERT (near(cr(0,0), cw(0,0)));
  ASSERT (near(ROScalarColumn<double>(spwin, "TOTAL_BANDWIDTH")(0),
               cw(0,0)));
  ASSERT (allEQ (ROScalarColumn<double>(spwin, "REF_FREQUENCY").getColumn(),
                 ROScalarColumn<double>(spwout,"REF_FREQUENCY").getColumn()));
  // Check the TIME_RANGE in the OBSERVATION table.
  Table obsout(tout.keywordSet().asTable("OBSERVATION"));
  Vector<double> timeRange
    (ROArrayColumn<double>(obsout, "TIME_RANGE").getColumn());
  ASSERT (near(timeRange(0), ROScalarColumn<double>(tout,"TIME")(0) - 300));
  ASSERT (near(timeRange(1), ROScalarColumn<double>(tout,"TIME")(0) + 295));
}

void testAvg1()
{
  cout << endl << "** testAvg1 **" << endl;
  {
    // Average in a single step.
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    // Give start and end time as actual, hence no missing timeslots.
    ostr << "msin.starttime=03-Aug-2000/13:22:20" << endl;
    ostr << "msin.endtime=03-Aug-2000/13:31:45" << endl;
    ostr << "msout=tNDPPP_tmp.MS2" << endl;
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
  cout << endl << "** testAvg2 **" << endl;
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
  cout << endl << "** testAvg3 **" << endl;
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

void testUpdate()
{
  cout << endl << "** testUpdate **" << endl;
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
    ASSERT (allEQ(ROArrayColumn<bool>(tab,"FLAG").getColumn(), flags));
  }
}

void checkFlags (const string& outName)
{
  // Only check the FULL_RES_FLAGS and table size.
  // The flags are created in various ways, but should be the same in all cases.
  // 3 time slots are averaged to 1.
  Table tout(outName);
  ASSERT (tout.nrow() == 6*4);
  // Check the full-res-flags.
  // Channels 0,2,6,7,8 are flagged everywhere.
  {
  // Input time slots 3,4 are inserted, thus flagged.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [0:6]",
       tout);
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    ASSERT (oflag(0).shape() == IPosition(2,2,6));
    Array<uChar> flags = oflag.getColumn();
    ASSERT (allEQ(flags(IPosition(3,0,0,0), IPosition(3,0,2,5)), uChar(0xc5)));
    ASSERT (allEQ(flags(IPosition(3,1,0,0), IPosition(3,1,2,5)), uChar(0x01)));
    ASSERT (allEQ(flags(IPosition(3,0,3,0), IPosition(3,0,4,5)), uChar(0xff)));
    ASSERT (allEQ(flags(IPosition(3,1,3,0), IPosition(3,1,4,5)), uChar(0x0f)));
    ASSERT (allEQ(flags(IPosition(3,0,5,0), IPosition(3,0,5,5)), uChar(0xc5)));
    ASSERT (allEQ(flags(IPosition(3,1,5,0), IPosition(3,1,5,5)), uChar(0x01)));
  }
  {
  // Input time slots 10,11 are flagged.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [6:12]",
       tout);
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    ASSERT (oflag(0).shape() == IPosition(2,2,6));
    Array<uChar> flags = oflag.getColumn();
    cout << flags<<endl;
    ASSERT (allEQ(flags(IPosition(3,0,0,0), IPosition(3,0,3,5)), uChar(0xc5)));
    ASSERT (allEQ(flags(IPosition(3,1,0,0), IPosition(3,1,3,5)), uChar(0x01)));
    ASSERT (allEQ(flags(IPosition(3,0,4,0), IPosition(3,0,5,5)), uChar(0xff)));
    ASSERT (allEQ(flags(IPosition(3,1,4,0), IPosition(3,1,5,5)), uChar(0x0f)));
  }
  {
    // Input time slots 12-17 are not flagged.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [12:18]",
       tout);
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    ASSERT (oflag(0).shape() == IPosition(2,2,6));
    Array<uChar> flags = oflag.getColumn();
    ASSERT (allEQ(flags(IPosition(3,0,0,0), IPosition(3,0,5,5)), uChar(0xc5)));
    ASSERT (allEQ(flags(IPosition(3,1,0,0), IPosition(3,1,5,5)), uChar(0x01)));
  }
  {
    // Input time slot 20-23 did not exist, thus flagged in average.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() in [18:24]",
       tout);
    ROArrayColumn<uChar> oflag(t1, "LOFAR_FULL_RES_FLAG");
    ASSERT (oflag(0).shape() == IPosition(2,2,6));
    Array<uChar> flags = oflag.getColumn();
    ASSERT (allEQ(flags(IPosition(3,0,0,0), IPosition(3,0,1,5)), uChar(0xc5)));
    ASSERT (allEQ(flags(IPosition(3,1,0,0), IPosition(3,1,1,5)), uChar(0x01)));
    ASSERT (allEQ(flags(IPosition(3,0,2,0), IPosition(3,0,5,5)), uChar(0xff)));
    ASSERT (allEQ(flags(IPosition(3,1,2,0), IPosition(3,1,5,5)), uChar(0x0f)));
  }
}

void testFlags1()
{
  cout << endl << "** testFlags1 **" << endl;
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
  cout << endl << "** testFlags2 **" << endl;
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
    ostr2 << "msin.startchan=2" << endl;    // output chan 0,2 are now flagged
    ostr2 << "msin.nchan=12" << endl;
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
  cout << endl << "** testFlags3 **" << endl;
  {
    // Even a bit more advanced, also in two NDPPP runs.
    // Input channels 6,7,8 are flagged by flagging their avergaed channel.
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


int main()
{
  try
  {
    testCopy();
    testAvg1();
    testAvg2();
    testAvg3();
    testUpdate();
    testFlags1();
    testFlags2();
  } catch (std::exception& err) {
    std::cerr << "Error detected: " << err.what() << std::endl;
    return 1;
  }
}
