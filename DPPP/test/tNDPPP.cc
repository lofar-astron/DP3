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
#include <tables/Tables.h>
#include <tables/Tables/TableIter.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayPartMath.h>
#include <Common/LofarLogger.h>
#include <iostream>
#include <stdexcept>

using namespace casa;

void testCopy()
{
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msout=tNDPPP_tmp.MS1" << endl;
    ostr << "steps=[]" << endl;
  }
  system ("rm -rf tNDPPP_tmp.MS1; ../src/NDPPP tNDPPP_tmp.parset");
  Table tin("tNDPPP_tmp.MS");
  Table tout("tNDPPP_tmp.MS1");
  ASSERT (tout.nrow() == 6*20);
  {
    // Two dummy time slots were inserted, so ignore those.
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() not between 3*6 and 5*6-1", tout);
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
    Table t1 = tableCommand
      ("using style python "
       "select from $1 where rownumber() between 3*6 and 5*6-1", tout);
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
    double time = ROScalarColumn<double>(tin,"TIME")(2*6);
    for (uint i=0; i<12; ++i) {
      double timec = time + (1 + i/6)*30.;
      ASSERT (near(ROScalarColumn<double>(t1,"TIME")(i), timec));
      ASSERT (near(ROScalarColumn<double>(t1,"TIME_CENTROID")(i), timec));
    }
  }
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
}

void testAvg1()
{
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msin.countflags = true" << endl;
    ostr << "msout=tNDPPP_tmp.MS2" << endl;
    ostr << "msout.countflags = true" << endl;
    ostr << "steps=[avg]" << endl;
    ostr << "avg.type=average" << endl;
    ostr << "avg.timestep=20" << endl;
    ostr << "avg.freqstep=100" << endl;
  }
  system ("rm -rf tNDPPP_tmp.MS2; ../src/NDPPP tNDPPP_tmp.parset");
  checkAvg ("tNDPPP_tmp.MS2");
}

void testAvg2()
{
  // Averaging in multiple steps should be the same as above.
  {
    ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS" << endl;
    ostr << "msout=tNDPPP_tmp.MS3" << endl;
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
  system ("rm -rf tNDPPP_tmp.MS3; ../src/NDPPP tNDPPP_tmp.parset");
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
    ostr1 << "steps=[avg1]" << endl;
    ostr1 << "avg1.type=average" << endl;
    ostr1 << "avg1.timestep=5" << endl;
    ostr1 << "avg1.freqstep=2" << endl;
    ostr2 << "msin=tNDPPP_tmp.MS4a" << endl;
    ostr2 << "msout=tNDPPP_tmp.MS4b" << endl;
    ostr2 << "steps=[avg2]" << endl;
    ostr2 << "avg2.type=average" << endl;
    ostr2 << "avg2.timestep=1" << endl;
    ostr2 << "avg2.freqstep=2" << endl;
    ostr3 << "msin=tNDPPP_tmp.MS4b" << endl;
    ostr3 << "msout=tNDPPP_tmp.MS4c" << endl;
    ostr3 << "steps=[avg3]" << endl;
    ostr3 << "avg3.type=average" << endl;
    ostr3 << "avg3.timestep=2" << endl;
    ostr3 << "avg3.freqstep=1" << endl;
    ostr4 << "msin=tNDPPP_tmp.MS4c" << endl;
    ostr4 << "msout=tNDPPP_tmp.MS4d" << endl;
    ostr4 << "steps=[avg4]" << endl;
    ostr4 << "avg4.type=average" << endl;
    ostr4 << "avg4.timestep=2" << endl;
    ostr4 << "avg4.freqstep=4" << endl;
  }
  system ("rm -rf tNDPPP_tmp.MS4[abcd]; "
          "../src/NDPPP tNDPPP_tmp.parset1  && "
          "../src/NDPPP tNDPPP_tmp.parset2  && "
          "../src/NDPPP tNDPPP_tmp.parset3  && "
          "../src/NDPPP tNDPPP_tmp.parset4");
  checkAvg ("tNDPPP_tmp.MS4d");
}

int main()
{
  try
  {
    testCopy();
    testAvg1();
    testAvg2();
    testAvg3();
  } catch (std::exception& err) {
    std::cerr << "Error detected: " << err.what() << std::endl;
    return 1;
  }
}
