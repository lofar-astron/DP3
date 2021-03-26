// tDPPP.cc: test program for DPPP
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
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

#include "../../DP3.h"

using casacore::ArrayColumn;
using casacore::Block;
using casacore::IPosition;
using casacore::Matrix;
using casacore::near;
using casacore::ScalarColumn;
using casacore::Slicer;
using casacore::Table;
using casacore::TableIterator;
using casacore::Vector;

using Complex = std::complex<float>;

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_SUITE(dppp)

// This test program uses the MS in tNDPPP.in_MS.tgz.
// The MS contains 4 corr, 16 freq, 6 baselines, 18 time slots of 30 sec.
// Two time slots are missing between time slot 2 and 3.

void checkCopy(const std::string& in, const std::string& out, int nms) {
  Table tin(in);
  Table tout(out);
  BOOST_CHECK_EQUAL(tout.nrow(), size_t{6 * 24});
  for (int j = 0; j < nms; ++j) {
    // A few dummy time slots were inserted, so ignore those.
    Table t1 = tableCommand(
        "using style python "
        "select from $1 where rownumber() not in [0:6, 4*6:6*6, 21*6:24*6]",
        tout);
    ArrayColumn<Complex> data(t1, "DATA");
    ArrayColumn<bool> flag(t1, "FLAG");
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    IPosition dshape(2, 4, 16);
    IPosition dst(2, 0, j * 16);
    Slicer dslicer(dst, dshape);
    BOOST_CHECK_EQUAL(data(0).shape(), IPosition(2, 4, 16 * nms));
    BOOST_CHECK_EQUAL(flag(0).shape(), IPosition(2, 4, 16 * nms));
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 16 * nms / 8, 1));
    BOOST_CHECK(allEQ(data.getColumn(dslicer),
                      ArrayColumn<Complex>(tin, "DATA").getColumn()));
    BOOST_CHECK(allEQ(flag.getColumn(), false));
    BOOST_CHECK(allEQ(oflag.getColumn(), (unsigned char)(0)));
    BOOST_CHECK(
        allEQ(ArrayColumn<float>(t1, "WEIGHT_SPECTRUM").getColumn(), float(1)));
    //    cout<<ArrayColumn<double>(t1,"UVW").getColumn()<<
    //                  ArrayColumn<double>(tin,"UVW").getColumn();
    //    BOOST_CHECK (allEQ(ArrayColumn<double>(t1,"UVW").getColumn(),
    //                  ArrayColumn<double>(tin,"UVW").getColumn()));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "TIME").getColumn(),
                      ScalarColumn<double>(tin, "TIME").getColumn()));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "TIME_CENTROID").getColumn(),
                      ScalarColumn<double>(tin, "TIME_CENTROID").getColumn()));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "INTERVAL").getColumn(),
                      ScalarColumn<double>(tin, "INTERVAL").getColumn()));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "EXPOSURE").getColumn(),
                      ScalarColumn<double>(tin, "EXPOSURE").getColumn()));
  }
  {
    // Check the inserted time slots.
    // The MS misses a few time slots (3 and 4).
    Table t1 = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [0:6, 4*6:6*6, 21*6:24*6]",
        tout);
    ArrayColumn<Complex> data(t1, "DATA");
    ArrayColumn<bool> flag(t1, "FLAG");
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL(data(0).shape(), IPosition(2, 4, 16 * nms));
    BOOST_CHECK_EQUAL(flag(0).shape(), IPosition(2, 4, 16 * nms));
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 16 * nms / 8, 1));
    BOOST_CHECK(allEQ(data.getColumn(), Complex()));
    BOOST_CHECK(allEQ(flag.getColumn(), true));
    BOOST_CHECK(allEQ(oflag.getColumn(), (unsigned char)(0xff)));
    BOOST_CHECK(
        allEQ(ArrayColumn<float>(t1, "WEIGHT_SPECTRUM").getColumn(), float(0)));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "INTERVAL").getColumn(), 30.));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "EXPOSURE").getColumn(), 30.));
    double time = ScalarColumn<double>(tin, "TIME")(0);
    for (unsigned int i = 0; i < 6; ++i) {
      double timec = time - 30;
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME")(i), timec));
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME_CENTROID")(i), timec));
    }
    time = ScalarColumn<double>(tin, "TIME")(2 * 6);
    for (unsigned int i = 6; i < 18; ++i) {
      double timec = time + (i / 6) * 30.;
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME")(i), timec));
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME_CENTROID")(i), timec));
    }
    time = ScalarColumn<double>(tin, "TIME")(17 * 6);
    for (unsigned int i = 18; i < 36; ++i) {
      double timec = time + (i / 6 - 2) * 30.;
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME")(i), timec));
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME_CENTROID")(i), timec));
    }
  }
  // Now check if the SPECTRAL_WINDOW table is fine.
  Table spwin(tin.keywordSet().asTable("SPECTRAL_WINDOW"));
  Table spwout(tout.keywordSet().asTable("SPECTRAL_WINDOW"));
  for (int j = 0; j < nms; ++j) {
    IPosition dshape(1, 16);
    IPosition dst(1, j * 16);
    Slicer dslicer(dst, dshape);
    BOOST_CHECK(
        allEQ(ArrayColumn<double>(spwin, "CHAN_FREQ").getColumn(),
              ArrayColumn<double>(spwout, "CHAN_FREQ").getColumn(dslicer)));
    BOOST_CHECK(
        allEQ(ArrayColumn<double>(spwin, "CHAN_WIDTH").getColumn(),
              ArrayColumn<double>(spwout, "CHAN_WIDTH").getColumn(dslicer)));
    BOOST_CHECK(
        allEQ(ArrayColumn<double>(spwin, "EFFECTIVE_BW").getColumn(),
              ArrayColumn<double>(spwout, "EFFECTIVE_BW").getColumn(dslicer)));
    BOOST_CHECK(
        allEQ(ArrayColumn<double>(spwin, "RESOLUTION").getColumn(),
              ArrayColumn<double>(spwout, "RESOLUTION").getColumn(dslicer)));
  }
  BOOST_CHECK(allEQ(
      double(nms) * ScalarColumn<double>(spwin, "TOTAL_BANDWIDTH").getColumn(),
      ScalarColumn<double>(spwout, "TOTAL_BANDWIDTH").getColumn()));
  if (nms == 1) {
    BOOST_CHECK(
        allEQ(ScalarColumn<double>(spwin, "REF_FREQUENCY").getColumn(),
              ScalarColumn<double>(spwout, "REF_FREQUENCY").getColumn()));
  }
  // Check the TIME_RANGE in the OBSERVATION table.
  Table obsout(tout.keywordSet().asTable("OBSERVATION"));
  casacore::Vector<double> timeRange(
      ArrayColumn<double>(obsout, "TIME_RANGE").getColumn());
  BOOST_CHECK(near(timeRange(0), ScalarColumn<double>(tout, "TIME")(0) - 15));
  BOOST_CHECK(near(timeRange(1), ScalarColumn<double>(tout, "TIME")(143) + 15));
}

void checkCopyColumn(const std::string& in) {
  Table tin(in);
  BOOST_CHECK_EQUAL(tin.nrow(), size_t{6 * 24});
  ArrayColumn<Complex> data1(tin, "DATA");
  ArrayColumn<Complex> data2(tin, "COPY_DATA");
  ArrayColumn<float> weight1(tin, "NEW_WEIGHT_SPECTRUM");
  ArrayColumn<float> weight2(tin, "COPY_NEW_WEIGHT_SPECTRUM");
  BOOST_CHECK(allEQ(data1.getColumn(), data2.getColumn()));
  BOOST_CHECK(allEQ(weight1.getColumn(), weight2.getColumn()));
}

void testCopy() {
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    // Give starttime 35 sec before actual, hence 1 missing timeslot.
    ostr << "msin.starttime=03-Aug-2000/13:21:45\n";
    // Give endtime 90 sec after actual, hence 3 missing timeslots.
    ostr << "msin.endtime=03-Aug-2000/13:33:15\n";
    ostr << "msout=tNDPPP_tmp.MS1\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  checkCopy("tNDPPP_tmp.MS", "tNDPPP_tmp.MS1", 1);
}

void testCopyColumn() {
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS1\n";
    ostr << "msout=.\n";
    ostr << "msout.datacolumn=COPY_DATA\n";
    ostr << "msout.weightcolumn=NEW_WEIGHT_SPECTRUM\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");

  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS1\n";
    ostr << "msin.datacolumn=COPY_DATA\n";
    ostr << "msin.weightcolumn=NEW_WEIGHT_SPECTRUM\n";
    ostr << "msout=.\n";
    ostr << "msout.datacolumn=DATA\n";
    ostr << "msout.weightcolumn=COPY_NEW_WEIGHT_SPECTRUM\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");

  checkCopyColumn("tNDPPP_tmp.MS1");
}

void testMultiIn() {
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=[tNDPPP_tmp.MS1, tNDPPP_tmp.MS1]\n";
    ostr << "msout=tNDPPP_tmp.MS1a\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  checkCopy("tNDPPP_tmp.MS", "tNDPPP_tmp.MS1a", 2);
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=[tNDPPP_tmp.MS1, tNDPPP_tmp.MS1]\n";
    ostr << "msin.datacolumn=CORRECTED_DATA\n";
    ostr << "msin.weightcolumn=NEW_WEIGHT_SPECTRUM\n";
    ostr << "msin.missingdata=true\n";
    ostr << "msin.baseline=0,2&6\n";
    ostr << "msout=tNDPPP_tmp.MS1a\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  Table tab("tNDPPP_tmp.MS1a");
  BOOST_CHECK_EQUAL(tab.nrow(), size_t{48});
  BOOST_CHECK(allEQ(ArrayColumn<Complex>(tab, "DATA").getColumn(), Complex()));
  BOOST_CHECK(tab.tableDesc().isColumn("WEIGHT_SPECTRUM"));
  BOOST_CHECK(allEQ(ArrayColumn<bool>(tab, "FLAG").getColumn(), true));
  BOOST_CHECK(allEQ(ScalarColumn<int>(tab, "ANTENNA2").getColumn(), 6));
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=[notexist, tNDPPP_tmp.MS1, notexist, notexist]\n";
    ostr << "msin.datacolumn=CORRECTED_DATA\n";
    ostr << "msin.missingdata=true\n";
    ostr << "msin.orderms=false\n";
    ostr << "msin.baseline=0,2&6\n";
    ostr << "msout=tNDPPP_tmp.MS1b\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  tab = Table("tNDPPP_tmp.MS1b");
  BOOST_CHECK_EQUAL(tab.nrow(), size_t{48});
  BOOST_CHECK(allEQ(ArrayColumn<Complex>(tab, "DATA").getColumn(), Complex()));
  BOOST_CHECK(allEQ(ArrayColumn<bool>(tab, "FLAG").getColumn(), true));
  BOOST_CHECK(allEQ(ScalarColumn<int>(tab, "ANTENNA2").getColumn(), 6));
}

void checkAvg(const std::string& outName) {
  Table tin("tNDPPP_tmp.MS");
  Table tout(outName);
  ScalarColumn<double> timeCol(tin, "TIME");
  double time = 0.5 * (timeCol(tin.nrow() - 1) + timeCol(0));
  BOOST_CHECK_EQUAL(tout.nrow(), size_t{6});
  Block<casacore::String> colNames(2);
  colNames[0] = "ANTENNA1";
  colNames[1] = "ANTENNA2";
  TableIterator iterin(tin, colNames);
  TableIterator iterout(tout, colNames);
  // Iterate over baseline to be able to average the input in an easy way.
  while (!iterin.pastEnd()) {
    Table t2(iterin.table());
    Table t1(iterout.table());
    ArrayColumn<Complex> data(t1, "DATA");
    ArrayColumn<bool> flag(t1, "FLAG");
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL(data(0).shape(), IPosition(2, 4, 1));
    BOOST_CHECK_EQUAL(flag(0).shape(), IPosition(2, 4, 1));
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 16 / 8, 20));
    // Average the original data over all channels and times.
    casacore::Array<Complex> dataAvg = partialMeans(
        ArrayColumn<Complex>(t2, "DATA").getColumn(), IPosition(2, 1, 2));
    casacore::Array<Complex> dataRes = dataAvg.reform(IPosition(3, 4, 1, 1));
    BOOST_CHECK(allNear(data.getColumn(), dataRes, 1e-5));
    BOOST_CHECK(allNear(ArrayColumn<float>(t1, "WEIGHT_SPECTRUM").getColumn(),
                        float(18 * 16), 1e-5));
    BOOST_CHECK(allEQ(flag.getColumn(), false));
    BOOST_CHECK(
        allNear(ScalarColumn<double>(t1, "TIME").getColumn(), time, 1e-5));
    BOOST_CHECK(allNear(ScalarColumn<double>(t1, "TIME_CENTROID").getColumn(),
                        time, 1e-5));
    BOOST_CHECK(
        allEQ(ScalarColumn<double>(t1, "INTERVAL").getColumn(), 20 * 30.));
    BOOST_CHECK(
        allEQ(ScalarColumn<double>(t1, "EXPOSURE").getColumn(), 20 * 30.));
    // Two time entries should be flagged.
    casacore::Array<unsigned char> of = oflag.getColumn();
    BOOST_CHECK(allEQ(of(Slicer(IPosition(3, 0, 0, 0), IPosition(3, 2, 3, 1))),
                      (unsigned char)(0)));
    BOOST_CHECK(allEQ(of(Slicer(IPosition(3, 0, 3, 0), IPosition(3, 2, 2, 1))),
                      (unsigned char)(0xff)));
    BOOST_CHECK(allEQ(of(Slicer(IPosition(3, 0, 5, 0), IPosition(3, 2, 15, 1))),
                      (unsigned char)(0)));
    iterin.next();
    iterout.next();
  }
  // Now check if the SPECTRAL_WINDOW table is fine.
  Table spwin(tin.keywordSet().asTable("SPECTRAL_WINDOW"));
  Table spwout(tout.keywordSet().asTable("SPECTRAL_WINDOW"));
  Matrix<double> cw = ArrayColumn<double>(spwout, "CHAN_WIDTH").getColumn();
  BOOST_CHECK_EQUAL(cw.size(), size_t{1});
  BOOST_CHECK(near(cw(0, 0),
                   sum(ArrayColumn<double>(spwin, "CHAN_WIDTH").getColumn())));
  Matrix<double> cfi = ArrayColumn<double>(spwin, "CHAN_FREQ").getColumn();
  Matrix<double> cfo = ArrayColumn<double>(spwout, "CHAN_FREQ").getColumn();
  BOOST_CHECK(near(cfo(0, 0), 0.5 * (cfi(0, 0) + cfi(15, 0))));
  Matrix<double> ce = ArrayColumn<double>(spwout, "EFFECTIVE_BW").getColumn();
  BOOST_CHECK_EQUAL(ce.size(), size_t{1});
  BOOST_CHECK(near(ce(0, 0), cw(0, 0)));
  Matrix<double> cr = ArrayColumn<double>(spwout, "RESOLUTION").getColumn();
  BOOST_CHECK_EQUAL(cr.size(), size_t{1});
  BOOST_CHECK(near(cr(0, 0), cw(0, 0)));
  BOOST_CHECK(
      near(ScalarColumn<double>(spwin, "TOTAL_BANDWIDTH")(0), cw(0, 0)));
  BOOST_CHECK(allEQ(ScalarColumn<double>(spwin, "REF_FREQUENCY").getColumn(),
                    ScalarColumn<double>(spwout, "REF_FREQUENCY").getColumn()));
  // Check the TIME_RANGE in the OBSERVATION table.
  Table obsout(tout.keywordSet().asTable("OBSERVATION"));
  Vector<double> timeRange(
      ArrayColumn<double>(obsout, "TIME_RANGE").getColumn());
  BOOST_CHECK(near(timeRange(0), ScalarColumn<double>(tout, "TIME")(0) - 300));
  BOOST_CHECK(near(timeRange(1), ScalarColumn<double>(tout, "TIME")(0) + 295));
}

void testAvg1() {
  {
    // Average in a single step.
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin.name=tNDPPP_tmp.MS\n";
    // Give start and end time as actual, hence no missing timeslots.
    ostr << "msin.starttime=03-Aug-2000/13:22:20\n";
    ostr << "msin.endtime=03-Aug-2000/13:31:45\n";
    ostr << "msout.name=tNDPPP_tmp.MS2\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[avg,count]\n";
    ostr << "avg.type=average\n";
    ostr << "avg.timestep=20\n";
    ostr << "avg.freqstep=100\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  checkAvg("tNDPPP_tmp.MS2");
}

void testAvg2() {
  // Averaging in multiple steps should be the same as above.
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=tNDPPP_tmp.MS3\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[avg1,avg2,avg3,avg4]\n";
    ostr << "avg1.type=average\n";
    ostr << "avg1.timestep=5\n";
    ostr << "avg1.freqstep=2\n";
    ostr << "avg2.type=average\n";
    ostr << "avg2.timestep=1\n";
    ostr << "avg2.freqstep=2\n";
    ostr << "avg3.type=average\n";
    ostr << "avg3.timestep=2\n";
    ostr << "avg3.freqstep=1\n";
    ostr << "avg4.type=average\n";
    ostr << "avg4.timestep=2\n";
    ostr << "avg4.freqstep=4\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  checkAvg("tNDPPP_tmp.MS3");
}

void testAvg3() {
  // Averaging in multiple steps with multiple outputs should be the same
  // as above.
  {
    std::ofstream ostr1("tNDPPP_tmp.parset1");
    std::ofstream ostr2("tNDPPP_tmp.parset2");
    std::ofstream ostr3("tNDPPP_tmp.parset3");
    std::ofstream ostr4("tNDPPP_tmp.parset4");
    ostr1 << "msin=tNDPPP_tmp.MS\n";
    ostr1 << "msout=tNDPPP_tmp.MS4a\n";
    ostr1 << "msout.overwrite=true\n";
    ostr1 << "steps=[avg1]\n";
    ostr1 << "avg1.type=average\n";
    ostr1 << "avg1.timestep=5\n";
    ostr1 << "avg1.freqstep=2\n";
    ostr2 << "msin=tNDPPP_tmp.MS4a\n";
    ostr2 << "msout=tNDPPP_tmp.MS4b\n";
    ostr2 << "msout.overwrite=true\n";
    ostr2 << "steps=[avg2]\n";
    ostr2 << "avg2.type=average\n";
    ostr2 << "avg2.timestep=1\n";
    ostr2 << "avg2.freqstep=2\n";
    ostr3 << "msin=tNDPPP_tmp.MS4b\n";
    ostr3 << "msout=tNDPPP_tmp.MS4c\n";
    ostr3 << "msout.overwrite=true\n";
    ostr3 << "steps=[avg3]\n";
    ostr3 << "avg3.type=average\n";
    ostr3 << "avg3.timestep=2\n";
    ostr3 << "avg3.freqstep=1\n";
    ostr4 << "msin=tNDPPP_tmp.MS4c\n";
    ostr4 << "msout=tNDPPP_tmp.MS4d\n";
    ostr4 << "msout.overwrite=true\n";
    ostr4 << "steps=[avg4]\n";
    ostr4 << "avg4.type=average\n";
    ostr4 << "avg4.timestep=2\n";
    ostr4 << "avg4.freqstep=4\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset1");
  dp3::base::DP3::execute("tNDPPP_tmp.parset2");
  dp3::base::DP3::execute("tNDPPP_tmp.parset3");
  dp3::base::DP3::execute("tNDPPP_tmp.parset4");
  checkAvg("tNDPPP_tmp.MS4d");
}

// This function tests if the correct start time is used when selecting times.
void testAvg4() {
  {
    // Average in a single step.
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    // Give start a few seconds after first one, hence skip first time slot.
    ostr << "msin.starttime=03-Aug-2000/13:22:25\n";
    ostr << "msout=tNDPPP_tmp.MS5\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[avg]\n";
    ostr << "avg.type=average\n";
    ostr << "avg.timestep=2\n";
    ostr << "avg.freqstep=100\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  {
    // Only check the times;
    Table tab("tNDPPP_tmp.MS");
    // First time to be used.
    double time = ScalarColumn<double>(tab, "TIME")(6) + 15;
    Table t2("tNDPPP_tmp.MS5");
    BOOST_CHECK_EQUAL(t2.nrow(), size_t{6 * 10});
    ScalarColumn<double> timeCol(t2, "TIME");
    for (unsigned int i = 0; i < t2.nrow(); ++i) {
      BOOST_CHECK(near(timeCol(i), time));
      if (i % 6 == 5) time += 60;
    }
  }
}

void testUpdate1() {
  // Test if update works fine.
  // In fact, it does not do anything apart from rewriting the current flags.
  // However, it should ignore the inserted time slots.
  casacore::Array<bool> flags;
  {
    Table tab("tNDPPP_tmp.MS");
    flags = ArrayColumn<bool>(tab, "FLAG").getColumn();
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=''\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  // Check that the flags did not change.
  {
    Table tab("tNDPPP_tmp.MS");
    BOOST_CHECK(allEQ(ArrayColumn<bool>(tab, "FLAG").getColumn(), flags));
  }
}

void testUpdate2() {
  // Test if update all flags works fine.
  {
    Table tab("tNDPPP_tmp.MS");
    tab.deepCopy("tNDPPP_tmp.MS_copy1", Table::New);
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS_copy1\n";
    ostr << "msout=.\n";
    ostr << "steps=[preflag]\n";
    ostr << "preflag.blmin=1e6\n";  // should flag all data
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  // Check that all flags are true.
  {
    Table tab("tNDPPP_tmp.MS_copy1");
    BOOST_CHECK(allEQ(ArrayColumn<bool>(tab, "FLAG").getColumn(), true));
  }
}

void testUpdateScale() {
  // Test if update data works fine.
  casacore::Array<Complex> data;
  casacore::Array<bool> flags;
  {
    Table tab("tNDPPP_tmp.MS");
    data = ArrayColumn<Complex>(tab, "DATA").getColumn();
    flags = ArrayColumn<bool>(tab, "FLAG").getColumn();
    tab.deepCopy("tNDPPP_tmp.MS_copy1", Table::New);
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS_copy1\n";
    ostr << "msout=tNDPPP_tmp.MS_copy1\n";  // same name means update
    ostr << "steps=[scaledata]\n";
    ostr << "scaledata.coeffs=2\n";
    ostr << "scaledata.stations=*\n";
    ostr << "scaledata.scalesize=false\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  // Check that all data is doubled.
  {
    Table tab("tNDPPP_tmp.MS_copy1");
    data *= Complex(2, 0);
    BOOST_CHECK(
        allNear(ArrayColumn<Complex>(tab, "DATA").getColumn(), data, 1e-5));
    BOOST_CHECK(allEQ(ArrayColumn<bool>(tab, "FLAG").getColumn(), flags));
  }
}

void checkFlags(const std::string& outName) {
  // Only check the FULL_RES_FLAGS and table size.
  // The flags are created in various ways, but should be the same in all cases.
  // 3 time slots are averaged to 1.
  Table tout(outName);
  BOOST_CHECK_EQUAL(tout.nrow(), size_t{6 * 4});
  // Check the full-res-flags.
  // Channels 0,2,6,7,8 are flagged everywhere.
  {
    // Input time slots 3,4 are inserted, thus flagged.
    Table t1 = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [0:6]",
        tout);
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 2, 6));
    casacore::Array<unsigned char> flags = oflag.getColumn();
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 0, 0), IPosition(3, 0, 2, 5)),
                      (unsigned char)(0xc5)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 0, 0), IPosition(3, 1, 2, 5)),
                      (unsigned char)(0x01)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 3, 0), IPosition(3, 0, 4, 5)),
                      (unsigned char)(0xff)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 3, 0), IPosition(3, 1, 4, 5)),
                      (unsigned char)(0x0f)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 5, 0), IPosition(3, 0, 5, 5)),
                      (unsigned char)(0xc5)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 5, 0), IPosition(3, 1, 5, 5)),
                      (unsigned char)(0x01)));
  }
  {
    // Input time slots 10,11 are flagged.
    Table t1 = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [6:12]",
        tout);
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 2, 6));
    casacore::Array<unsigned char> flags = oflag.getColumn();
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 0, 0), IPosition(3, 0, 3, 5)),
                      (unsigned char)(0xc5)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 0, 0), IPosition(3, 1, 3, 5)),
                      (unsigned char)(0x01)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 4, 0), IPosition(3, 0, 5, 5)),
                      (unsigned char)(0xff)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 4, 0), IPosition(3, 1, 5, 5)),
                      (unsigned char)(0x0f)));
  }
  {
    // Input time slots 12-17 are not flagged.
    Table t1 = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [12:18]",
        tout);
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 2, 6));
    casacore::Array<unsigned char> flags = oflag.getColumn();
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 0, 0), IPosition(3, 0, 5, 5)),
                      (unsigned char)(0xc5)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 0, 0), IPosition(3, 1, 5, 5)),
                      (unsigned char)(0x01)));
  }
  {
    // Input time slot 20-23 did not exist, thus flagged in average.
    Table t1 = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [18:24]",
        tout);
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 2, 6));
    casacore::Array<unsigned char> flags = oflag.getColumn();
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 0, 0), IPosition(3, 0, 1, 5)),
                      (unsigned char)(0xc5)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 0, 0), IPosition(3, 1, 1, 5)),
                      (unsigned char)(0x01)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 2, 0), IPosition(3, 0, 5, 5)),
                      (unsigned char)(0xff)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 2, 0), IPosition(3, 1, 5, 5)),
                      (unsigned char)(0x0f)));
  }
}

void testFlags1() {
  {
    // Most simple case.
    // Just flag some channels and time stamps.
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msin.startchan=1\n";
    ostr << "msin.nchan=12\n";
    ostr << "msout=tNDPPP_tmp.MS5\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[preflag,average]\n";
    ostr << "preflag.expr='flag1 or flag2'\n";
    ostr << "preflag.flag1.timeslot=[10,11]\n";
    ostr << "preflag.flag2.chan=[0,2,6..8]\n";
    ostr << "average.timestep=6\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  checkFlags("tNDPPP_tmp.MS5");
}

void testFlags2() {
  {
    // A more advanced case in two NDPPP steps.
    // Flag the channels, but shifted 2. In the next step the first 2 channels
    // are skipped, thus the channel numbers shift back 2.
    // An averaged time slot is flagged, so the original time slots should
    // come out flagged.
    std::ofstream ostr1("tNDPPP_tmp.parset1");
    ostr1 << "msin=tNDPPP_tmp.MS\n";
    ostr1 << "msin.nchan=15\n";
    ostr1 << "msout=tNDPPP_tmp.MS6a\n";
    ostr1 << "msout.overwrite=true\n";
    ostr1 << "steps=[preflag,average]\n";
    ostr1 << "preflag.chan=[2,4,8..10]\n";
    ostr1 << "average.timestep=2\n";
    std::ofstream ostr2("tNDPPP_tmp.parset2");
    ostr2 << "msin=tNDPPP_tmp.MS6a\n";
    ostr2 << "msin.startchan=2*1\n";  // output chan 0,2 are now flagged
    ostr2 << "msin.nchan=nchan-3\n";
    ostr2 << "msout=tNDPPP_tmp.MS6b\n";
    ostr2 << "msout.overwrite=true\n";
    ostr2 << "steps=[preflag,average]\n";
    ostr2 << "preflag.timeslot=5\n";  // is 10,11 in input
    ostr2 << "average.timestep=3\n";
    ostr2 << "average.freqstep=2\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset1");
  dp3::base::DP3::execute("tNDPPP_tmp.parset2");
  checkFlags("tNDPPP_tmp.MS6b");
}

void testFlags3() {
  {
    // Even a bit more advanced, also in two NDPPP runs.
    // Input channels 6,7,8 are flagged by flagging their averaged channel.
    // This is done at the end of run 1, so the averager of run 2 should pick
    // up those flags.
    std::ofstream ostr1("tNDPPP_tmp.parset1");
    ostr1 << "msin=tNDPPP_tmp.MS\n";
    ostr1 << "msin.nchan=15\n";
    ostr1 << "msout=tNDPPP_tmp.MS7a\n";
    ostr1 << "msout.overwrite=true\n";
    ostr1 << "steps=[preflag,average,pre2]\n";
    ostr1 << "preflag.chan=[0,2]\n";
    ostr1 << "average.timestep=2\n";
    ostr1 << "average.freqstep=3\n";
    ostr1 << "pre2.type=preflag\n";
    ostr1 << "pre2.chan=2\n";  // is input channel 6,7,8
    std::ofstream ostr2("tNDPPP_tmp.parset2");
    ostr2 << "msin=tNDPPP_tmp.MS7a\n";
    ostr2 << "msin.nchan=4\n";
    ostr2 << "msout=tNDPPP_tmp.MS7b\n";
    ostr2 << "msout.overwrite=true\n";
    ostr2 << "steps=[preflag,average]\n";
    ostr2 << "preflag.timeslot=5\n";  // is 10,11 in input
    ostr2 << "average.timestep=3\n";
    ostr2 << "average.freqstep=2\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset1");
  dp3::base::DP3::execute("tNDPPP_tmp.parset2");
  checkFlags("tNDPPP_tmp.MS7b");
}

void testStationAdd() {
  // Add station RT0, 1 and 2.
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=tNDPPP_tmp.MSa\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[stationadd]\n";
    ostr << "stationadd.stations={RTnew:[RT0..2]}\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  Table t1("tNDPPP_tmp.MS/ANTENNA");
  Table t2("tNDPPP_tmp.MSa/ANTENNA");
  BOOST_CHECK_EQUAL(t2.nrow(), t1.nrow() + 1);  // 1 antenna has been added
  BOOST_CHECK_EQUAL(ScalarColumn<casacore::String>(t2, "NAME")(t2.nrow() - 1),
                    "RTnew");
  int oldNant = t1.nrow();
  t1 = Table("tNDPPP_tmp.MS/FEED");
  t2 = Table("tNDPPP_tmp.MSa/FEED");
  BOOST_CHECK_EQUAL(t2.nrow(), t1.nrow() + 1);  // 1 antenna has been added
  BOOST_CHECK_EQUAL(ScalarColumn<int>(t2, "ANTENNA_ID")(t2.nrow() - 1),
                    oldNant);
  t1 = Table("tNDPPP_tmp.MS");
  t2 = Table("tNDPPP_tmp.MSa");
  BOOST_CHECK_EQUAL(t2.nrow(),
                    t1.nrow() + 40 + 12);  // 2 baselines and 2 time slots added
}

void testFilter1() {
  // Remove all baselines containing station RT1 or 6.
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=tNDPPP_tmp.MSa\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[filter]\n";
    ostr << "filter.baseline=!RT[16]&&*\n";
    ostr << "filter.remove=true\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  Table t1("tNDPPP_tmp.MS/ANTENNA");
  Table t2("tNDPPP_tmp.MSa/ANTENNA");
  // Note: the ANTENNA table also contained RT8, RT9, etc., but they do not
  // have baselines. So these were removed as well meaning only 0,2,7 are left.
  Vector<unsigned> rownrs(3);
  rownrs[0] = 0;
  rownrs[1] = 2;
  rownrs[2] = 7;
  Table t1s = t1(rownrs);
  BOOST_CHECK_EQUAL(t2.nrow(), t1s.nrow());
  BOOST_CHECK(allEQ(ScalarColumn<casacore::String>(t2, "NAME").getColumn(),
                    ScalarColumn<casacore::String>(t1s, "NAME").getColumn()));
  t1 = Table("tNDPPP_tmp.MS/FEED");
  t2 = Table("tNDPPP_tmp.MSa/FEED");
  t1s = t1(rownrs);
  BOOST_CHECK_EQUAL(t2.nrow(), t1s.nrow());
  // The ANTENNA_IDs in the FEED table must be 0,1,2.
  Vector<int> ids(t2.nrow());
  indgen(ids);
  BOOST_CHECK(allEQ(ScalarColumn<int>(t2, "ANTENNA_ID").getColumn(), ids));
  // Check the main table.
  t1 = Table("tNDPPP_tmp.MS");
  t2 = Table("tNDPPP_tmp.MSa");
  BOOST_CHECK_EQUAL(
      t2.nrow(), t1.nrow() - 72 + 4);  // 4 baselines removed, 2 timeslots added
  t1s = t1((t1.col("ANTENNA1") == 0 || t1.col("ANTENNA1") == 2) &&
           t1.col("ANTENNA2") == 7);
  // A few dummy time slots were inserted, so ignore those.
  Table t2s = t2(t2.nodeRownr() < 6 || t2.nodeRownr() >= 10);
  BOOST_CHECK(allEQ(ArrayColumn<Complex>(t2s, "DATA").getColumn(),
                    ArrayColumn<Complex>(t1s, "DATA").getColumn()));
  t2s = t2(t2.nodeRownr() % 2 == 0);
  BOOST_CHECK(allEQ(ScalarColumn<int>(t2s, "ANTENNA1").getColumn(), 0));
  BOOST_CHECK(allEQ(ScalarColumn<int>(t2s, "ANTENNA2").getColumn(), 2));
  t2s = t2(t2.nodeRownr() % 2 == 1);
  BOOST_CHECK(allEQ(ScalarColumn<int>(t2s, "ANTENNA1").getColumn(), 1));
  BOOST_CHECK(allEQ(ScalarColumn<int>(t2s, "ANTENNA2").getColumn(), 2));
}

void testFilter2() {
  // Keep all baselines.
  // First by not specifying baseline selection, second by all baselines.
  // Also alter between remove and !remove.
  for (int iter = 0; iter < 4; ++iter) {
    {
      std::ofstream ostr("tNDPPP_tmp.parset");
      ostr << "msin=tNDPPP_tmp.MS\n";
      ostr << "msout=tNDPPP_tmp.MSa\n";
      ostr << "msout.overwrite=true\n";
      ostr << "steps=[filter]\n";
      if (iter % 2 == 1) {
        ostr << "filter.baseline=*&&*\n";
      }
      if (iter / 2 == 1) {
        ostr << "filter.remove=true\n";
      }
    }
    dp3::base::DP3::execute("tNDPPP_tmp.parset");
    // cout << "check ANTENNA"<<'\n';
    Table t1("tNDPPP_tmp.MS/ANTENNA");
    Table t2("tNDPPP_tmp.MSa/ANTENNA");
    // Note: the ANTENNA table also contained RT8, RT9, etc., but they do not
    // have baselines. So these were removed meaning only 0,1,2,6,7 are left.
    Vector<unsigned> rownrs(5);
    rownrs[0] = 0;
    rownrs[1] = 1;
    rownrs[2] = 2;
    rownrs[3] = 6;
    rownrs[4] = 7;
    Table t1s(t1);
    if (iter / 2 == 1) {
      t1s = t1(rownrs);
    }
    BOOST_CHECK_EQUAL(t2.nrow(), t1s.nrow());
    BOOST_CHECK(allEQ(ScalarColumn<casacore::String>(t2, "NAME").getColumn(),
                      ScalarColumn<casacore::String>(t1s, "NAME").getColumn()));
    // cout << "check FEED"<<'\n';
    t1 = Table("tNDPPP_tmp.MS/FEED");
    t2 = Table("tNDPPP_tmp.MSa/FEED");
    t1s = t1;
    if (iter / 2 == 1) {
      t1s = t1(rownrs);
    }
    BOOST_CHECK_EQUAL(t2.nrow(), t1s.nrow());
    // The ANTENNA_IDs in the FEED table must be 0,1,2.
    Vector<int> ids(t2.nrow());
    indgen(ids);
    BOOST_CHECK(allEQ(ScalarColumn<int>(t2, "ANTENNA_ID").getColumn(), ids));
    // Check the main table.
    t1 = Table("tNDPPP_tmp.MS");
    t2 = Table("tNDPPP_tmp.MSa");
    BOOST_CHECK_EQUAL(t2.nrow(), t1.nrow() + 12);  // 2 timeslots added
    // A few dummy time slots were inserted, so ignore those.
    Table t2s = t2(t2.nodeRownr() < 18 || t2.nodeRownr() >= 30);
    BOOST_CHECK(allEQ(ArrayColumn<Complex>(t2s, "DATA").getColumn(),
                      ArrayColumn<Complex>(t1, "DATA").getColumn()));
    int ant1[] = {0, 0, 1, 1, 2, 2};
    int ant2[] = {6, 7, 6, 7, 6, 7};
    int sub = (iter / 2 == 0 ? 0 : 3);  // if remove, ant2 6->3 and 7->4
    for (int i = 0; i < 6; ++i) {
      t2s = t2(t2.nodeRownr() % 6 == i);
      BOOST_CHECK(
          allEQ(ScalarColumn<int>(t2s, "ANTENNA1").getColumn(), ant1[i]));
      BOOST_CHECK(
          allEQ(ScalarColumn<int>(t2s, "ANTENNA2").getColumn(), ant2[i] - sub));
    }
  }
}

void testFilter3() {
  // Remove some baselines, update original file with different data column
  // This test justs tests if it runs without throwing exceptions
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=.\n";
    ostr << "msout.datacolumn=DATA_FILTER\n";
    ostr << "steps=[filter]\n";
    ostr << "filter.baseline=!RT[16]&&*\n";
    ostr << "filter.remove=False\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
}

void testClear() {
  casacore::Array<bool> flags;
  // First flag in the same way as testFlags1.
  testFlags1();
  // Get the resulting flags.
  {
    Table tab("tNDPPP_tmp.MS5");
    flags.reference(ArrayColumn<bool>(tab, "FLAG").getColumn());
  }
  // Flag all data.
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS5\n";
    ostr << "msout=tNDPPP_tmp.MS5a\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[preflag]\n";
    ostr << "preflag.baseline=[[*]]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  checkFlags("tNDPPP_tmp.MS5a");
  {
    Table tab("tNDPPP_tmp.MS5a");
    BOOST_CHECK(allEQ(ArrayColumn<bool>(tab, "FLAG").getColumn(), true));
  }
  // Clear the flags.
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS5a\n";
    ostr << "msout=tNDPPP_tmp.MS5b\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[preflag]\n";
    ostr << "preflag.mode=clear\n";
    ostr << "preflag.baseline=[[*]]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  checkFlags("tNDPPP_tmp.MS5b");
  {
    Table tab("tNDPPP_tmp.MS5b");
    BOOST_CHECK(allEQ(ArrayColumn<bool>(tab, "FLAG").getColumn(), false));
  }
}

void testMultiOut() {
  {
    // First make the reference output MS.
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=tNDPPP_tmp.MS_copy\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  // Test if update data works fine with multiple outputs:
  // read from tNDPPP_tmp.MS, write to copy3, update to copy3
  casacore::Array<Complex> data;
  casacore::Array<bool> flags;
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "steps=[scaledata,out1,scaledata,out2]\n";
    ostr << "scaledata.coeffs=2\n";
    ostr << "scaledata.stations=*\n";
    ostr << "scaledata.scalesize=false\n";
    ostr << "out1.type=out\n";
    ostr << "out1.name=tNDPPP_tmp.MS_copy3\n";
    ostr << "out1.overwrite=true\n";
    ostr << "out2.type=out\n";
    ostr << "out2.name=.\n";  // Defaults to the previous out, so _copy3
    ostr << "out2.datacolumn=DATA_2\n";
    ostr << "msout=tNDPPP_tmp.MS_copy4\n";  // same name means update
    ostr << "msout.overwrite=true\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset");
  // Check that tables exist, contain the specified columns
  {
    Table tab1("tNDPPP_tmp.MS_copy");
    Table tab2("tNDPPP_tmp.MS_copy3");
    /// cout<<ArrayColumn<Complex>(tab2,"DATA_2").getColumn();
    /// cout<<Complex(2,0)*ArrayColumn<Complex>(tab1,"DATA").getColumn();
    BOOST_CHECK(allNear(
        ArrayColumn<Complex>(tab2, "DATA").getColumn(),
        Complex(2, 0) * ArrayColumn<Complex>(tab1, "DATA").getColumn(), 1e-5));
    BOOST_CHECK(allNear(
        ArrayColumn<Complex>(tab2, "DATA_2").getColumn(),
        Complex(4, 0) * ArrayColumn<Complex>(tab1, "DATA").getColumn(), 1e-5));
    BOOST_CHECK(allNear(ArrayColumn<float>(tab2, "WEIGHT_SPECTRUM").getColumn(),
                        ArrayColumn<float>(tab1, "WEIGHT_SPECTRUM").getColumn(),
                        1e-5));
    BOOST_CHECK(allEQ(ArrayColumn<bool>(tab2, "FLAG").getColumn(),
                      ArrayColumn<bool>(tab1, "FLAG").getColumn()));
    Table tab3("tNDPPP_tmp.MS_copy4");
    BOOST_CHECK(allNear(
        ArrayColumn<Complex>(tab3, "DATA").getColumn(),
        Complex(4, 0) * ArrayColumn<Complex>(tab1, "DATA").getColumn(), 1e-5));
    BOOST_CHECK(allNear(ArrayColumn<float>(tab3, "WEIGHT_SPECTRUM").getColumn(),
                        ArrayColumn<float>(tab1, "WEIGHT_SPECTRUM").getColumn(),
                        1e-5));
    BOOST_CHECK(allEQ(ArrayColumn<bool>(tab3, "FLAG").getColumn(),
                      ArrayColumn<bool>(tab1, "FLAG").getColumn()));
  }
}

void tryErr(const std::string& parsetName) {
  bool err = false;
  try {
    dp3::base::DP3::execute(parsetName);
  } catch (const std::exception& x) {
    err = true;
  }
  BOOST_CHECK(err);
}

void testErrorOut() {
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "steps=[filter,out1,average,out2]\n";
    ostr << "out1.type=out\n";
    ostr << "out1.name=''\n";
    ostr << "out2.type=out\n";
    ostr << "out2.name=.\n";  // update not possible when avg
    ostr << "msout=''\n";
    tryErr("tNDPPP_tmp.parset");
  }
  {
    std::ofstream ostr("tNDPPP_tmp.parset");
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "steps=[average,out1,filter,out2]\n";
    ostr << "out1.type=out\n";
    ostr << "out1.name=tNDPPP_tmp.MSx\n";
    ostr << "out1.overwrite=true\n";
    ostr << "filter.remove=true\n";
    ostr << "out2.type=out\n";
    ostr << "out2.name=./tNDPPP_tmp.MSx"
         << '\n';  // update not possible (filter)
    ostr << "msout=''\n";
    tryErr("tNDPPP_tmp.parset");
  }
}

BOOST_AUTO_TEST_CASE(test_copy, *utf::disabled()) { testCopy(); }

BOOST_AUTO_TEST_CASE(test_copy_column, *utf::depends_on("dppp/test_copy")) {
  testCopyColumn();
}

BOOST_AUTO_TEST_CASE(test_multi_in, *utf::depends_on("dppp/test_copy_column")) {
  testMultiIn();
}

BOOST_AUTO_TEST_CASE(test_avg1, *utf::depends_on("dppp/test_multi_in")) {
  testAvg1();
}

BOOST_AUTO_TEST_CASE(test_avg2, *utf::depends_on("dppp/test_avg1")) {
  testAvg2();
}

BOOST_AUTO_TEST_CASE(test_avg3, *utf::depends_on("dppp/test_avg2")) {
  testAvg3();
}

BOOST_AUTO_TEST_CASE(test_avg4, *utf::depends_on("dppp/test_avg3")) {
  testAvg4();
}

BOOST_AUTO_TEST_CASE(test_update1, *utf::depends_on("dppp/test_avg4")) {
  testUpdate1();
}

BOOST_AUTO_TEST_CASE(test_update2, *utf::depends_on("dppp/test_update1")) {
  testUpdate2();
}

BOOST_AUTO_TEST_CASE(test_update_scale, *utf::depends_on("dppp/test_update2")) {
  testUpdateScale();
}

BOOST_AUTO_TEST_CASE(test_flags1, *utf::depends_on("dppp/test_update_scale")) {
  testFlags1();
}

BOOST_AUTO_TEST_CASE(test_flags2, *utf::depends_on("dppp/test_flags1")) {
  testFlags2();
}

BOOST_AUTO_TEST_CASE(test_flags3, *utf::depends_on("dppp/test_flags2")) {
  testFlags3();
}

BOOST_AUTO_TEST_CASE(test_station_add, *utf::depends_on("dppp/test_flags3")) {
  testStationAdd();
}

BOOST_AUTO_TEST_CASE(test_filter1, *utf::depends_on("dppp/test_station_add")) {
  testFilter1();
}

BOOST_AUTO_TEST_CASE(test_filter2, *utf::depends_on("dppp/test_filter1")) {
  testFilter2();
}

BOOST_AUTO_TEST_CASE(test_filter3, *utf::depends_on("dppp/test_filter2")) {
  testFilter3();
}

BOOST_AUTO_TEST_CASE(test_clear, *utf::depends_on("dppp/test_filter3")) {
  testClear();
}

BOOST_AUTO_TEST_CASE(test_multi_out, *utf::depends_on("dppp/test_clear")) {
  testMultiOut();
}

BOOST_AUTO_TEST_CASE(test_error_out, *utf::depends_on("dppp/test_multi_out")) {
  testErrorOut();
}

BOOST_AUTO_TEST_SUITE_END()
