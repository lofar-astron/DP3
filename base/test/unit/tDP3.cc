// tDP3.cc: test program for DPPP
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../DP3.h"

#include "../../../common/test/unit/fixtures/fDirectory.h"

#include <casacore/tables/Tables.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayPartMath.h>
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <stdexcept>

using casacore::ArrayColumn;
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

BOOST_AUTO_TEST_SUITE(dp3)

// This test program uses the MS in tNDPPP.in_MS.tgz.
// The MS contains 4 corr, 16 freq, 6 baselines, 18 time slots of 30 sec.
// Two time slots are missing between time slot 2 and 3.

namespace {
const std::string kInputMs = "../tNDPPP_tmp.MS";
const std::string kCopyMs = "tNDPPP_tmp.copy.MS";
const std::string kFlaggedMs = "tNDPPP_tmp.flag.MS";
const std::string kParsetFile = "tDP3.parset";
const std::size_t kNBaselines = 6;
const double kInterval = 30;

// The copy has six extra time slots:
// - Two slots are added since the input has two missing time slots.
// - One time slot is prepended to the ms. See the test_copy test.
// - Three time slots are appended to the ms. See the test_copy test.
const size_t kInputMsTimeSlots = 18;
const size_t kNCopyMsTimeSlots = kInputMsTimeSlots + 6;

void CheckCopy(const std::string& out_ms, std::vector<bool> ms_flagged) {
  const std::size_t n_ms = ms_flagged.size();
  Table table_in{kInputMs};
  Table table_out{out_ms};
  BOOST_CHECK_EQUAL(table_out.nrow(), kNCopyMsTimeSlots * kNBaselines);
  for (std::size_t j = 0; j < n_ms; ++j) {
    // A few dummy time slots were inserted, so ignore those.
    Table t1 = tableCommand(
        "using style python "
        "select from $1 where rownumber() not in [0:6, 4*6:6*6, 21*6:24*6]",
        table_out);
    ArrayColumn<Complex> data(t1, "DATA");
    ArrayColumn<bool> flag(t1, "FLAG");
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    ArrayColumn<float> weight(t1, "WEIGHT_SPECTRUM");
    IPosition dshape(2, 4, 16);
    IPosition dst(2, 0, j * 16);
    Slicer dslicer(dst, dshape);
    IPosition oflag_shape(2, 16 / 8, 1);
    IPosition oflag_dst(2, j * 16 / 8, 0);
    Slicer oflag_slicer(oflag_dst, oflag_shape);
    BOOST_CHECK_EQUAL(data(0).shape(), IPosition(2, 4, 16 * n_ms));
    BOOST_CHECK_EQUAL(flag(0).shape(), IPosition(2, 4, 16 * n_ms));
    BOOST_CHECK_EQUAL(weight(0).shape(), IPosition(2, 4, 16 * n_ms));
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 16 * n_ms / 8, 1));
    if (ms_flagged[j]) {
      BOOST_CHECK(allEQ(data.getColumn(dslicer), Complex()));
      BOOST_CHECK(allEQ(flag.getColumn(dslicer), true));
      BOOST_CHECK(allEQ(oflag.getColumn(oflag_slicer), (unsigned char)(255)));
      BOOST_CHECK(allEQ(weight.getColumn(dslicer), 0.0f));
    } else {
      BOOST_CHECK(allEQ(data.getColumn(dslicer),
                        ArrayColumn<Complex>(table_in, "DATA").getColumn()));
      BOOST_CHECK(allEQ(flag.getColumn(dslicer), false));
      BOOST_CHECK(allEQ(oflag.getColumn(oflag_slicer), (unsigned char)(0)));
      BOOST_CHECK(allEQ(weight.getColumn(dslicer), 1.0f));
    }

    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "TIME").getColumn(),
                      ScalarColumn<double>(table_in, "TIME").getColumn()));
    BOOST_CHECK(
        allEQ(ScalarColumn<double>(t1, "TIME_CENTROID").getColumn(),
              ScalarColumn<double>(table_in, "TIME_CENTROID").getColumn()));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "INTERVAL").getColumn(),
                      ScalarColumn<double>(table_in, "INTERVAL").getColumn()));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "EXPOSURE").getColumn(),
                      ScalarColumn<double>(table_in, "EXPOSURE").getColumn()));
    BOOST_CHECK(allEQ(ArrayColumn<double>(t1, "UVW").getColumn(),
                      ArrayColumn<double>(table_in, "UVW").getColumn()));
  }
  {
    // Check the inserted time slots.
    // The MS misses a few time slots (3 and 4).
    Table t1 = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [0:6, 4*6:6*6, 21*6:24*6]",
        table_out);
    ArrayColumn<Complex> data(t1, "DATA");
    ArrayColumn<bool> flag(t1, "FLAG");
    ArrayColumn<unsigned char> oflag(t1, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL(data(0).shape(), IPosition(2, 4, 16 * n_ms));
    BOOST_CHECK_EQUAL(flag(0).shape(), IPosition(2, 4, 16 * n_ms));
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 16 * n_ms / 8, 1));
    BOOST_CHECK(allEQ(data.getColumn(), Complex()));
    BOOST_CHECK(allEQ(flag.getColumn(), true));
    BOOST_CHECK(allEQ(oflag.getColumn(), (unsigned char)(0xff)));
    BOOST_CHECK(
        allEQ(ArrayColumn<float>(t1, "WEIGHT_SPECTRUM").getColumn(), float(0)));
    BOOST_CHECK(
        allEQ(ScalarColumn<double>(t1, "INTERVAL").getColumn(), kInterval));
    BOOST_CHECK(
        allEQ(ScalarColumn<double>(t1, "EXPOSURE").getColumn(), kInterval));
    double time = ScalarColumn<double>(table_in, "TIME")(0);
    for (unsigned int i = 0; i < kNBaselines; ++i) {
      double timec = time - kInterval;
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME")(i), timec));
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME_CENTROID")(i), timec));
    }
    time = ScalarColumn<double>(table_in, "TIME")(2 * 6);
    for (unsigned int i = kNBaselines; i < kNBaselines * 3; ++i) {
      double timec = time + (i / kNBaselines) * kInterval;
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME")(i), timec));
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME_CENTROID")(i), timec));
    }
    time = ScalarColumn<double>(table_in, "TIME")(17 * 6);
    for (unsigned int i = kNBaselines * 3; i < kNBaselines * 6; ++i) {
      double timec = time + (i / kNBaselines - 2) * kInterval;
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME")(i), timec));
      BOOST_CHECK(near(ScalarColumn<double>(t1, "TIME_CENTROID")(i), timec));
    }
  }
  // Now check if the SPECTRAL_WINDOW table is fine.
  Table spwin(table_in.keywordSet().asTable("SPECTRAL_WINDOW"));
  Table spwout(table_out.keywordSet().asTable("SPECTRAL_WINDOW"));
  for (std::size_t j = 0; j < n_ms; ++j) {
    IPosition dshape(1, 16);
    IPosition dst(1, j * 16);
    Slicer dslicer(dst, dshape);
    if (!ms_flagged[j]) {
      BOOST_CHECK(
          allEQ(ArrayColumn<double>(spwin, "CHAN_FREQ").getColumn(),
                ArrayColumn<double>(spwout, "CHAN_FREQ").getColumn(dslicer)));
    }
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
      double(n_ms) * ScalarColumn<double>(spwin, "TOTAL_BANDWIDTH").getColumn(),
      ScalarColumn<double>(spwout, "TOTAL_BANDWIDTH").getColumn()));
  if (n_ms == 1) {
    BOOST_CHECK(
        allEQ(ScalarColumn<double>(spwin, "REF_FREQUENCY").getColumn(),
              ScalarColumn<double>(spwout, "REF_FREQUENCY").getColumn()));
  }
  // Check the TIME_RANGE in the OBSERVATION table.
  Table obsout(table_out.keywordSet().asTable("OBSERVATION"));
  casacore::Vector<double> timeRange(
      ArrayColumn<double>(obsout, "TIME_RANGE").getColumn());
  BOOST_CHECK(near(timeRange(0), ScalarColumn<double>(table_out, "TIME")(0) -
                                     (kInterval / 2)));

  // In test_copy (when n_ms == 1), the end time is the start time plus a
  // multiple of the time interval. The start time is 13:22:20, the (rounded)
  // end time is 13:33:20 so the time range ends at 13:33:35.
  // The last time slot has time 13:33:15 and is thus 20 seconds before that.
  // In test_multi_in (when n_ms == 2), the parset does not specify an end time.
  // No rounding occurs then and the end time is the last time plus half an
  // interval, which is 30/2 = 15 seconds.
  const size_t end_time_increment = (n_ms == 1) ? 20 : (kInterval / 2);
  BOOST_CHECK(near(timeRange(1), ScalarColumn<double>(table_out, "TIME")(143) +
                                     end_time_increment));
}

}  // namespace

// Common part of test_copy and test_multi_in, which allows reusing CheckCopy().
void CreateCopyMs() {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    // Give starttime 35 sec before the first time, hence 1 missing timeslot.
    ostr << "msin.starttime=03-Aug-2000/13:21:45\n";
    // Give endtime 120 sec after the last time. MSReader rounds the end time
    // to the start time of the MS + a multiple of the time interval, downwards.
    // The start time is 13:22:20, so the rounded time becomes 13:33:20,
    // hence 3 missing timeslots, and not 4!
    ostr << "msin.endtime=03-Aug-2000/13:33:45\n";
    ostr << "msout=" << kCopyMs << '\n';
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);
}

BOOST_FIXTURE_TEST_CASE(test_copy, FixtureDirectory) {
  CreateCopyMs();
  CheckCopy(kCopyMs, {false});
}

BOOST_FIXTURE_TEST_CASE(test_copy_column, FixtureDirectory) {
  const std::string kCopyColumnMS = "tNDPPP_tmp.copy_column.MS";

  // Copying columns is only possible when updating an MS, not when creating
  // a new MS. -> First copy the MS, then copy columns.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "msout=" << kCopyColumnMS << '\n';
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);

  // First copy standard data and weight columns
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kCopyColumnMS << '\n';
    ostr << "msout=.\n";
    ostr << "msout.datacolumn=COPY_DATA\n";
    ostr << "msout.weightcolumn=COPY_WEIGHT\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);

  // Copy custom data and weight columns.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kCopyColumnMS << '\n';
    ostr << "msin.datacolumn=COPY_DATA\n";
    ostr << "msin.weightcolumn=COPY_WEIGHT\n";
    ostr << "msout=.\n";
    ostr << "msout.datacolumn=COPY2_DATA\n";
    ostr << "msout.weightcolumn=COPY2_WEIGHT\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);

  Table table_copy(kCopyColumnMS);
  BOOST_CHECK_GT(table_copy.nrow(), 0u);
  const ArrayColumn<std::complex<float>> data(table_copy, "DATA");
  const ArrayColumn<std::complex<float>> data_copy(table_copy, "COPY_DATA");
  const ArrayColumn<std::complex<float>> data_copy2(table_copy, "COPY2_DATA");
  const ArrayColumn<float> weight(table_copy, "WEIGHT_SPECTRUM");
  const ArrayColumn<float> weight_copy(table_copy, "COPY_WEIGHT");
  const ArrayColumn<float> weight_copy2(table_copy, "COPY2_WEIGHT");
  BOOST_CHECK(allEQ(data.getColumn(), data_copy.getColumn()));
  BOOST_CHECK(allEQ(data.getColumn(), data_copy2.getColumn()));
  BOOST_CHECK(allEQ(weight.getColumn(), weight_copy.getColumn()));
  BOOST_CHECK(allEQ(weight.getColumn(), weight_copy2.getColumn()));
}

BOOST_FIXTURE_TEST_CASE(test_multi_in_basic, FixtureDirectory) {
  const std::string kMultiMS = "tNDPPP_tmp.plaincopy.MS";

  // Creating the same MS as test_copy allows re-using CheckCopy here.
  // Copying is necessary anyway, since the MultiMSReader cannot handle
  // missing time slots.
  CreateCopyMs();

  // Test basic reading of two inputs MS's into one output MS.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=[" << kCopyMs << ", " << kCopyMs << "]\n";
    ostr << "msout=" << kMultiMS << '\n';
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);
  CheckCopy(kMultiMS, {false, false});
}

BOOST_FIXTURE_TEST_CASE(test_multi_in_missing_data, FixtureDirectory) {
  const std::string kMissingDataMS = "tNDPPP_tmp.missingdata.MS";
  const size_t kNBaselinesRemaining = 2;

  CreateCopyMs();

  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=[" << kCopyMs << ", " << kCopyMs << "]\n";
    ostr << "msin.datacolumn=MISSING_DATA\n";
    ostr << "msin.weightcolumn=MISSING_WEIGHT_SPECTRUM\n";
    ostr << "msin.missingdata=true\n";
    ostr << "msin.baseline=0,2&6\n";
    ostr << "msout=" << kMissingDataMS << '\n';
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);
  const Table tab(kMissingDataMS);
  BOOST_CHECK_EQUAL(tab.nrow(), kNCopyMsTimeSlots * kNBaselinesRemaining);
  BOOST_CHECK(allEQ(ArrayColumn<Complex>(tab, "DATA").getColumn(), Complex()));
  BOOST_CHECK(
      allEQ(ArrayColumn<float>(tab, "WEIGHT_SPECTRUM").getColumn(), 1.0f));
  BOOST_CHECK(allEQ(ArrayColumn<bool>(tab, "FLAG").getColumn(), true));
  BOOST_CHECK(allEQ(ScalarColumn<int>(tab, "ANTENNA2").getColumn(), 6));
}

BOOST_FIXTURE_TEST_CASE(test_multi_in_missing_ms, FixtureDirectory) {
  const std::string kMultiMS = "tNDPPP_tmp.multi.MS";

  CreateCopyMs();

  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=[notexist, " << kCopyMs << ", notexist, notexist]\n";
    ostr << "msin.orderms=false\n";
    ostr << "msout=" << kMultiMS << '\n';
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);
  CheckCopy(kMultiMS, {true, false, true, true});
}

void CheckAvg(const std::string& out_ms) {
  Table table_in{kInputMs};
  Table table_out{out_ms};
  ScalarColumn<double> timeCol(table_in, "TIME");
  double time = 0.5 * (timeCol(table_in.nrow() - 1) + timeCol(0));
  BOOST_CHECK_EQUAL(table_out.nrow(), kNBaselines);
  casacore::Block<casacore::String> colNames(2);
  colNames[0] = "ANTENNA1";
  colNames[1] = "ANTENNA2";
  TableIterator iterin(table_in, colNames);
  TableIterator iterout(table_out, colNames);
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
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "INTERVAL").getColumn(),
                      20 * kInterval));
    BOOST_CHECK(allEQ(ScalarColumn<double>(t1, "EXPOSURE").getColumn(),
                      20 * kInterval));
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
  Table spwin(table_in.keywordSet().asTable("SPECTRAL_WINDOW"));
  Table spwout(table_out.keywordSet().asTable("SPECTRAL_WINDOW"));
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
  Table obsout(table_out.keywordSet().asTable("OBSERVATION"));
  Vector<double> timeRange(
      ArrayColumn<double>(obsout, "TIME_RANGE").getColumn());
  BOOST_CHECK(
      near(timeRange(0), ScalarColumn<double>(table_out, "TIME")(0) - 300));
  BOOST_CHECK(
      near(timeRange(1), ScalarColumn<double>(table_out, "TIME")(0) + 295));
}

// Test averaging in a single step.
BOOST_FIXTURE_TEST_CASE(test_avg_single, FixtureDirectory) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin.name=" << kInputMs << '\n';
    ostr << "msout.name=tNDPPP_tmp.avg.MS\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[avg,count]\n";
    ostr << "avg.type=average\n";
    ostr << "avg.timestep=20\n";
    ostr << "avg.freqstep=100\n";
  }
  dp3::base::DP3::execute(kParsetFile);
  CheckAvg("tNDPPP_tmp.avg.MS");
}

// Averaging in multiple steps should be the same as above.
BOOST_FIXTURE_TEST_CASE(test_avg_multiple, FixtureDirectory) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "msout=tNDPPP_tmp.avg.MS\n";
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
  dp3::base::DP3::execute(kParsetFile);
  CheckAvg("tNDPPP_tmp.avg.MS");
}

// Averaging in multiple steps with multiple outputs should be the same, too.
BOOST_FIXTURE_TEST_CASE(test_avg_multiple_outputs, FixtureDirectory) {
  {
    std::ofstream ostr1("tNDPPP_tmp.parset1");
    std::ofstream ostr2("tNDPPP_tmp.parset2");
    std::ofstream ostr3("tNDPPP_tmp.parset3");
    std::ofstream ostr4("tNDPPP_tmp.parset4");
    ostr1 << "checkparset=1\n";
    ostr1 << "msin=" << kInputMs << '\n';
    ostr1 << "msout=tNDPPP_tmp.avg1.MS\n";
    ostr1 << "msout.overwrite=true\n";
    ostr1 << "steps=[avg1]\n";
    ostr1 << "avg1.type=average\n";
    ostr1 << "avg1.timestep=5\n";
    ostr1 << "avg1.freqstep=2\n";

    ostr2 << "checkparset=1\n";
    ostr2 << "msin=tNDPPP_tmp.avg1.MS\n";
    ostr2 << "msout=tNDPPP_tmp.avg2.MS\n";
    ostr2 << "msout.overwrite=true\n";
    ostr2 << "steps=[avg2]\n";
    ostr2 << "avg2.type=average\n";
    ostr2 << "avg2.timestep=1\n";
    ostr2 << "avg2.freqstep=2\n";

    ostr3 << "checkparset=1\n";
    ostr3 << "msin=tNDPPP_tmp.avg2.MS\n";
    ostr3 << "msout=tNDPPP_tmp.avg3.MS\n";
    ostr3 << "msout.overwrite=true\n";
    ostr3 << "steps=[avg3]\n";
    ostr3 << "avg3.type=average\n";
    ostr3 << "avg3.timestep=2\n";
    ostr3 << "avg3.freqstep=1\n";

    ostr4 << "checkparset=1\n";
    ostr4 << "msin=tNDPPP_tmp.avg3.MS\n";
    ostr4 << "msout=tNDPPP_tmp.avg4.MS\n";
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
  CheckAvg("tNDPPP_tmp.avg4.MS");
}

// This function tests if the correct start time is used when selecting times.
BOOST_FIXTURE_TEST_CASE(test_avg_start_time, FixtureDirectory) {
  {
    // Average in a single step.
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    // Give start a few seconds after first one, hence skip first time slot.
    ostr << "msin.starttime=03-Aug-2000/13:22:25\n";
    ostr << "msout=tNDPPP_tmp.avg.MS\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[avg]\n";
    ostr << "avg.type=average\n";
    ostr << "avg.timestep=2\n";
    ostr << "avg.freqstep=100\n";
  }
  dp3::base::DP3::execute(kParsetFile);

  // Only check the times in the averaged MS.
  Table table_input(kInputMs);
  // Take the first time to be used from the second time slot.
  double expected_time =
      ScalarColumn<double>(table_input, "TIME")(kNBaselines) + kInterval / 2;
  Table table_avg("tNDPPP_tmp.avg.MS");

  // -1: Because this test skips the first time slot.
  // +1: For rounding (kInputMsTimeSlots - 1) / 2 up.
  const std::size_t kExpectedTimeSlots = (kInputMsTimeSlots - 1 + 1) / 2;
  BOOST_CHECK_EQUAL(table_avg.nrow(), kNBaselines * kExpectedTimeSlots);

  ScalarColumn<double> time_column(table_avg, "TIME");
  for (unsigned int i = 0; i < table_avg.nrow(); ++i) {
    BOOST_CHECK_CLOSE(time_column(i), expected_time, 1e-6);
    if ((i % kNBaselines) == (kNBaselines - 1)) expected_time += kInterval * 2;
  }
}

namespace {

/**
 * Fixture for update tests.
 *
 * Besides creating a temporary directory, it also copies kInputMs into the
 * temporary directory as kCopyMs
 */
class FixtureCopyInput : public FixtureDirectory {
 public:
  FixtureCopyInput() : FixtureDirectory() {
    Table(kInputMs).deepCopy(kCopyMs, Table::New);
  }
};

}  // namespace

// Test if updating works fine.
// In fact, it does not do anything apart from rewriting the current flags.
// However, it should ignore the inserted time slots.
BOOST_FIXTURE_TEST_CASE(test_update_basic, FixtureCopyInput) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kCopyMs << '\n';
    ostr << "msout=''\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);

  // Check that the flags did not change.
  const Table table_input(kInputMs);
  const Table table_copy(kCopyMs);
  const ArrayColumn<bool> expected_flags(table_input, "FLAG");
  const ArrayColumn<bool> flags(table_copy, "FLAG");
  BOOST_CHECK(allEQ(flags.getColumn(), expected_flags.getColumn()));
}

// Test if updating all flags works fine.
BOOST_FIXTURE_TEST_CASE(test_update_flags, FixtureCopyInput) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kCopyMs << '\n';
    ostr << "msout=.\n";
    ostr << "steps=[preflag]\n";
    ostr << "preflag.blmin=1e6\n";  // should flag all data
  }
  dp3::base::DP3::execute(kParsetFile);

  // Check that all flags are true.
  const Table table_copy(kCopyMs);
  const ArrayColumn<bool> flags(table_copy, "FLAG");
  BOOST_CHECK(allEQ(flags.getColumn(), true));
}

// Test if updating data works fine.
BOOST_FIXTURE_TEST_CASE(test_update_scale, FixtureCopyInput) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kCopyMs << '\n';
    ostr << "msout=" << kCopyMs << '\n';  // same name means update
    ostr << "steps=[scaledata]\n";
    ostr << "scaledata.coeffs=2\n";
    ostr << "scaledata.stations=*\n";
    ostr << "scaledata.scalesize=false\n";
  }
  dp3::base::DP3::execute(kParsetFile);

  // Check that all data is doubled and that the flags remain equal.
  const Table table_input(kInputMs);
  const ArrayColumn<Complex> data_input(table_input, "DATA");
  const ArrayColumn<bool> flags_input(table_input, "FLAG");

  const Table table_output(kCopyMs);
  const ArrayColumn<Complex> data_output(table_output, "DATA");
  const ArrayColumn<bool> flags_output(table_output, "FLAG");

  BOOST_CHECK(allNear(data_input.getColumn() * Complex(2, 0),
                      data_output.getColumn(), 1e-6));
  BOOST_CHECK(allEQ(flags_input.getColumn(), flags_output.getColumn()));
}

namespace {

void CheckFullResFlags(const std::string& out_ms) {
  // Only check the FULL_RES_FLAGS and table size.
  // The flags are created in various ways, but should be the same in all cases.
  // 3 time slots are averaged to 1.
  const std::size_t kNOutputTimeSlots = 4;
  Table table_out(out_ms);
  BOOST_CHECK_EQUAL(table_out.nrow(), kNBaselines * kNOutputTimeSlots);
  // Check the full-res-flags.
  // Channels 0,2,6,7,8 are flagged everywhere.
  {
    // Input time slots 3,4 are inserted, thus flagged.
    Table table = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [0:6]",
        table_out);
    ArrayColumn<unsigned char> oflag(table, "LOFAR_FULL_RES_FLAG");
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
    Table table = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [6:12]",
        table_out);
    ArrayColumn<unsigned char> oflag(table, "LOFAR_FULL_RES_FLAG");
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
    Table table = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [12:18]",
        table_out);
    ArrayColumn<unsigned char> oflag(table, "LOFAR_FULL_RES_FLAG");
    BOOST_CHECK_EQUAL(oflag(0).shape(), IPosition(2, 2, 6));
    casacore::Array<unsigned char> flags = oflag.getColumn();
    BOOST_CHECK(allEQ(flags(IPosition(3, 0, 0, 0), IPosition(3, 0, 5, 5)),
                      (unsigned char)(0xc5)));
    BOOST_CHECK(allEQ(flags(IPosition(3, 1, 0, 0), IPosition(3, 1, 5, 5)),
                      (unsigned char)(0x01)));
  }
  {
    // Input time slot 20-23 did not exist, thus flagged in average.
    Table table = tableCommand(
        "using style python "
        "select from $1 where rownumber() in [18:24]",
        table_out);
    ArrayColumn<unsigned char> oflag(table, "LOFAR_FULL_RES_FLAG");
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

// Common part of test_flags_basic() and test_clear().
void CreateFlaggedMs() {
  // Just flag some channels and time stamps.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "msin.startchan=1\n";
    ostr << "msin.nchan=12\n";
    ostr << "msout=" << kFlaggedMs << '\n';
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[preflag,average]\n";
    ostr << "preflag.expr='flag1 or flag2'\n";
    ostr << "preflag.flag1.timeslot=[10,11]\n";
    ostr << "preflag.flag2.chan=[0,2,6..8]\n";
    ostr << "average.timestep=6\n";
  }
  dp3::base::DP3::execute(kParsetFile);
}

}  // namespace

BOOST_FIXTURE_TEST_CASE(test_flags_basic, FixtureDirectory) {
  CreateFlaggedMs();
  CheckFullResFlags(kFlaggedMs);
}

BOOST_FIXTURE_TEST_CASE(test_flags_shifted, FixtureDirectory) {
  const std::string kIntermediateMs = "tNDPPP_tmp.intermediate.MS";
  {
    // A more advanced case in two DP3 executions.
    // Flag the channels, but shifted 2. In the next step the first 2 channels
    // are skipped, thus the channel numbers shift back 2.
    // An averaged time slot is flagged, so the original time slots should
    // come out flagged.
    std::ofstream ostr1("tNDPPP_tmp.parset1");
    ostr1 << "checkparset=1\n";
    ostr1 << "msin=" << kInputMs << '\n';
    ostr1 << "msin.nchan=15\n";
    ostr1 << "msout=" << kIntermediateMs << '\n';
    ostr1 << "msout.overwrite=true\n";
    ostr1 << "steps=[preflag,average]\n";
    ostr1 << "preflag.chan=[2,4,8..10]\n";
    ostr1 << "average.timestep=2\n";
    std::ofstream ostr2("tNDPPP_tmp.parset2");
    ostr2 << "checkparset=1\n";
    ostr2 << "msin=" << kIntermediateMs << '\n';
    ostr2 << "msin.startchan=2*1\n";  // output chan 0,2 are now flagged
    ostr2 << "msin.nchan=nchan-3\n";
    ostr2 << "msout=" << kFlaggedMs << '\n';
    ostr2 << "msout.overwrite=true\n";
    ostr2 << "steps=[preflag,average]\n";
    ostr2 << "preflag.timeslot=5\n";  // is 10,11 in input
    ostr2 << "average.timestep=3\n";
    ostr2 << "average.freqstep=2\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset1");
  dp3::base::DP3::execute("tNDPPP_tmp.parset2");
  CheckFullResFlags(kFlaggedMs);
}

BOOST_FIXTURE_TEST_CASE(test_flags_averaged_channel, FixtureDirectory) {
  const std::string kIntermediateMs = "tNDPPP_tmp.intermediate.MS";
  {
    // Even a bit more advanced, also in two DP3 executions.
    // Input channels 6,7,8 are flagged by flagging their averaged channel.
    // This is done at the end of run 1, so the averager of run 2 should pick
    // up those flags.
    std::ofstream ostr1("tNDPPP_tmp.parset1");
    ostr1 << "checkparset=1\n";
    ostr1 << "msin=" << kInputMs << '\n';
    ostr1 << "msin.nchan=15\n";
    ostr1 << "msout=" << kIntermediateMs << '\n';
    ostr1 << "msout.overwrite=true\n";
    ostr1 << "steps=[preflag,average,pre2]\n";
    ostr1 << "preflag.chan=[0,2]\n";
    ostr1 << "average.timestep=2\n";
    ostr1 << "average.freqstep=3\n";
    ostr1 << "pre2.type=preflag\n";
    ostr1 << "pre2.chan=2\n";  // is input channel 6,7,8
    std::ofstream ostr2("tNDPPP_tmp.parset2");
    ostr2 << "checkparset=1\n";
    ostr2 << "msin=" << kIntermediateMs << '\n';
    ostr2 << "msin.nchan=4\n";
    ostr2 << "msout=" << kFlaggedMs << '\n';
    ostr2 << "msout.overwrite=true\n";
    ostr2 << "steps=[preflag,average]\n";
    ostr2 << "preflag.timeslot=5\n";  // is 10,11 in input
    ostr2 << "average.timestep=3\n";
    ostr2 << "average.freqstep=2\n";
  }
  dp3::base::DP3::execute("tNDPPP_tmp.parset1");
  dp3::base::DP3::execute("tNDPPP_tmp.parset2");
  CheckFullResFlags(kFlaggedMs);
}

BOOST_FIXTURE_TEST_CASE(test_flags_set_clear, FixtureDirectory) {
  const std::string kAllFlaggedMs = "tNDPPP_tmp.allflagged.MS";
  const std::string kAllClearMs = "tNDPPP_tmp.allclear.MS";

  // First flag in the same way as test_flags_basic.
  CreateFlaggedMs();

  // Check that the flags are not already set.
  BOOST_CHECK(
      !allEQ(ArrayColumn<bool>(Table(kFlaggedMs), "FLAG").getColumn(), true));

  // Flag all data.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kFlaggedMs << '\n';
    ostr << "msout=" << kAllFlaggedMs << '\n';
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[preflag]\n";
    ostr << "preflag.baseline=[[*]]\n";
  }
  dp3::base::DP3::execute(kParsetFile);
  CheckFullResFlags(kAllFlaggedMs);  // Full res flags shouldn't change.
  BOOST_CHECK(
      allEQ(ArrayColumn<bool>(Table(kAllFlaggedMs), "FLAG").getColumn(), true));

  // Clear the flags.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kAllFlaggedMs << '\n';
    ostr << "msout=" << kAllClearMs << '\n';
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[preflag]\n";
    ostr << "preflag.mode=clear\n";
    ostr << "preflag.baseline=[[*]]\n";
  }
  dp3::base::DP3::execute(kParsetFile);
  CheckFullResFlags(kAllClearMs);  // Full res flags shouldn't change.
  BOOST_CHECK(
      allEQ(ArrayColumn<bool>(Table(kAllClearMs), "FLAG").getColumn(), false));
}

BOOST_AUTO_TEST_CASE(test_station_add, *utf::disabled()) {
  // Add station RT0, 1 and 2.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=tNDPPP_tmp.MSa\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[stationadd]\n";
    ostr << "stationadd.stations={RTnew:[RT0..2]}\n";
  }
  dp3::base::DP3::execute(kParsetFile);
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

BOOST_AUTO_TEST_CASE(test_filter1, *utf::disabled()) {
  // Remove all baselines containing station RT1 or 6.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=tNDPPP_tmp.MSa\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[filter]\n";
    ostr << "filter.baseline=!RT[16]&&*\n";
    ostr << "filter.remove=true\n";
  }
  dp3::base::DP3::execute(kParsetFile);
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

BOOST_AUTO_TEST_CASE(test_filter2, *utf::disabled()) {
  // Keep all baselines.
  // First by not specifying baseline selection, second by all baselines.
  // Also alter between remove and !remove.
  for (int iter = 0; iter < 4; ++iter) {
    {
      std::ofstream ostr(kParsetFile);
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
    dp3::base::DP3::execute(kParsetFile);
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

BOOST_AUTO_TEST_CASE(test_filter3, *utf::disabled()) {
  // Remove some baselines, update original file with different data column
  // This test justs tests if it runs without throwing exceptions
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=tNDPPP_tmp.MS\n";
    ostr << "msout=.\n";
    ostr << "msout.datacolumn=DATA_FILTER\n";
    ostr << "steps=[filter]\n";
    ostr << "filter.baseline=!RT[16]&&*\n";
    ostr << "filter.remove=False\n";
  }
  dp3::base::DP3::execute(kParsetFile);
}

BOOST_FIXTURE_TEST_CASE(test_multi_out, FixtureDirectory) {
  {
    // First create an MS with the reference output.
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "msout=tNDPPP_tmp.MS_ref\n";
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[]\n";
  }
  dp3::base::DP3::execute(kParsetFile);

  // Test if update data works fine with multiple outputs:
  // Read from tNDPPP_tmp.MS, write to MS_out, update to MS_out.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "steps=[scaledata,out1,scaledata,out2]\n";
    ostr << "scaledata.coeffs=2\n";
    ostr << "scaledata.stations=*\n";
    ostr << "scaledata.scalesize=false\n";
    ostr << "out1.type=out\n";
    ostr << "out1.name=tNDPPP_tmp.MS_out\n";
    ostr << "out1.overwrite=true\n";
    ostr << "out2.type=out\n";
    ostr << "out2.name=.\n";  // Defaults to the previous out, so .MS_out.
    ostr << "out2.datacolumn=DATA_2\n";
  }
  dp3::base::DP3::execute(kParsetFile);

  // Check that the tables exist and that they contain the specified columns.
  Table table_ref("tNDPPP_tmp.MS_ref");
  Table table_out1("tNDPPP_tmp.MS_out");
  BOOST_CHECK(allNear(
      ArrayColumn<Complex>(table_out1, "DATA").getColumn(),
      Complex(2, 0) * ArrayColumn<Complex>(table_ref, "DATA").getColumn(),
      1e-5));
  BOOST_CHECK(allNear(
      ArrayColumn<Complex>(table_out1, "DATA_2").getColumn(),
      Complex(4, 0) * ArrayColumn<Complex>(table_ref, "DATA").getColumn(),
      1e-5));
  BOOST_CHECK(allNear(
      ArrayColumn<float>(table_out1, "WEIGHT_SPECTRUM").getColumn(),
      ArrayColumn<float>(table_ref, "WEIGHT_SPECTRUM").getColumn(), 1e-5));
  BOOST_CHECK(allEQ(ArrayColumn<bool>(table_out1, "FLAG").getColumn(),
                    ArrayColumn<bool>(table_ref, "FLAG").getColumn()));
}

BOOST_FIXTURE_TEST_CASE(test_error_out_avg, FixtureDirectory) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "steps=[filter,out1,average,out2]\n";
    ostr << "out1.type=out\n";
    ostr << "out1.name=''\n";
    ostr << "out2.type=out\n";
    ostr << "out2.name=.\n";  // update not possible when avg
    ostr << "msout=''\n";
  }
  BOOST_CHECK_THROW(dp3::base::DP3::execute(kParsetFile), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(test_error_out_filter, FixtureDirectory) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "steps=[average,out1,filter,out2]\n";
    ostr << "out1.type=out\n";
    ostr << "out1.name=tNDPPP_tmp.MSx\n";
    ostr << "out1.overwrite=true\n";
    ostr << "filter.remove=true\n";
    ostr << "out2.type=out\n";
    ostr << "out2.name=./tNDPPP_tmp.MSx\n";  // update not possible (filter)
    ostr << "msout=''\n";
  }
  BOOST_CHECK_THROW(dp3::base::DP3::execute(kParsetFile), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
