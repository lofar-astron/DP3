// tDP3.cc: test program for DP3
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "base/DP3.h"

#include "common/test/unit/fixtures/fDirectory.h"

#include <casacore/tables/Tables.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayPartMath.h>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>

#include "../LoggerFixture.h"

#include "steps/NullStep.h"
#include "steps/test/unit/mock/MockInput.h"
#include "steps/test/unit/mock/ThrowStep.h"

using casacore::ArrayColumn;
using casacore::Complex;
using casacore::IPosition;
using casacore::Matrix;
using casacore::near;
using casacore::ScalarColumn;
using casacore::Slicer;
using casacore::Table;
using casacore::TableIterator;

using dp3::common::Fields;
using dp3::common::test::FixtureDirectory;

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_SUITE(dp3, *utf::fixture<dp3::base::test::LoggerFixture>())

// This test program uses the MS in tNDPPP.in_MS.tgz.
// The MS contains 4 corr, 16 freq, 6 baselines, and 18 time slots of 30 sec.
// Two time slots are missing between time slot 2 and 3.

namespace {

// MockInput class with the getRequiredFields() and getProvidedFields() function
// implementations.
class TestInput : public steps::MockInput {
 public:
  common::Fields getRequiredFields() const override { return {}; };
  common::Fields getProvidedFields() const override { return {}; };
};

// Dummy step class where the required/provided fields can be set via the
// constructor.
class TestStep : public steps::test::ThrowStep {
 public:
  TestStep(Fields required_fields, Fields provided_fields)
      : required_fields_(required_fields), provided_fields_(provided_fields){};
  Fields getRequiredFields() const override { return required_fields_; };
  Fields getProvidedFields() const override { return provided_fields_; };

 private:
  Fields required_fields_;
  Fields provided_fields_;
};

class TestOutput : public steps::OutputStep {
 public:
  common::Fields getRequiredFields() const override {
    throw std::runtime_error("Unexpected getRequiredFields call");
  }
  common::Fields getProvidedFields() const override {
    throw std::runtime_error("Unexpected getProvidedFields call");
  }
  void updateInfo(const base::DPInfo&) override {
    BOOST_ERROR("Unexpected updateInfo() call");
  }
  void finish() override { BOOST_ERROR("Unexpected finish() call"); }
  void show(std::ostream&) const override {
    BOOST_ERROR("Unexpected show() call");
  }
};

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
  // AST-1186: Check that there are no symbolic links in the copy.
  for (const std::filesystem::directory_entry& entry :
       std::filesystem::recursive_directory_iterator(out_ms)) {
    BOOST_CHECK(!entry.is_symlink());
  }

  const std::size_t n_ms = ms_flagged.size();
  Table table_in{kInputMs};
  Table table_out{out_ms};
  BOOST_CHECK_EQUAL(table_out.nrow(), kNCopyMsTimeSlots * kNBaselines);
  for (std::size_t j = 0; j < n_ms; ++j) {
    // A few dummy time slots were inserted, so ignore those.
    Table t1 =
        tableCommand(
            "using style python "
            "select from $1 where rownumber() not in [0:6, 4*6:6*6, 21*6:24*6]",
            table_out)
            .table();
    ArrayColumn<Complex> data(t1, "DATA");
    ArrayColumn<bool> flag(t1, "FLAG");
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
    if (ms_flagged[j]) {
      BOOST_CHECK(allEQ(data.getColumn(dslicer), Complex()));
      BOOST_CHECK(allEQ(flag.getColumn(dslicer), true));
      BOOST_CHECK(allEQ(weight.getColumn(dslicer), 0.0f));
    } else {
      BOOST_CHECK(allEQ(data.getColumn(dslicer),
                        ArrayColumn<Complex>(table_in, "DATA").getColumn()));
      BOOST_CHECK(allEQ(flag.getColumn(dslicer), false));
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
    Table t1 =
        tableCommand(
            "using style python "
            "select from $1 where rownumber() in [0:6, 4*6:6*6, 21*6:24*6]",
            table_out)
            .table();
    ArrayColumn<Complex> data(t1, "DATA");
    ArrayColumn<bool> flag(t1, "FLAG");
    BOOST_CHECK_EQUAL(data(0).shape(), IPosition(2, 4, 16 * n_ms));
    BOOST_CHECK_EQUAL(flag(0).shape(), IPosition(2, 4, 16 * n_ms));
    BOOST_CHECK(allEQ(data.getColumn(), Complex()));
    BOOST_CHECK(allEQ(flag.getColumn(), true));
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
void CreateCopyMs(const std::string& input_ms) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << input_ms << '\n';
    // Give starttime 35 sec before the first time, hence 1 missing timeslot.
    ostr << "msin.starttime=03-Aug-2000/13:21:45\n";
    // Give endtime 120 sec after the last time. MsReader rounds the end time
    // to the start time of the MS + a multiple of the time interval, downwards.
    // The start time is 13:22:20, so the rounded time becomes 13:33:20,
    // hence 3 missing timeslots, and not 4!
    ostr << "msin.endtime=03-Aug-2000/13:33:45\n";
    ostr << "msout=" << kCopyMs << '\n';
    ostr << "msout.overwrite=true\n";
    ostr << "verbosity=QUIET\n";
    ostr << "steps=[]\n";
  }
  dp3::base::Execute(kParsetFile);
}

BOOST_FIXTURE_TEST_CASE(test_copy, FixtureDirectory) {
  CreateCopyMs(kInputMs);
  CheckCopy(kCopyMs, {false});
}

BOOST_FIXTURE_TEST_CASE(test_copy_symlinks, FixtureDirectory) {
  // AST-1186: Create a copy of the MS with symbolic links to read-only data,
  // and verify that DP3 can copy that MS correctly.
  const std::string kSymbolicMs = "tNDPPP_tmp.symbolic.MS";
  const std::string kReadOnlyMs = "tNDPPP_tmp.read-only.MS";

  std::filesystem::create_directory(kSymbolicMs);
  std::filesystem::create_directory(kReadOnlyMs);

  std::filesystem::copy(kInputMs, kSymbolicMs,
                        std::filesystem::copy_options::recursive);
  std::filesystem::copy(kInputMs, kReadOnlyMs,
                        std::filesystem::copy_options::recursive);

  const std::filesystem::path kReadOnlyPath =
      std::filesystem::current_path() / kReadOnlyMs;

  const std::filesystem::perms kWritePerms =
      std::filesystem::perms::owner_write |
      std::filesystem::perms::group_write |
      std::filesystem::perms::others_write;

  // Make everything in the read only directory read-only.
  for (const std::filesystem::directory_entry& entry :
       std::filesystem::recursive_directory_iterator(kReadOnlyPath)) {
    std::filesystem::permissions(entry.path(), kWritePerms,
                                 std::filesystem::perm_options::remove);
  }

  // Replace files in kSymbolicMs by symbolic links to the read-only directory.
  for (const std::filesystem::directory_entry& entry :
       std::filesystem::directory_iterator(kSymbolicMs)) {
    if (entry.is_regular_file()) {
      std::filesystem::remove(entry.path());
      std::filesystem::create_symlink(kReadOnlyPath / entry.path().filename(),
                                      entry.path());
    } else if (entry.is_directory()) {
      for (const std::filesystem::directory_entry& sub_entry :
           std::filesystem::directory_iterator(entry.path())) {
        // kInputMs has no sub-sub-directories.
        assert(sub_entry.is_regular_file());

        std::filesystem::remove(sub_entry.path());
        std::filesystem::create_symlink(kReadOnlyPath /
                                            entry.path().filename() /
                                            sub_entry.path().filename(),
                                        sub_entry.path());
      }
    }
  }

  CreateCopyMs(kSymbolicMs);
  CheckCopy(kCopyMs, {false});

  // Make everything in the read only directory writable again, otherwise
  // the fixture cannot remove it.
  for (const std::filesystem::directory_entry& entry :
       std::filesystem::recursive_directory_iterator(kReadOnlyMs)) {
    std::filesystem::permissions(entry.path(), kWritePerms,
                                 std::filesystem::perm_options::add);
  }
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
  dp3::base::Execute(kParsetFile);

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
  dp3::base::Execute(kParsetFile);

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
  dp3::base::Execute(kParsetFile);

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
  // Copying is necessary anyway, since the MultiMsReader cannot handle
  // missing time slots.
  CreateCopyMs(kInputMs);

  // Test basic reading of two inputs MS's into one output MS.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=[" << kCopyMs << ", " << kCopyMs << "]\n";
    ostr << "msout=" << kMultiMS << '\n';
    ostr << "steps=[]\n";
  }
  dp3::base::Execute(kParsetFile);
  CheckCopy(kMultiMS, {false, false});
}

BOOST_FIXTURE_TEST_CASE(test_multi_in_missing_data, FixtureDirectory) {
  const std::string kMissingDataMS = "tNDPPP_tmp.missingdata.MS";
  const size_t kNBaselinesRemaining = 2;

  CreateCopyMs(kInputMs);

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
  dp3::base::Execute(kParsetFile);
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

  CreateCopyMs(kInputMs);

  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=[notexist, " << kCopyMs << ", notexist, notexist]\n";
    ostr << "msin.orderms=false\n";
    ostr << "msout=" << kMultiMS << '\n';
    ostr << "steps=[]\n";
  }
  dp3::base::Execute(kParsetFile);
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
    BOOST_CHECK_EQUAL(data(0).shape(), IPosition(2, 4, 1));
    BOOST_CHECK_EQUAL(flag(0).shape(), IPosition(2, 4, 1));
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
  casacore::Vector<double> timeRange(
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
  dp3::base::Execute(kParsetFile);
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
  dp3::base::Execute(kParsetFile);
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
  dp3::base::Execute("tNDPPP_tmp.parset1");
  dp3::base::Execute("tNDPPP_tmp.parset2");
  dp3::base::Execute("tNDPPP_tmp.parset3");
  dp3::base::Execute("tNDPPP_tmp.parset4");
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
  dp3::base::Execute(kParsetFile);

  // Only check the times in the averaged MS.
  Table table_input(kInputMs);
  // Take the first time to be used from the second time slot.
  double expected_time =
      ScalarColumn<double>(table_input, "TIME")(kNBaselines) + kInterval / 2;
  Table table_avg("tNDPPP_tmp.avg.MS");

  // -1: Because this test skips the first time slot.
  // +2: Because the input has two missing time slots.
  // +1: For rounding (kInputMsTimeSlots - 1 + 2) / 2 up.
  const std::size_t kExpectedTimeSlots = (kInputMsTimeSlots - 1 + 2 + 1) / 2;
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
  dp3::base::Execute(kParsetFile);

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
  dp3::base::Execute(kParsetFile);

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
  dp3::base::Execute(kParsetFile);

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

// Test updating a subset of the channels, using startchan and nchan.
BOOST_FIXTURE_TEST_CASE(test_update_channel_subset, FixtureCopyInput) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kCopyMs << '\n';
    ostr << "msin.startchan=1\n";
    ostr << "msin.nchan=2\n";
    ostr << "msout=" << kCopyMs << '\n';  // same name means update
    ostr << "steps=[scaledata]\n";
    ostr << "scaledata.coeffs=2\n";
    ostr << "scaledata.stations=*\n";
    ostr << "scaledata.scalesize=false\n";
  }
  dp3::base::Execute(kParsetFile);

  const casacore::Array<Complex> data_input =
      ArrayColumn<Complex>(Table(kInputMs), "DATA").getColumn();
  const casacore::Array<Complex> data_output =
      ArrayColumn<Complex>(Table(kCopyMs), "DATA").getColumn();

  IPosition start = IPosition(3, 0, 0, 0);
  IPosition end = data_input.shape();
  const int n_channels = end[1];
  // casacore::Array includes the last element.
  --end[0];
  --end[2];

  // First channel should be unchanged.
  end[1] = 0;
  BOOST_CHECK(allNear(data_input(start, end), data_output(start, end), 1.0e-6));

  // Channels 1 and 2 should be scaled by a factor of 2.
  start[1] = 1;
  end[1] = 2;
  BOOST_CHECK(allNear(data_input(start, end) * Complex(2, 0),
                      data_output(start, end), 1.0e-6));

  // All remaining channels should be unchanged.
  start[1] = 3;
  end[1] = n_channels - 1;
  BOOST_CHECK(allNear(data_input(start, end), data_output(start, end), 1.0e-6));
}

// Test using msin.startchan/nchan together with filter.startchan/nchan, and
// then updating the input ms.
BOOST_FIXTURE_TEST_CASE(test_filter_update_channel_subset, FixtureCopyInput) {
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kCopyMs << '\n';
    ostr << "msin.startchan=4\n";
    ostr << "msin.nchan=4\n";
    ostr << "msout=" << kCopyMs << '\n';  // same name means update
    ostr << "steps=[filter,scaledata]\n";
    ostr << "filter.startchan=1\n";
    ostr << "filter.nchan=2\n";
    ostr << "scaledata.coeffs=2\n";
    ostr << "scaledata.stations=*\n";
    ostr << "scaledata.scalesize=false\n";
  }
  dp3::base::Execute(kParsetFile);

  const casacore::Array<Complex> data_input =
      ArrayColumn<Complex>(Table(kInputMs), "DATA").getColumn();
  const casacore::Array<Complex> data_output =
      ArrayColumn<Complex>(Table(kCopyMs), "DATA").getColumn();

  IPosition start = IPosition(3, 0, 0, 0);
  IPosition end = data_input.shape();
  const int n_channels = end[1];
  // casacore::Array includes the last element.
  --end[0];
  --end[2];

  // First 5 channels (indices 0-4) should be unchanged.
  end[1] = 4;
  BOOST_CHECK(allNear(data_input(start, end), data_output(start, end), 1.0e-6));

  // Channels 5 and 6 should be scaled by a factor of 2.
  start[1] = 5;
  end[1] = 6;
  BOOST_CHECK(allNear(data_input(start, end) * Complex(2, 0),
                      data_output(start, end), 1.0e-6));

  // All remaining channels should be unchanged.
  start[1] = 7;
  end[1] = n_channels - 1;
  BOOST_CHECK(allNear(data_input(start, end), data_output(start, end), 1.0e-6));
}

namespace {

casacore::Array<bool> GetFlags(const std::string& out_ms) {
  const std::size_t kNOutputTimeSlots = 4;
  Table table_out(out_ms);
  BOOST_CHECK_EQUAL(table_out.nrow(), kNBaselines * kNOutputTimeSlots);
  return ArrayColumn<bool>(table_out, "FLAG").getColumn();
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
    ostr << "preflag.flag1.timeslot=[10..17]\n";
    ostr << "preflag.flag2.chan=[0,2,6..8]\n";
    ostr << "average.timestep=6\n";
  }
  dp3::base::Execute(kParsetFile);
}

}  // namespace

BOOST_FIXTURE_TEST_CASE(test_flags_basic, FixtureDirectory) {
  CreateFlaggedMs();

  const casacore::Array<bool> flags = GetFlags(kFlaggedMs);
  BOOST_REQUIRE_EQUAL(flags.shape(), IPosition(3, 4, 12, 24));

  // Channels 0,2,6,7,8 are flagged everywhere.
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 0, 0), IPosition(3, 3, 0, 23)), true));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 2, 0), IPosition(3, 3, 2, 23)), true));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 6, 0), IPosition(3, 3, 8, 23)), true));

  // Input time slots 10-17 are flagged. After averaging, row 12-17 are flagged.
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 0, 12), IPosition(3, 3, 11, 17)), true));

  // All other flags are false.
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 1, 0), IPosition(3, 3, 1, 11)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 3, 0), IPosition(3, 3, 5, 11)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 9, 0), IPosition(3, 3, 11, 11)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 1, 18), IPosition(3, 3, 1, 23)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 3, 18), IPosition(3, 3, 5, 23)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 9, 18), IPosition(3, 3, 11, 23)), false));
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
    ostr2 << "preflag.timeslot=[5..8]\n";  // is 10..17 in input
    ostr2 << "average.timestep=3\n";
    ostr2 << "average.freqstep=2\n";
  }
  dp3::base::Execute("tNDPPP_tmp.parset1");
  dp3::base::Execute("tNDPPP_tmp.parset2");

  const casacore::Array<bool> flags = GetFlags(kFlaggedMs);
  BOOST_REQUIRE_EQUAL(flags.shape(), IPosition(3, 4, 6, 24));

  // Channels 0,2,6,7,8 are flagged everywhere. After averaging by a factor
  // of 2, only channel 3 (which was 6,7) is flagged.
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 3, 0), IPosition(3, 3, 3, 23)), true));

  // Input time slots 10-17 are flagged. After averaging, row 12-17 are flagged.
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 0, 12), IPosition(3, 3, 5, 17)), true));

  // All other flags are false.
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 0, 0), IPosition(3, 3, 2, 11)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 4, 0), IPosition(3, 3, 5, 11)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 0, 18), IPosition(3, 3, 2, 23)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 4, 18), IPosition(3, 3, 5, 23)), false));
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
    ostr2 << "preflag.timeslot=[5..8]\n";  // is [10..17] in input
    ostr2 << "average.timestep=3\n";
    ostr2 << "average.freqstep=2\n";
  }
  dp3::base::Execute("tNDPPP_tmp.parset1");
  dp3::base::Execute("tNDPPP_tmp.parset2");

  const casacore::Array<bool> flags = GetFlags(kFlaggedMs);
  BOOST_REQUIRE_EQUAL(flags.shape(), IPosition(3, 4, 2, 24));

  // Channels 0,2,6,7,8 are flagged everywhere. After averaging by a factor
  // of 6 (3 * 2), no channels are flagged anymore.

  // Input time slots 10-17 are flagged. After averaging, row 12-17 are flagged.
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 0, 12), IPosition(3, 3, 1, 17)), true));

  // All other flags are false.
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 0, 0), IPosition(3, 3, 1, 11)), false));
  BOOST_CHECK(
      allEQ(flags(IPosition(3, 0, 0, 18), IPosition(3, 3, 1, 23)), false));
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
  dp3::base::Execute(kParsetFile);
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
  dp3::base::Execute(kParsetFile);
  BOOST_CHECK(
      allEQ(ArrayColumn<bool>(Table(kAllClearMs), "FLAG").getColumn(), false));
}

BOOST_FIXTURE_TEST_CASE(test_station_add, FixtureDirectory) {
  const std::string kOutputMs = "tNDPPP_tmp.stationadd.MS";
  // Add station RT0, 1 and 2.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "msout=" << kOutputMs << '\n';
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[stationadd]\n";
    ostr << "stationadd.stations={RTnew:[RT0..2]}\n";
  }
  dp3::base::Execute(kParsetFile);

  Table antenna_in(kInputMs + "/ANTENNA");
  Table antenna_out(kOutputMs + "/ANTENNA");
  // Check that DP3 added 1 antenna.
  BOOST_CHECK_EQUAL(antenna_out.nrow(), antenna_in.nrow() + 1);
  BOOST_CHECK_EQUAL(ScalarColumn<casacore::String>(
                        antenna_out, "NAME")(antenna_out.nrow() - 1),
                    "RTnew");
  const int old_n_antennas = antenna_in.nrow();

  Table feed_in(kInputMs + "/FEED");
  Table feed_out(kOutputMs + "/FEED");
  BOOST_CHECK_EQUAL(feed_out.nrow(), feed_in.nrow() + 1);  // 1 antenna added
  BOOST_CHECK_EQUAL(
      ScalarColumn<int>(feed_out, "ANTENNA_ID")(feed_out.nrow() - 1),
      old_n_antennas);

  // Check that DP3 added:
  // - For each existing time slot, 2 baselines (kInputMsTimeSlots * 2).
  // - For each of the the two missing time slots in the MS, (kNBaselines + 2)
  //   baselines.
  BOOST_CHECK_EQUAL(
      Table(kOutputMs).nrow(),
      Table(kInputMs).nrow() + (kInputMsTimeSlots * 2) + (kNBaselines + 2) * 2);
}

BOOST_FIXTURE_TEST_CASE(test_filter_baseline, FixtureDirectory) {
  const std::string kOutputMs = "tNDPPP_tmp.filtered.MS";

  // Remove all baselines containing station RT1 or 6.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kInputMs << '\n';
    ostr << "msout=" << kOutputMs << '\n';
    ostr << "msout.overwrite=true\n";
    ostr << "steps=[filter]\n";
    ostr << "filter.baseline=!RT[16]&&*\n";
    ostr << "filter.remove=true\n";
  }
  dp3::base::Execute(kParsetFile);

  // Note: the ANTENNA table also contained RT8, RT9, etc., but they do not
  // have baselines. So these were removed as well meaning only 0,2,7 are left.
  casacore::Vector<dp3::common::rownr_t> rownrs(3);
  rownrs[0] = 0;
  rownrs[1] = 2;
  rownrs[2] = 7;

  Table antenna_in(kInputMs + "/ANTENNA");
  Table antenna_out(kOutputMs + "/ANTENNA");
  antenna_in = antenna_in(rownrs);
  BOOST_CHECK_EQUAL(antenna_out.nrow(), antenna_in.nrow());
  BOOST_CHECK(
      allEQ(ScalarColumn<casacore::String>(antenna_out, "NAME").getColumn(),
            ScalarColumn<casacore::String>(antenna_in, "NAME").getColumn()));

  Table feed_in(kInputMs + "/FEED");
  Table feed_out(kOutputMs + "/FEED");
  feed_in = feed_in(rownrs);
  BOOST_CHECK_EQUAL(feed_out.nrow(), feed_in.nrow());
  // The ANTENNA_IDs in the FEED table must be 0,1,2.
  casacore::Vector<int> ids(feed_out.nrow());
  indgen(ids);
  BOOST_CHECK(
      allEQ(ScalarColumn<int>(feed_out, "ANTENNA_ID").getColumn(), ids));

  // Check the main table.
  Table table_in(kInputMs);
  Table table_out(kOutputMs);
  // This test removes 4 baselines, thus kInputMsTimeSlots * 4 rows.
  // It adds 2 missing timeslots to the ms. Since there are 2 remaining
  // baselines, it thus adds 2 * 2 rows.
  BOOST_CHECK_EQUAL(table_out.nrow(),
                    table_in.nrow() - kInputMsTimeSlots * 4 + 2 * 2);

  Table table_in_rows = table_in(
      (table_in.col("ANTENNA1") == 0 || table_in.col("ANTENNA1") == 2) &&
      table_in.col("ANTENNA2") == 7);
  // The test inserts a few dummy time slots, so ignore those.
  Table table_out_rows =
      table_out(table_out.nodeRownr() < 6 || table_out.nodeRownr() >= 10);
  BOOST_CHECK(allEQ(ArrayColumn<Complex>(table_out_rows, "DATA").getColumn(),
                    ArrayColumn<Complex>(table_in_rows, "DATA").getColumn()));

  Table table_out_even = table_out(table_out.nodeRownr() % 2 == 0);
  BOOST_CHECK(
      allEQ(ScalarColumn<int>(table_out_even, "ANTENNA1").getColumn(), 0));
  BOOST_CHECK(
      allEQ(ScalarColumn<int>(table_out_even, "ANTENNA2").getColumn(), 2));

  Table table_out_odd = table_out(table_out.nodeRownr() % 2 == 1);
  BOOST_CHECK(
      allEQ(ScalarColumn<int>(table_out_odd, "ANTENNA1").getColumn(), 1));
  BOOST_CHECK(
      allEQ(ScalarColumn<int>(table_out_odd, "ANTENNA2").getColumn(), 2));
}

BOOST_FIXTURE_TEST_CASE(test_filter_keep_baselines, FixtureDirectory) {
  // Keep all baselines.
  // First by not specifying baseline selection, second by all baselines.
  // Also alter between remove and !remove.
  for (int iter = 0; iter < 4; ++iter) {
    const std::string kOutputMs =
        "tNDPPP_tmp.filtered." + std::to_string(iter) + ".MS";
    const bool filter_baseline = (iter % 2) == 1;
    const bool filter_remove = (iter / 2) == 1;
    {
      std::ofstream ostr(kParsetFile);
      ostr << "msin=" << kInputMs << '\n';
      ostr << "msout=" << kOutputMs << '\n';
      ostr << "msout.overwrite=true\n";
      ostr << "steps=[filter]\n";
      if (filter_baseline) {
        ostr << "filter.baseline=*&&*\n";
      }
      if (filter_remove) {
        ostr << "filter.remove=true\n";
      }
    }
    dp3::base::Execute(kParsetFile);

    // Note: the ANTENNA table also contained RT8, RT9, etc., but they do not
    // have baselines. So these were removed meaning only 0,1,2,6,7 are left.
    casacore::Vector<dp3::common::rownr_t> rownrs(5);
    rownrs[0] = 0;
    rownrs[1] = 1;
    rownrs[2] = 2;
    rownrs[3] = 6;
    rownrs[4] = 7;

    Table antenna_in(kInputMs + "/ANTENNA");
    Table antenna_out(kOutputMs + "/ANTENNA");
    if (filter_remove) {
      antenna_in = antenna_in(rownrs);
    }
    BOOST_CHECK_EQUAL(antenna_out.nrow(), antenna_in.nrow());
    BOOST_CHECK(
        allEQ(ScalarColumn<casacore::String>(antenna_out, "NAME").getColumn(),
              ScalarColumn<casacore::String>(antenna_in, "NAME").getColumn()));

    Table feed_in(kInputMs + "/FEED");
    Table feed_out(kOutputMs + "/FEED");
    if (filter_remove) {
      feed_in = feed_in(rownrs);
    }
    BOOST_CHECK_EQUAL(feed_out.nrow(), feed_in.nrow());
    // The ANTENNA_IDs in the FEED table must be 0,1,2.
    casacore::Vector<int> ids(feed_out.nrow());
    indgen(ids);
    BOOST_CHECK(
        allEQ(ScalarColumn<int>(feed_out, "ANTENNA_ID").getColumn(), ids));

    // Check the main table.
    Table main_in(kInputMs);
    Table main_out(kOutputMs);
    // The input MS has two missing timeslots, which should be added.
    BOOST_CHECK_EQUAL(main_out.nrow(), main_in.nrow() + 2 * kNBaselines);
    // A few dummy time slots were inserted, so ignore those.
    Table main_rows =
        main_out(main_out.nodeRownr() < 18 || main_out.nodeRownr() >= 30);
    BOOST_CHECK(allEQ(ArrayColumn<Complex>(main_rows, "DATA").getColumn(),
                      ArrayColumn<Complex>(main_in, "DATA").getColumn()));
    int ant1[] = {0, 0, 1, 1, 2, 2};
    int ant2[] = {6, 7, 6, 7, 6, 7};
    int sub = (filter_remove ? 3 : 0);  // if remove, ant2 6->3 and 7->4
    for (int i = 0; i < 6; ++i) {
      main_rows = main_out(main_out.nodeRownr() % 6 == i);
      BOOST_CHECK(
          allEQ(ScalarColumn<int>(main_rows, "ANTENNA1").getColumn(), ant1[i]));
      BOOST_CHECK(allEQ(ScalarColumn<int>(main_rows, "ANTENNA2").getColumn(),
                        ant2[i] - sub));
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_filter_different_data_column, FixtureCopyInput) {
  // Remove some baselines, update original file with different data column.
  // This test justs tests if it runs without throwing exceptions.
  {
    std::ofstream ostr(kParsetFile);
    ostr << "checkparset=1\n";
    ostr << "msin=" << kCopyMs << '\n';
    ostr << "msout=.\n";
    ostr << "msout.datacolumn=DATA_FILTER\n";
    ostr << "steps=[filter]\n";
    ostr << "filter.baseline=!RT[16]&&*\n";
    ostr << "filter.remove=False\n";
  }
  BOOST_CHECK_NO_THROW(dp3::base::Execute(kParsetFile));
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
  dp3::base::Execute(kParsetFile);

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
  dp3::base::Execute(kParsetFile);

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
  BOOST_CHECK_THROW(dp3::base::Execute(kParsetFile), std::runtime_error);
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
  BOOST_CHECK_THROW(dp3::base::Execute(kParsetFile), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_required_fields) {
  std::shared_ptr<TestInput> mock_input = std::make_shared<TestInput>();

  // step which requires data, uvw, provides data, uvw,
  std::shared_ptr<TestStep> step_read_write = std::make_shared<TestStep>(
      Fields(Fields::Single::kData) | Fields(Fields::Single::kUvw),
      Fields(Fields::Single::kData) | Fields(Fields::Single::kUvw));

  // step which requires flags, provides none
  std::shared_ptr<TestStep> step_read =
      std::make_shared<TestStep>(Fields(Fields::Single::kFlags), Fields());

  // step which requires uvw, provides data
  std::shared_ptr<TestStep> step_write = std::make_shared<TestStep>(
      Fields(Fields::Single::kUvw), Fields(Fields::Single::kData));

  // [step_write - step_read - step_read_write]
  // needs flags, uvw. Data is provided by the step_write step
  mock_input->setNextStep(step_write);
  step_write->setNextStep(step_read);
  step_read->setNextStep(step_read_write);
  step_read_write->setNextStep(std::make_shared<dp3::steps::NullStep>());

  const Fields kFlagsUvwFields =
      (Fields(Fields::Single::kFlags) | Fields(Fields::Single::kUvw));
  Fields overall_fields = dp3::base::GetChainRequiredFields(mock_input);
  BOOST_TEST(overall_fields == kFlagsUvwFields);

  // [step_read_write - step_read - step_write]
  // needs data, flags, uvw
  mock_input->setNextStep(step_read_write);
  step_read_write->setNextStep(step_read);
  step_read->setNextStep(step_write);
  step_write->setNextStep(std::make_shared<dp3::steps::NullStep>());

  overall_fields = dp3::base::GetChainRequiredFields(mock_input);
  BOOST_TEST(overall_fields ==
             (kFlagsUvwFields | Fields(Fields::Single::kData)));
}

BOOST_AUTO_TEST_CASE(test_provided_fields_single_step) {
  const Fields kProvidedFields =
      Fields(Fields::Single::kData) | Fields(Fields::Single::kFlags);
  const Fields kExtraField(Fields::Single::kUvw);

  auto output = std::make_shared<TestOutput>();
  BOOST_TEST(dp3::base::SetChainProvidedFields(output) == Fields());
  BOOST_TEST(output->GetFieldsToWrite() == Fields());

  BOOST_TEST(dp3::base::SetChainProvidedFields(output, kProvidedFields) ==
             Fields());
  BOOST_TEST(output->GetFieldsToWrite() == kProvidedFields);

  auto provides = std::make_shared<TestStep>(Fields(), kProvidedFields);
  BOOST_TEST(dp3::base::SetChainProvidedFields(provides) == kProvidedFields);
  BOOST_TEST(dp3::base::SetChainProvidedFields(provides, kExtraField) ==
             (kProvidedFields | kExtraField));
}

BOOST_AUTO_TEST_CASE(test_provided_fields_two_provides) {
  const Fields kProvidedFields1(Fields::Single::kWeights);
  const Fields kProvidedFields2(Fields::Single::kFlags);
  const Fields kAllProvidedFields = kProvidedFields1 | kProvidedFields2;

  auto provides1 = std::make_shared<TestStep>(Fields(), kProvidedFields1);
  auto provides2 = std::make_shared<TestStep>(Fields(), kProvidedFields2);

  provides1->setNextStep(provides2);
  BOOST_TEST(dp3::base::SetChainProvidedFields(provides1) ==
             kAllProvidedFields);

  auto output = std::make_shared<TestOutput>();
  provides2->setNextStep(output);
  BOOST_TEST(dp3::base::SetChainProvidedFields(provides1) == Fields());
  BOOST_TEST(output->GetFieldsToWrite() == kAllProvidedFields);
}

BOOST_AUTO_TEST_CASE(test_provided_fields_intermediate_output) {
  const Fields kProvidedFields1(Fields::Single::kWeights);
  const Fields kProvidedFields2(Fields::Single::kFlags);

  auto provides1 = std::make_shared<TestStep>(Fields(), kProvidedFields1);
  auto output1 = std::make_shared<TestOutput>();
  auto provides2 = std::make_shared<TestStep>(Fields(), kProvidedFields2);
  auto output2 = std::make_shared<TestOutput>();

  provides1->setNextStep(output1);
  output1->setNextStep(provides2);
  provides2->setNextStep(output2);

  BOOST_TEST(dp3::base::SetChainProvidedFields(provides1) == Fields());
  BOOST_TEST(output1->GetFieldsToWrite() == kProvidedFields1);
  BOOST_TEST(output2->GetFieldsToWrite() == kProvidedFields2);
}

BOOST_AUTO_TEST_SUITE_END()
