// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MSOper/MSSummary.h>
#include <casacore/tables/TaQL/TableParse.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Inputs/Input.h>

#include <iostream>
#include <sstream>

namespace {
void RunMSOverview(int argc, char* argv[]) {
  // enable input in no-prompt mode
  casacore::Input inputs(1);
  // define the input structure
  inputs.version("20110407GvD");
  inputs.create("in", "", "Name of input MeasurementSet", "string");
  inputs.create("verbose", "False", "Make a verbose listing?", "bool");
  // Fill the input structure from the command line.
  inputs.readArguments(argc, argv);

  // Get and check the input specification.
  const std::string msin(inputs.getString("in"));
  if (msin.empty()) {
    throw std::runtime_error("An input MeasurementSet must be given");
  }
  const bool verbose = inputs.getBool("verbose");

  // Show the MS summary.
  // Note that class MSSummary uses LogIO. Use that on a stringstream.
  const casacore::MeasurementSet ms(msin);
  casacore::MSSummary summ(ms);
  std::ostringstream ostr;
  casacore::LogSink logsink(casacore::LogMessage::NORMAL, &ostr, false);
  casacore::LogIO logio(logsink);
  summ.listTitle(logio);
  // Show if the MS is a reference or concatenation.
  casacore::Block<casacore::String> partNames = ms.getPartNames();
  if (partNames.size() == 1) {
    if (partNames[0] != ms.tableName()) {
      casacore::Table tab(partNames[0]);
      logio << casacore::LogIO::NORMAL << "           "
            << "The MS references " << ms.nrow() << " out of " << tab.nrow()
            << " rows in " << partNames[0] << casacore::LogIO::POST;
    } else {
      // Show if it is a raw LOFAR MS (stored with LofarStMan).
      casacore::Record dminfo = ms.dataManagerInfo();
      for (size_t i = 0; i < dminfo.nfields(); ++i) {
        casacore::Record subrec = dminfo.subRecord(i);
        if (subrec.asString("TYPE") == "LofarStMan") {
          logio << casacore::LogIO::NORMAL << "           "
                << "This is a raw LOFAR MS (stored with LofarStMan)"
                << casacore::LogIO::POST;
          break;
        }
      }
    }
  } else if (!partNames.empty()) {
    logio << casacore::LogIO::NORMAL << "           "
          << "The MS is a concatenation of: \n";
    for (const casacore::String& name : partNames) {
      casacore::Table tab(name);
      logio << "               " << name << "  (" << tab.nrow() << " rows) \n";
    }
    logio << casacore::LogIO::POST;
  }
  logio << casacore::LogIO::NORMAL << '\n' << casacore::LogIO::POST;
  summ.listWhere(logio, true);
  // If possible, show the AntennaSet.
  casacore::Table obsTab(ms.keywordSet().asTable("OBSERVATION"));
  if (obsTab.tableDesc().isColumn("LOFAR_ANTENNA_SET")) {
    logio << casacore::LogIO::NORMAL << "Antenna-set: "
          << casacore::ScalarColumn<casacore::String>(obsTab,
                                                      "LOFAR_ANTENNA_SET")(0)
          << casacore::LogIO::POST;
  }
  logio << casacore::LogIO::NORMAL << '\n' << casacore::LogIO::POST;
  summ.listMain(logio, false);
  logio << casacore::LogIO::NORMAL << '\n' << casacore::LogIO::POST;
  summ.listField(logio, false);
  logio << casacore::LogIO::NORMAL << '\n' << casacore::LogIO::POST;
  summ.listSpectralAndPolInfo(logio, verbose, false);  // new version
  if (verbose) {
    logio << casacore::LogIO::NORMAL << '\n' << casacore::LogIO::POST;
    summ.listAntenna(logio, true);
  }

  // Remove the extra fields (time, severity) from the output string.
  casacore::String str(ostr.str());
  str.gsub(casacore::Regex(".*\tINFO\t[+]?\t"), "");
  std::cout << str;

  // Test if the MS is regular.
  if (verbose) {
    size_t nrdd = ms.dataDescription().nrow();
    // An MS is regular if all times have same nr of baselines.
    const size_t nrtime =
        tableCommand("select from $1 orderby unique TIME", ms).table().nrow();
    casacore::Table blTab =
        tableCommand("select from $1 orderby unique ANTENNA1,ANTENNA2", ms)
            .table();
    const size_t nrbl = blTab.nrow();
    const size_t nrauto =
        tableCommand("select from $1 where ANTENNA1=ANTENNA2", blTab)
            .table()
            .nrow();
    tableCommand("select from $1 orderby unique ANTENNA1,ANTENNA2", ms)
        .table()
        .nrow();
    if (nrdd > 1) {
      // Get actual nr of bands.
      nrdd = tableCommand("select from $1 orderby unique DATA_DESC_ID", ms)
                 .table()
                 .nrow();
    }
    std::cout << '\n';
    if (ms.nrow() == nrtime * nrbl * nrdd) {
      std::cout << "The MS is fully regular.\n";
    } else {
      std::cout << "The MS is not regular.\n"
                   "  use msregularize in python casacore.tables to make it "
                   "regular\n";
    }
    std::cout << "   nrows=" << ms.nrow() << "   ntimes=" << nrtime
              << "   nbands=" << nrdd << "   nbaselines=" << nrbl << " ("
              << nrauto << " autocorr)\n";
  }
  std::cout << '\n';
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    RunMSOverview(argc, argv);
  } catch (std::exception& x) {
    std::cerr << "Error: " << x.what() << '\n';
    return 1;
  }
  return 0;
}
