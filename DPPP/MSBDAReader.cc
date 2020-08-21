// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
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

#include "MSBDAReader.h"
#include "BDABuffer.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Exceptions.h"

#include <EveryBeam/load.h>
#include <EveryBeam/lofarreadutils.h>

#include "../Common/ParameterSet.h"
#include "../Common/BaselineSelect.h"
#include "../AOFlaggerStep/AOFlaggerStep.h"

#include <casacore/casa/OS/HostInfo.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSAntennaParse.h>
#include <casacore/ms/MSSel/MSSelectionErrorHandler.h>
#include <casacore/casa/version.h>

#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/OS/Conversion.h>

#include <iostream>
#include <boost/make_unique.hpp>

using namespace casacore;

namespace DP3 {
namespace DPPP {

MSBDAReader::MSBDAReader()
    : readVisData_(False), lastMSTime_(0), nread_(0), maxChanWidth_(0) {}

MSBDAReader::MSBDAReader(const string& msName, const ParameterSet& parset,
                         const string& prefix)
    : msName_(msName),
      readVisData_(False),
      lastMSTime_(0),
      nread_(0),
      maxChanWidth_(0) {
  NSTimer::StartStop sstime(timer_);
  spw_ = parset.getInt(prefix + "band", -1);
  dataColName_ = parset.getString(prefix + "datacolumn", "DATA");
  weightColName_ = parset.getString(prefix + "weightcolumn", "WEIGHT_SPECTRUM");
  memory_ = parset.getUint(prefix + "memorymax", 0);
  memory_percentage_ = parset.getUint(prefix + "memoryperc", 0);
}

MSBDAReader::~MSBDAReader() {}

DPStep::MSType MSBDAReader::outputs() const { return BDA; };

void MSBDAReader::updateInfo(const DPInfo& dpInfo) {
  info().setNThreads(dpInfo.nThreads());
  DetermineAvailableMemory();

  if (!Table::isReadable(msName_)) {
    throw std::invalid_argument("No such MS: " + msName_);
  }

  ms_ = MeasurementSet(msName_, TableLock::AutoNoReadLocking);

  if (!ms_.keywordSet().isDefined("BDA_FACTORS") ||
      ms_.keywordSet().asTable("BDA_FACTORS").nrow() == 0) {
    throw std::invalid_argument(
        "Input MS does not contain BDA data. Table BDA_FACTORS is missing "
        "or not filled");
  }

  // Find the nr of corr, chan, and baseline.
  IPosition shp(ArrayColumn<Complex>(ms_, "DATA").shape(0));
  ncorr_ = shp[0];

  // Read the meta data tables and store required values.
  FillInfoMetaData();

  // Read the antenna set.
  Table obstab(ms_.keywordSet().asTable("OBSERVATION"));
  string antennaSet;
  if (obstab.nrow() > 0 && obstab.tableDesc().isColumn("LOFAR_ANTENNA_SET")) {
    antennaSet = ScalarColumn<String>(obstab, "LOFAR_ANTENNA_SET")(0);
  }

  // Determine the pool size for the bda buffers
  std::size_t expected_row_size =
      (sizeof(float) * 2 + sizeof(casacore::Complex) + sizeof(double) * 2) *
          ncorr_ * maxChanWidth_ +
      3 * sizeof(double) + 4 * sizeof(std::size_t) + 3 * sizeof(double);
  std::size_t rowsPerBuffer = memory_avail_ / expected_row_size;
  poolSize_ = std::max(std::size_t(1), rowsPerBuffer);

  double startTime = 0;
  unsigned int startChan = 0;
  unsigned int ntimeAprox = std::ceil(ms_.nrow() / (float)rowsPerBuffer);

  info().init(ncorr_, startChan, maxChanWidth_, ntimeAprox, startTime,
              interval_, msName(), antennaSet);
  info().setDataColName(dataColName_);
  info().setWeightColName(weightColName_);
}

std::string MSBDAReader::msName() const { return ms_.tableName(); }

void MSBDAReader::setReadVisData(bool readVisData) {
  readVisData_ = readVisData || readVisData_;
}

bool MSBDAReader::process(const DPBuffer&) {
  NSTimer::StartStop sstime(timer_);
  std::unique_ptr<BDABuffer> buffer = boost::make_unique<BDABuffer>(poolSize_);

  ScalarColumn<int> ant1Col(ms_, "ANTENNA1");
  ScalarColumn<int> ant2Col(ms_, "ANTENNA2");
  ArrayColumn<Complex> dataCol(ms_, "DATA");
  ArrayColumn<float> weightsCol(ms_, "WEIGHT_SPECTRUM");
  ArrayColumn<double> uvwCol(ms_, "UVW");

  while (nread_ < ms_.nrow() && buffer->GetRemainingCapacity() > 0) {
    // Take time from row 0 in subset.
    double msTime = ScalarColumn<double>(ms_, "TIME")(0);

    if (msTime < lastMSTime_) {
      DPLOG_WARN_STR("Time at rownr " + std::to_string(nread_) + " of MS " +
                     msName() + " is less than previous time slot");
      ++nread_;
      continue;
    }

    Cube<Complex> data = dataCol.get(nread_);
    Cube<float> weights = weightsCol.get(nread_);
    Cube<double> uvw = uvwCol.get(nread_);
    double interval = ScalarColumn<double>(ms_, "INTERVAL")(nread_);
    double exposure = ScalarColumn<double>(ms_, "EXPOSURE")(nread_);
    int dataDescId = ScalarColumn<int>(ms_, "DATA_DESC_ID")(nread_);

    size_t blNr = blToBLId_[std::make_pair(ant1Col(nread_), ant2Col(nread_))];

    buffer->AddRow(msTime, interval, exposure, blNr, descIdToNchan_[dataDescId],
                   info().ncorr(), data.tovector().data(), nullptr,
                   weights.tovector().data(), nullptr, uvw.tovector().data());

    lastMSTime_ = msTime;
    ++nread_;
  }

  getNextStep()->process(std::move(buffer));

  // Return true while there are still items remaining
  return nread_ < ms_.nrow();
}

void MSBDAReader::finish() { getNextStep()->finish(); }

void MSBDAReader::FillInfoMetaData() {
  Table factors = ms_.keywordSet().asTable("BDA_FACTORS");
  Table axis = ms_.keywordSet().asTable("BDA_TIME_AXIS");
  Table spw = ms_.keywordSet().asTable("SPECTRAL_WINDOW");

  interval_ = axis.col("UNIT_TIME_INTERVAL").getDouble(0);
  nbl_ = factors.nrow();

  // Required columns to read
  ScalarColumn<int> factorCol(factors, "FACTOR");
  ScalarColumn<int> ant1Col(factors, "ANTENNA1");
  ScalarColumn<int> ant2Col(factors, "ANTENNA2");
  ScalarColumn<int> idsCol(factors, "SPECTRAL_WINDOW_ID");
  ArrayColumn<double> freqsCol(spw, "CHAN_FREQ");
  ArrayColumn<double> widthsCol(spw, "CHAN_WIDTH");

  // Fill info with the data required to repopulate BDA_FACTORS
  std::vector<std::vector<double>> freqs(nbl_);
  std::vector<std::vector<double>> widths(nbl_);
  std::vector<unsigned int> baseline_factors(nbl_);
  for (unsigned int i = 0; i < nbl_; ++i) {
    unsigned int spwId = idsCol(i);
    baseline_factors[i] = factorCol(i);
    freqs[i] = freqsCol.get(spwId).tovector();
    widths[i] = widthsCol.get(spwId).tovector();

    maxChanWidth_ = std::max(maxChanWidth_, freqs[i].size());
    descIdToNchan_[idsCol(i)] = freqs[i].size();
    blToBLId_[std::make_pair(ant1Col(i), ant2Col(i))] = i;
  }

  Table anttab(ms_.keywordSet().asTable("ANTENNA"));
  ScalarColumn<String> nameCol(anttab, "NAME");
  ScalarColumn<Double> diamCol(anttab, "DISH_DIAMETER");
  ROScalarMeasColumn<MPosition> antCol(anttab, "POSITION");
  vector<MPosition> antPos;
  antPos.reserve(anttab.nrow());
  for (unsigned int i = 0; i < anttab.nrow(); ++i) {
    antPos.push_back(antCol(i));
  }

  // Set antenna/baseline info.
  info().set(nameCol.getColumn(), diamCol.getColumn(), antPos,
             ant1Col.getColumn(), ant2Col.getColumn());

  info().update(baseline_factors);
  info().set(std::move(freqs), std::move(widths));
}

void MSBDAReader::DetermineAvailableMemory() {
  // Determine available memory.
  double availMemory = casacore::HostInfo::memoryTotal() * 1024.;
  // Determine how much memory can be used.
  double memoryMax = memory_ * 1024 * 1024 * 1024;
  double memory = memoryMax;
  if (memory_percentage_ > 0) {
    memory = memory_percentage_ * availMemory / 100.;
    if (memoryMax > 0 && memory > memoryMax) {
      memory = memoryMax;
    }
  } else if (memory_ <= 0) {
    // Nothing given, so use available memory on this machine.
    // Set 50% (max 2 GB) aside for other purposes.
    memory = availMemory - std::min(0.5 * availMemory, 2. * 1024 * 1024 * 1024);
  }
  memory_avail_ = memory;
}

void MSBDAReader::show(std::ostream& os) const {
  os << "MSReader\n";
  os << "  input MS:       " << msName_ << '\n';
  if (ms_.isNull()) {
    os << "    *** MS does not exist ***\n";
  } else {
    os << "  band            " << spw_ << '\n';
    os << "  startchan:      " << 0 << '\n';
    os << "  nchan:          " << getInfo().nchan() << '\n';
    os << "  ncorrelations:  " << getInfo().ncorr() << '\n';
    os << "  nbaselines:     " << nbl_ << '\n';
    os << "  first time:     " << MVTime::Format(MVTime::YMD) << MVTime(0)
       << '\n';
    os << "  last time:      " << MVTime::Format(MVTime::YMD)
       << MVTime(lastMSTime_ / (24 * 3600.)) << '\n';
    os << "  ntimes:         " << getInfo().ntime() << '\n';
    os << "  time interval:  " << getInfo().timeInterval() << '\n';
    os << "  DATA column:    " << dataColName_;
    os << '\n';
    os << "  WEIGHT column:  " << weightColName_ << '\n';
    os << "  memory ";
    AOFlaggerStep::formatBytes(os, memory_avail_);
  }
}

void MSBDAReader::showCounts(std::ostream& os) const {
  os << endl << "NaN/infinite data flagged in reader";
  os << endl << "===================================" << endl;
  os << 0 << " missing time slots were inserted" << endl;
}

void MSBDAReader::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " MSBDAReader" << endl;
}

}  // namespace DPPP
}  // namespace DP3
