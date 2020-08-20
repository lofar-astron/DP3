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
    : readVisData_(False), lastMSTime_(0), nread_(0), ninserted_(0) {}

MSBDAReader::MSBDAReader(const string& msName, const ParameterSet& parset,
                         const string& prefix)
    : msName_(msName),
      readVisData_(False),
      lastMSTime_(0),
      nread_(0),
      ninserted_(0) {
  NSTimer::StartStop sstime(timer_);
  spw_ = parset.getInt(prefix + "band", -1);
}

MSBDAReader::~MSBDAReader() {}

void MSBDAReader::updateInfo(const DPInfo& dpInfo) {
  info().setNThreads(dpInfo.nThreads());

  if (!Table::isReadable(msName_)) {
    throw std::invalid_argument("No such MS: " + msName_);
  }

  if (!ms_.keywordSet().isDefined("BDA_TIME_AXIS")) {
    throw std::invalid_argument(
        "Input MS does not contain BDA data. Table BDA_TIME_AXIS is missing");
  }

  ms_ = MeasurementSet(msName_, TableLock::AutoNoReadLocking);

  // Create iterator over time. Do not sort again.
  iter_ = TableIterator(ms_, Block<String>(1, "TIME"), TableIterator::Ascending,
                        TableIterator::NoSort);

  // Find the nr of corr, chan, and baseline.
  IPosition shp(ArrayColumn<Complex>(ms_, "DATA").shape(0));
  ncorr_ = shp[0];
  // TODO move to other part
  // itsNrChan = shp[1];
  // TODO if this is correct, then the iter is used wrong in the process
  nbl_ = iter_.table().nrow();

  // TODO initialize info
  // // Set antenna/baseline info.
  // info().set(nameCol.getColumn(), diamCol.getColumn(), antPos,
  //            ant1col.getColumn(), ant2col.getColumn());

  // if (itsAutoWeight) {
  //   info().setNeedVisData();
  //   info().setWriteWeights();
  // }
  // info().set(arrayPos, phaseCenter, delayCenter, tileBeamDir);
  // info().init(ncorr_, itsStartChan, itsNrChan, ntime, itsStartTime,
  //             interval_, msName(), antennaSet);
  // info().setDataColName(itsDataColName);
  // info().setWeightColName(itsWeightColName);
  // if (itsUseAllChan) {
  //   info().set(std::move(chanFreqs), std::move(chanWidths),
  //              std::move(resolutions), std::move(effectiveBW), refFreq);
  // } else {
  //   auto freqBegin = chanFreqs.begin() + itsStartChan;
  //   auto widthBegin = chanWidths.begin() + itsStartChan;
  //   auto resolBegin = resolutions.begin() + itsStartChan;
  //   auto effbwBegin = effectiveBW.begin() + itsStartChan;
  //   info().set(std::vector<double>(freqBegin, freqBegin + itsNrChan),
  //              std::vector<double>(widthBegin, widthBegin + itsNrChan),
  //              std::vector<double>(resolBegin, resolBegin + itsNrChan),
  //              std::vector<double>(effbwBegin, effbwBegin + itsNrChan),
  //              refFreq);
  // }
  FillMetaData();
}

std::string MSBDAReader::msName() const { return ms_.tableName(); }

void MSBDAReader::setReadVisData(bool readVisData) {
  readVisData_ = readVisData || readVisData_;
}

bool MSBDAReader::process(std::unique_ptr<BDABuffer>) {
  // TODO determine actual size
  // TODO cannot use info().nchan(), maybe use max nchan?
  size_t poolSize = ncorr_ * info().nchan();
  std::unique_ptr<BDABuffer> buffer = boost::make_unique<BDABuffer>(poolSize);

  while (!iter_.pastEnd() && buffer->GetRemainingCapacity() > 0) {
    // Take time from row 0 in subset.
    double msTime = ScalarColumn<double>(iter_.table(), "TIME")(0);

    if (msTime < lastMSTime_) {
      DPLOG_WARN_STR("Time at rownr " +
                     std::to_string(iter_.table().rowNumbers(ms_)[0]) +
                     " of MS " + msName() + " is less than previous time slot");
      continue;
    }

    casacore::Cube<casacore::Complex> data;
    casacore::Cube<float> weights;
    casacore::Cube<double> uvw;
    // TODO dont use iter
    double interval = ScalarColumn<double>(iter_.table(), "INTERVAL")(0);
    double exposure = ScalarColumn<double>(iter_.table(), "EXPOSURE")(0);
    int dataDescId = ScalarColumn<int>(iter_.table(), "DATA_DESC_ID")(0);
    ArrayColumn<casacore::Complex>(iter_.table(), "DATA").getColumn(data);
    ArrayColumn<float>(iter_.table(), "WEIGHT_SPECTRUM").getColumn(weights);
    ArrayColumn<double>(iter_.table(), "UVW").getColumn(uvw);

    // TODO for i in range(info().nbaselines) if ant1 == ANTENNA1 and ant2 ==
    // ANTENNA2
    size_t baselineNr = 0;
    // TODO make map from DATA_DESC_ID to nchan // from SPECTRAL_WINDOW and
    // DATA_DESCRIPTION
    buffer->AddRow(msTime, interval, exposure, baselineNr,
                   descIdToNchan_[dataDescId], info().ncorr(),
                   data.tovector().data(), nullptr, weights.tovector().data(),
                   nullptr, uvw.tovector().data());

    lastMSTime_ = msTime;
    iter_.next();
  }

  getNextStep()->process(std::move(buffer));

  // Return true while there are still items remaining
  return !iter_.pastEnd();
}

void MSBDAReader::finish() { getNextStep()->finish(); }

void MSBDAReader::show(std::ostream& os) const {}

void MSBDAReader::showCounts(std::ostream& os) const {}

void MSBDAReader::showTimings(std::ostream& os, double duration) const {}

void MSBDAReader::getUVW(const RefRows& rowNrs, double time, DPBuffer& buf) {}

void MSBDAReader::getWeights(const RefRows& rowNrs, DPBuffer& buf) {}

bool MSBDAReader::getFullResFlags(const RefRows& rowNrs, DPBuffer& buf) {
  return true;
}

void MSBDAReader::getModelData(const casacore::RefRows& rowNrs,
                               casacore::Cube<casacore::Complex>& arr) {}

void MSBDAReader::fillBeamInfo(
    vector<everybeam::Station::Ptr>& vec,
    const casacore::Vector<casacore::String>& antNames) {}

void MSBDAReader::FillMetaData() {
  // Fill info with the data required repopulate BDA_TIME_AXIS
  casacore::Table t = ms_.keywordSet().asTable("BDA_TIME_AXIS");
  if (t.nrow() > 1) {
    // TODO throw if the factors table is empty instead!
    throw std::runtime_error(
        "DP3 cannot handle multiple BDA_TIME_AXIS entries");
  }
  interval_ = t.col("UNIT_TIME_INTERVAL").getDouble(0);

  // TODO read SPECTRAL_WINDOW and make maps for id -> freqs and widths

  // Fill info with the data required to repopulate BDA_TIME_FACTOR
  t = ms_.keywordSet().asTable("BDA_TIME_FACTOR");
  std::vector<unsigned int> baseline_factors;
  baseline_factors.reserve(nbl_);
  for (unsigned int i = 0; i < t.nrow(); ++i) {
    baseline_factors.emplace_back(t.col("FACTOR").getDouble(i));
    // TODO add to freqs and widths
  }
  info().update(baseline_factors);

  // TODO
  // Fill info with the data required repopulate the main table
  // info().set(std::move(freqs), std::move(widths));
}

}  // namespace DPPP
}  // namespace DP3
