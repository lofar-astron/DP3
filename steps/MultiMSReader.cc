// MultiMSReader.cc: DPPP step reading from multiple MSs
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "MultiMSReader.h"

#include <iostream>

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/version.h>
#include <casacore/casa/Utilities/GenSort.h>
#include <casacore/casa/OS/Conversion.h>

#include "../base/DPLogger.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

using casacore::Cube;
using casacore::IPosition;
using casacore::RefRows;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

MultiMSReader::MultiMSReader(const std::vector<std::string>& msNames,
                             const common::ParameterSet& parset,
                             const string& prefix)
    : itsFirst(-1), itsNMissing(0), itsMSNames(msNames) {
  if (msNames.empty())
    throw std::runtime_error("No names of MeasurementSets given");
  itsStartChanStr = parset.getString(prefix + "startchan", "0");
  itsNrChanStr = parset.getString(prefix + "nchan", "0");
  itsUseFlags = parset.getBool(prefix + "useflag", true);
  itsDataColName = parset.getString(prefix + "datacolumn", "DATA");
  itsFlagColName = parset.getString(prefix + "flagcolumn", "FLAG");
  itsWeightColName =
      parset.getString(prefix + "weightcolumn", "WEIGHT_SPECTRUM"),
  itsMissingData = parset.getBool(prefix + "missingdata", false);
  itsAutoWeight = parset.getBool(prefix + "autoweight", false);
  itsNeedSort = parset.getBool(prefix + "sort", false);
  itsOrderMS = parset.getBool(prefix + "orderms", true);
  // Open all MSs.
  itsReaders.reserve(msNames.size());
  for (const std::string& name : msNames) {
    if (!casacore::Table::isReadable(name)) {
      // Ignore if the MS is missing.
      itsReaders.push_back(nullptr);
      itsResults.push_back(nullptr);
      itsNMissing++;
    } else {
      const casacore::MeasurementSet ms(name,
                                        casacore::TableLock::AutoNoReadLocking);
      if (HasBda(ms)) {
        throw std::invalid_argument(name +
                                    " contains BDA data. DP3 does not support "
                                    "multiple input MS with BDA data.");
      }
      auto reader =
          std::make_shared<MSReader>(ms, parset, prefix, itsMissingData);
      // Add a result step for the reader.
      auto result = std::make_shared<ResultStep>();
      reader->setNextStep(result);
      itsReaders.push_back(std::move(reader));
      itsResults.push_back(std::move(result));
      if (itsFirst < 0) {
        itsFirst = itsReaders.size() - 1;
      }
    }
  }

  // TODO: check if frequencies are regular, insert some empy readers
  // if necessary

  if (itsFirst < 0)
    throw std::runtime_error("All input MeasurementSets do not exist");
}

MultiMSReader::~MultiMSReader() {}

std::string MultiMSReader::msName() const { return itsMSNames.front(); }

void MultiMSReader::setFieldsToRead(const dp3::common::Fields& fields) {
  InputStep::setFieldsToRead(fields);

  // Read all fields except UVW in the MSReaders.
  // Since the UVW values are equal for all MSReaders, read those once,
  // directly into the target buffer.
  dp3::common::Fields reader_fields;
  if (fields.Data()) reader_fields |= kDataField;
  if (fields.Flags()) reader_fields |= kFlagsField;
  if (fields.Weights()) reader_fields |= kWeightsField;
  for (std::shared_ptr<MSReader>& reader : itsReaders) {
    if (reader) {
      reader->setFieldsToRead(reader_fields);
    }
  }
}

void MultiMSReader::handleBands() {
  if (itsNMissing > 0) {
    fillBands();
    return;
  }
  if (itsOrderMS) {
    sortBands();
  }

  // Collect the channel info of all MSs.
  std::vector<double> chanFreqs(itsNrChan);
  std::vector<double> chanWidths(itsNrChan);
  std::vector<double> resolutions(itsNrChan);
  std::vector<double> effectiveBW(itsNrChan);
  unsigned int inx = 0;
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    unsigned int nchan = itsReaders[i]->getInfo().nchan();
    casacore::objcopy(chanFreqs.data() + inx,
                      itsReaders[i]->getInfo().chanFreqs().data(), nchan);
    casacore::objcopy(chanWidths.data() + inx,
                      itsReaders[i]->getInfo().chanWidths().data(), nchan);
    casacore::objcopy(resolutions.data() + inx,
                      itsReaders[i]->getInfo().resolutions().data(), nchan);
    casacore::objcopy(effectiveBW.data() + inx,
                      itsReaders[i]->getInfo().effectiveBW().data(), nchan);
    inx += nchan;
  }
  info().setChannels(std::move(chanFreqs), std::move(chanWidths),
                     std::move(resolutions), std::move(effectiveBW),
                     itsReaders[itsFirst]->getInfo().refFreq(),
                     itsReaders[itsFirst]->getInfo().spectralWindow());
}

void MultiMSReader::sortBands() {
  // Order the bands (MSs) in order of frequency.
  int nband = itsReaders.size();
  casacore::Vector<double> freqs(nband);
  for (int i = 0; i < nband; ++i) {
    freqs[i] = itsReaders[i]->getInfo().chanFreqs().data()[0];
  }
  casacore::Vector<common::rownr_t> index;

#if CASACORE_MAJOR_VERSION < 3 || \
    (CASACORE_MAJOR_VERSION == 3 && CASACORE_MINOR_VERSION < 4)
  casacore::GenSortIndirect<double>::sort(index, freqs);
#else
  casacore::GenSortIndirect<double, common::rownr_t>::sort(index, freqs);
#endif
  std::vector<std::shared_ptr<MSReader>> oldReaders(itsReaders);
  for (int i = 0; i < nband; ++i) {
    itsReaders[i] = oldReaders[index[i]];
  }
}

void MultiMSReader::fillBands() {
  if (itsOrderMS)
    throw std::runtime_error("Cannot order the MSs if some are missing");
  // Get channel width (which should be the same for all bands).
  double chanw = itsReaders[itsFirst]->getInfo().chanWidths().data()[0];
  // Get frequency for first subband.
  double freq = itsReaders[itsFirst]->getInfo().chanFreqs().data()[0];
  freq -= itsFirst * itsFillNChan * chanw;
  // Add missing channels to the total nr.
  itsNrChan += itsNMissing * itsFillNChan;
  // Collect the channel info of all MSs.
  std::vector<double> chanFreqs(itsNrChan);
  std::vector<double> chanWidths(itsNrChan);
  unsigned int inx = 0;
  // Data for a missing MS can only be inserted if all other MSs have
  // the same nr of channels and are in increasing order of freq.
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      if (itsReaders[i]->getInfo().nchan() != itsFillNChan)
        throw std::runtime_error(
            "An MS is missing; the others should have equal nchan");
      // Check if all channels have the same width and are consecutive.
      const std::vector<double>& freqs = itsReaders[i]->getInfo().chanFreqs();
      const std::vector<double>& width = itsReaders[i]->getInfo().chanWidths();
      if (freqs[0] < freq && !casacore::near(freqs[0], freq, 1e-5))
        throw std::runtime_error(
            "Subbands should be in increasing order of frequency; found " +
            std::to_string(freqs[0]) + ", expected " + std::to_string(freq) +
            " (diff=" + std::to_string(freqs[0] - freq) + ')');
      freq = freqs[itsFillNChan - 1] + width[itsFillNChan - 1];
      casacore::objcopy(chanFreqs.data() + inx, freqs.data(), itsFillNChan);
      casacore::objcopy(chanWidths.data() + inx, width.data(), itsFillNChan);
      inx += itsFillNChan;
    } else {
      // Insert channel info for missing MSs.
      for (unsigned int j = 0; j < itsFillNChan; ++j) {
        chanFreqs[inx] = freq;
        chanWidths[inx] = chanw;
        freq += chanw;
        inx++;
      }
    }
  }

  info().setChannels(std::move(chanFreqs), std::move(chanWidths));
}

bool MultiMSReader::process(std::unique_ptr<DPBuffer> buffer) {
  // Reduce memory allocation overhead by reusing the DPBuffer from the
  // previous iteration.
  std::unique_ptr<DPBuffer> recycled_buffer = itsResults[itsFirst]->extract();
  if (!recycled_buffer) recycled_buffer = std::make_unique<DPBuffer>();

  // Stop if at end.
  if (!itsReaders[itsFirst]->process(std::move(recycled_buffer))) {
    return false;  // end of input
  }
  const DPBuffer& buf1 = itsResults[itsFirst]->get();
  buffer->setTime(buf1.getTime());
  buffer->setExposure(buf1.getExposure());
  buffer->setRowNrs(buf1.getRowNrs());
  // Size the buffers if they should be read.
  if (getFieldsToRead().Data()) {
    buffer->ResizeData({itsNrBl, itsNrChan, itsNrCorr});
  }
  if (getFieldsToRead().Flags()) {
    buffer->ResizeFlags({itsNrBl, itsNrChan, itsNrCorr});
  }
  // Loop through all readers and get data and flags.
  IPosition start(3, 0, 0, 0);
  IPosition end(3, itsNrCorr - 1, 0, itsNrBl - 1);
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      if (int(i) != itsFirst) {
        // Reduce memory allocation overhead by reusing the DPBuffer from the
        // previous iteration.
        recycled_buffer = itsResults[i]->extract();
        if (!recycled_buffer) recycled_buffer = std::make_unique<DPBuffer>();
        itsReaders[i]->process(std::move(recycled_buffer));
      }
      const DPBuffer& msBuf = itsResults[i]->get();
      if (msBuf.getRowNrs().empty())
        throw std::runtime_error(
            "When using multiple MSs, the times in all MSs have to be "
            "consecutive; this is not the case for MS " +
            std::to_string(i));
      // Copy data and flags.
      end[1] = start[1] + itsReaders[i]->getInfo().nchan() - 1;
      if (getFieldsToRead().Data()) {
        buffer->GetCasacoreData()(start, end) = msBuf.GetCasacoreData();
      }
      if (getFieldsToRead().Flags()) {
        buffer->GetCasacoreFlags()(start, end) = msBuf.GetCasacoreFlags();
      }
    } else {
      end[1] = start[1] + itsFillNChan - 1;
      if (getFieldsToRead().Data()) {
        buffer->GetCasacoreData()(start, end) = casacore::Complex();
      }
      if (getFieldsToRead().Flags()) {
        buffer->GetCasacoreFlags()(start, end) = true;
      }
    }
    start[1] = end[1] + 1;
  }

  if (getFieldsToRead().Uvw()) {
    // All Measurement Sets have the same UVWs, so use the first one.
    itsReaders[itsFirst]->getUVW(buffer->getRowNrs(), buffer->getTime(),
                                 *buffer);
  }
  if (getFieldsToRead().Weights()) getWeights(*buffer);

  getNextStep()->process(std::move(buffer));
  return true;
}

void MultiMSReader::finish() {
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      itsReaders[i]->finish();
    }
  }
  getNextStep()->finish();
}

void MultiMSReader::updateInfo(const DPInfo& infoIn) {
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      itsReaders[i]->updateInfo(infoIn);
    }
  }
  info() = itsReaders[itsFirst]->getInfo();
  // Use the first valid MS as the standard MS (for meta data)
  // Get meta data and check they are equal for all MSs.
  itsMS = itsReaders[itsFirst]->table();
  itsFirstTime = getInfo().firstTime();
  itsLastTime = getInfo().lastTime();
  itsTimeInterval = getInfo().timeInterval();
  itsSelBL = itsReaders[itsFirst]->baselineSelection();
  itsNrCorr = getInfo().ncorr();
  itsNrBl = getInfo().nbaselines();
  itsNrChan = 0;
  itsFillNChan = getInfo().nchan();
  itsStartChan = itsReaders[itsFirst]->startChan();
  itsBaseRowNrs = itsReaders[itsFirst]->getBaseRowNrs();
  for (std::size_t i = 0; i < itsMSNames.size(); ++i) {
    if (itsReaders[i]) {
      const DPInfo& rdinfo = itsReaders[i]->getInfo();
      if (!casacore::near(getInfo().firstTime(), rdinfo.firstTime()))
        throw std::runtime_error("First time of MS " + itsMSNames[i] +
                                 " differs from " + itsMSNames[itsFirst]);
      if (!casacore::near(getInfo().lastTime(), rdinfo.lastTime()))
        throw std::runtime_error("Last time of MS " + itsMSNames[i] +
                                 " differs from " + itsMSNames[itsFirst]);
      if (!casacore::near(getInfo().timeInterval(), rdinfo.timeInterval()))
        throw std::runtime_error("Time interval of MS " + itsMSNames[i] +
                                 " differs from " + itsMSNames[itsFirst]);
      if (itsNrCorr != rdinfo.ncorr())
        throw std::runtime_error("Number of correlations of MS " +
                                 itsMSNames[i] + " differs from " +
                                 itsMSNames[itsFirst]);
      if (itsNrBl != rdinfo.nbaselines())
        throw std::runtime_error("Number of baselines of MS " + itsMSNames[i] +
                                 " differs from " + itsMSNames[itsFirst]);
      if (getInfo().antennaSet() != rdinfo.antennaSet())
        throw std::runtime_error("Antenna set of MS " + itsMSNames[i] +
                                 " differs from " + itsMSNames[itsFirst]);
      if (getInfo().getAnt1() != rdinfo.getAnt1())
        throw std::runtime_error("Baseline order (ant1) of MS " +
                                 itsMSNames[i] + " differs from " +
                                 itsMSNames[itsFirst]);
      if (getInfo().getAnt2() != rdinfo.getAnt2())
        throw std::runtime_error("Baseline order (ant2) of MS " +
                                 itsMSNames[i] + " differs from " +
                                 itsMSNames[itsFirst]);
      itsNrChan += rdinfo.nchan();
    }
  }
  // Handle the bands and take care of missing MSs.
  // Sort them if needed.
  handleBands();

  // Initialize the flag counters.
  itsFlagCounter.init(getInfo());
}

void MultiMSReader::show(std::ostream& os) const {
  os << "MultiMSReader" << '\n';
  os << "  input MSs:      " << itsMSNames[0] << '\n';
  for (std::size_t i = 1; i < itsMSNames.size(); ++i) {
    os << "                  " << itsMSNames[i] << '\n';
  }
  if (!itsSelBL.empty()) {
    os << "  baseline:       " << itsSelBL << '\n';
  }
  os << "  band            " << getInfo().spectralWindow() << '\n';
  os << "  startchan:      " << itsStartChan << "  (" << itsStartChanStr << ')'
     << '\n';
  os << "  nchan:          " << itsNrChan << "  (" << itsNrChanStr << ')';
  if (getInfo().channelsAreRegular()) {
    os << " (regularly spaced)" << '\n';
  } else {
    os << " (NOT regularly spaced)" << '\n';
  }
  os << "  ncorrelations:  " << itsNrCorr << '\n';
  os << "  nbaselines:     " << itsNrBl << '\n';
  os << "  ntimes:         " << itsMS.nrow() / itsNrBl << '\n';
  os << "  time interval:  " << itsTimeInterval << '\n';
  os << "  DATA column:    " << itsDataColName << '\n';
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      if (itsReaders[i]->missingData()) {
        os << "      column missing in  " << itsMSNames[i] << '\n';
      }
    } else {
      os << "      MS missing         " << itsMSNames[i] << '\n';
    }
  }
  os << "  WEIGHT column:  " << itsWeightColName << '\n';
  os << "  autoweight:     " << std::boolalpha << itsAutoWeight << '\n';
}

void MultiMSReader::showCounts(std::ostream& os) const {
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      itsReaders[i]->showCounts(os);
    }
  }
}

void MultiMSReader::showTimings(std::ostream& os, double duration) const {
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      itsReaders[i]->showTimings(os, duration);
    }
  }
}

void MultiMSReader::getWeights(DPBuffer& buffer) {
  const RefRows& rowNrs = buffer.getRowNrs();
  buffer.ResizeWeights({itsNrBl, itsNrChan, itsNrCorr});
  Cube<float>& weights = buffer.GetCasacoreWeights();
  IPosition start(3, 0, 0, 0);
  IPosition end(3, itsNrCorr - 1, 0, itsNrBl - 1);
  for (const std::shared_ptr<ResultStep>& result : itsResults) {
    if (result) {
      const unsigned int n_channels = result->get().GetWeights().shape(1);
      end[1] = start[1] + n_channels - 1;
      weights(start, end) = result->get().GetCasacoreWeights();
    } else {  // Use zero weights for dummy bands.
      end[1] = start[1] + itsFillNChan - 1;
      weights(start, end) = float(0);
    }
    start[1] = end[1] + 1;
  }
}

}  // namespace steps
}  // namespace dp3
