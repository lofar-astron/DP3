// MultiMSReader.cc: DP3 step reading from multiple MSs
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

#include <xtensor/xview.hpp>

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

MultiMSReader::MultiMSReader(const std::vector<std::string>& msNames,
                             const common::ParameterSet& parset,
                             const std::string& prefix)
    : itsOrderMS(parset.getBool(prefix + "orderms", true)),
      itsFirst(-1),
      itsNMissing(0),
      itsFillNChan(0) {
  if (msNames.empty())
    throw std::runtime_error("No names of MeasurementSets given");
  // Inherited members from MSReader must be initialized here.
  itsStartChanStr = parset.getString(prefix + "startchan", "0");
  itsNrChanStr = parset.getString(prefix + "nchan", "0");
  itsUseFlags = parset.getBool(prefix + "useflag", true);
  itsDataColName = parset.getString(prefix + "datacolumn", "DATA");
  itsExtraDataColNames = parset.getStringVector(prefix + "extradatacolumns",
                                                std::vector<std::string>());
  itsFlagColName = parset.getString(prefix + "flagcolumn", "FLAG");
  itsWeightColName =
      parset.getString(prefix + "weightcolumn", "WEIGHT_SPECTRUM"),
  itsMissingData = parset.getBool(prefix + "missingdata", false);
  itsAutoWeight = parset.getBool(prefix + "autoweight", false);
  itsNeedSort = parset.getBool(prefix + "sort", false);

  // Open all MSs.
  readers_.reserve(msNames.size());
  for (const std::string& name : msNames) {
    Reader reader;
    reader.name = name;
    if (!casacore::Table::isReadable(name)) {
      // Ignore if the MS is missing.
      itsNMissing++;
    } else {
      const casacore::MeasurementSet ms(name,
                                        casacore::TableLock::AutoNoReadLocking);
      if (HasBda(ms)) {
        throw std::invalid_argument(name +
                                    " contains BDA data. DP3 does not support "
                                    "multiple input MSs with BDA data.");
      }
      reader.ms_reader =
          std::make_shared<MSReader>(ms, parset, prefix, itsMissingData);
      // Add a result step for the reader.
      reader.result = std::make_shared<ResultStep>();
      reader.ms_reader->setNextStep(reader.result);
      if (itsFirst < 0) {
        itsFirst = readers_.size();
      }
    }
    readers_.emplace_back(std::move(reader));
  }

  // TODO: check if frequencies are regular, insert some empty readers
  // if necessary

  if (itsFirst < 0)
    throw std::runtime_error("All input MeasurementSets do not exist");
}

MultiMSReader::~MultiMSReader() {}

std::string MultiMSReader::msName() const { return readers_.front().name; }

void MultiMSReader::setFieldsToRead(const dp3::common::Fields& fields) {
  InputStep::setFieldsToRead(fields);

  // Read all fields except UVW in the MSReaders.
  // Since the UVW values are equal for all MSReaders, read those once,
  // directly into the target buffer.
  dp3::common::Fields reader_fields;
  if (fields.Data()) reader_fields |= kDataField;
  if (fields.Flags()) reader_fields |= kFlagsField;
  if (fields.Weights()) reader_fields |= kWeightsField;
  for (Reader& reader : readers_) {
    if (reader.ms_reader) {
      reader.ms_reader->setFieldsToRead(reader_fields);
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
  for (unsigned int i = 0; i < readers_.size(); ++i) {
    unsigned int nchan = readers_[i].ms_reader->getInfo().nchan();
    casacore::objcopy(chanFreqs.data() + inx,
                      readers_[i].ms_reader->getInfo().chanFreqs().data(),
                      nchan);
    casacore::objcopy(chanWidths.data() + inx,
                      readers_[i].ms_reader->getInfo().chanWidths().data(),
                      nchan);
    casacore::objcopy(resolutions.data() + inx,
                      readers_[i].ms_reader->getInfo().resolutions().data(),
                      nchan);
    casacore::objcopy(effectiveBW.data() + inx,
                      readers_[i].ms_reader->getInfo().effectiveBW().data(),
                      nchan);
    inx += nchan;
  }
  info().setChannels(std::move(chanFreqs), std::move(chanWidths),
                     std::move(resolutions), std::move(effectiveBW),
                     readers_[itsFirst].ms_reader->getInfo().refFreq(),
                     readers_[itsFirst].ms_reader->getInfo().spectralWindow());
}

void MultiMSReader::sortBands() {
  // Order the bands (MSs) in order of frequency.
  int nband = readers_.size();
  casacore::Vector<double> freqs(nband);
  for (int i = 0; i < nband; ++i) {
    freqs[i] = readers_[i].ms_reader->getInfo().chanFreqs().data()[0];
  }
  casacore::Vector<common::rownr_t> index;

#if CASACORE_MAJOR_VERSION < 3 || \
    (CASACORE_MAJOR_VERSION == 3 && CASACORE_MINOR_VERSION < 4)
  casacore::GenSortIndirect<double>::sort(index, freqs);
#else
  casacore::GenSortIndirect<double, common::rownr_t>::sort(index, freqs);
#endif
  std::vector<Reader> oldReaders(readers_);
  for (int i = 0; i < nband; ++i) {
    readers_[i] = std::move(oldReaders[index[i]]);
  }
}

void MultiMSReader::fillBands() {
  if (itsOrderMS)
    throw std::runtime_error("Cannot order the MSs if some are missing");
  // Get channel width (which should be the same for all bands).
  double chanw = readers_[itsFirst].ms_reader->getInfo().chanWidths().data()[0];
  // Get frequency for first subband.
  double freq = readers_[itsFirst].ms_reader->getInfo().chanFreqs().data()[0];
  freq -= itsFirst * itsFillNChan * chanw;
  // Add missing channels to the total nr.
  itsNrChan += itsNMissing * itsFillNChan;
  // Collect the channel info of all MSs.
  std::vector<double> chanFreqs(itsNrChan);
  std::vector<double> chanWidths(itsNrChan);
  unsigned int inx = 0;
  // Data for a missing MS can only be inserted if all other MSs have
  // the same nr of channels and are in increasing order of freq.
  for (Reader& reader : readers_) {
    std::shared_ptr<MSReader>& ms_reader = reader.ms_reader;
    if (ms_reader) {
      if (ms_reader->getInfo().nchan() != itsFillNChan)
        throw std::runtime_error(
            "An MS is missing; the others should have equal nchan");
      // Check if all channels have the same width and are consecutive.
      const std::vector<double>& freqs = ms_reader->getInfo().chanFreqs();
      const std::vector<double>& width = ms_reader->getInfo().chanWidths();
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
  // previous process() call.
  std::unique_ptr<DPBuffer> recycled_buffer = readers_[itsFirst].result->take();
  if (!recycled_buffer) recycled_buffer = std::make_unique<DPBuffer>();

  // Read first MS, stop if at end.
  if (!readers_[itsFirst].ms_reader->process(std::move(recycled_buffer))) {
    return false;  // end of input
  }
  const DPBuffer& buf1 = readers_[itsFirst].result->get();
  buffer->SetTime(buf1.GetTime());
  buffer->SetExposure(buf1.GetExposure());
  buffer->SetRowNumbers(buf1.GetRowNumbers());
  // Size the buffers if they should be read.
  if (getFieldsToRead().Data()) {
    buffer->GetData().resize({itsNrBl, itsNrChan, itsNrCorr});
    for (std::string columnName : itsExtraDataColNames) {
      buffer->AddData(columnName);
      buffer->GetData(columnName).resize({itsNrBl, itsNrChan, itsNrCorr});
    }
  }
  if (getFieldsToRead().Flags()) {
    buffer->GetFlags().resize({itsNrBl, itsNrChan, itsNrCorr});
  }

  // Loop through all readers and get data and flags.
  int first_channel = 0;
  int last_channel = 0;
  for (unsigned int i = 0; i < readers_.size(); ++i) {
    std::shared_ptr<MSReader>& ms_reader = readers_[i].ms_reader;
    if (ms_reader) {
      if (int(i) != itsFirst) {
        // Reduce memory allocation overhead by reusing the DPBuffer from the
        // previous process() call.
        recycled_buffer = readers_[i].result->take();
        if (!recycled_buffer) recycled_buffer = std::make_unique<DPBuffer>();
        ms_reader->process(std::move(recycled_buffer));
      }
      const DPBuffer& msBuf = readers_[i].result->get();
      if (msBuf.GetRowNumbers().empty())
        throw std::runtime_error(
            "When using multiple MSs, the times in all MSs have to be "
            "consecutive; this is not the case for MS " +
            std::to_string(i) + ": " + ms_reader->msName());
      // Copy data and flags.
      last_channel = first_channel + ms_reader->getInfo().nchan();
      auto channel_range = xt::range(first_channel, last_channel);
      if (getFieldsToRead().Data()) {
        xt::view(buffer->GetData(), xt::all(), channel_range, xt::all()) =
            msBuf.GetData();
        for (std::string columnName : itsExtraDataColNames) {
          xt::view(buffer->GetData(columnName), xt::all(), channel_range,
                   xt::all()) = msBuf.GetData(columnName);
        }
      }
      if (getFieldsToRead().Flags()) {
        xt::view(buffer->GetFlags(), xt::all(), channel_range, xt::all()) =
            msBuf.GetFlags();
      }
    } else {
      // Corresponding MS is missing.
      last_channel = first_channel + itsFillNChan;
      auto channel_range = xt::range(first_channel, last_channel);
      if (getFieldsToRead().Data()) {
        xt::view(buffer->GetData(), xt::all(), channel_range, xt::all())
            .fill(std::complex<float>());
        for (std::string columnName : itsExtraDataColNames) {
          xt::view(buffer->GetData(columnName), xt::all(), channel_range,
                   xt::all())
              .fill(std::complex<float>());
        }
      }
      if (getFieldsToRead().Flags()) {
        xt::view(buffer->GetFlags(), xt::all(), channel_range, xt::all())
            .fill(true);
      }
    }
    first_channel = last_channel;
  }

  if (getFieldsToRead().Uvw()) {
    // All Measurement Sets have the same UVWs, so use the first one.
    readers_[itsFirst].ms_reader->getUVW(buffer->GetRowNumbers(),
                                         buffer->GetTime(), *buffer);
  }
  if (getFieldsToRead().Weights()) getWeights(buffer);

  getNextStep()->process(std::move(buffer));
  return true;
}

void MultiMSReader::finish() {
  for (const Reader& reader : readers_) {
    if (reader.ms_reader) {
      reader.ms_reader->finish();
    }
  }
  getNextStep()->finish();
}

void MultiMSReader::updateInfo(const DPInfo& infoIn) {
  for (const Reader& reader : readers_) {
    if (reader.ms_reader) {
      reader.ms_reader->updateInfo(infoIn);
    }
  }
  std::shared_ptr<MSReader>& first_reader = readers_[itsFirst].ms_reader;
  info() = first_reader->getInfo();
  // Use the first valid MS as the standard MS (for meta data)
  // Get meta data and check they are equal for all MSs.
  itsMS = first_reader->table();
  itsFirstTime = getInfo().firstTime();
  itsMaximumTime = getInfo().lastTime();
  itsTimeInterval = getInfo().timeInterval();
  itsSelBL = first_reader->baselineSelection();
  itsNrCorr = getInfo().ncorr();
  itsNrBl = getInfo().nbaselines();
  itsNrChan = 0;
  itsFillNChan = getInfo().nchan();
  itsStartChan = first_reader->startChan();
  itsBaseRowNrs = first_reader->getBaseRowNrs();
  for (const Reader& reader : readers_) {
    const std::shared_ptr<MSReader>& ms_reader = reader.ms_reader;
    if (ms_reader) {
      const DPInfo& rdinfo = ms_reader->getInfo();
      const std::string& name = reader.name;
      const std::string& first_name = readers_.front().name;
      if (!casacore::near(getInfo().firstTime(), rdinfo.firstTime()))
        throw std::runtime_error("First time of MS " + name + " differs from " +
                                 first_name);
      if (!casacore::near(getInfo().lastTime(), rdinfo.lastTime()))
        throw std::runtime_error("Last time of MS " + name + " differs from " +
                                 first_name);
      if (!casacore::near(getInfo().timeInterval(), rdinfo.timeInterval()))
        throw std::runtime_error("Time interval of MS " + name +
                                 " differs from " + first_name);
      if (itsNrCorr != rdinfo.ncorr())
        throw std::runtime_error("Number of correlations of MS " + name +
                                 " differs from " + first_name);
      if (itsNrBl != rdinfo.nbaselines())
        throw std::runtime_error("Number of baselines of MS " + name +
                                 " differs from " + first_name);
      if (getInfo().antennaSet() != rdinfo.antennaSet())
        throw std::runtime_error("Antenna set of MS " + name +
                                 " differs from " + first_name);
      if (getInfo().getAnt1() != rdinfo.getAnt1())
        throw std::runtime_error("Baseline order (ant1) of MS " + name +
                                 " differs from " + first_name);
      if (getInfo().getAnt2() != rdinfo.getAnt2())
        throw std::runtime_error("Baseline order (ant2) of MS " + name +
                                 " differs from " + first_name);
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
  os << "  input MSs:      " << readers_.front().name << '\n';
  for (std::size_t i = 1; i < readers_.size(); ++i) {
    os << "                  " << readers_[i].name << '\n';
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
  for (const Reader& reader : readers_) {
    if (reader.ms_reader) {
      if (reader.ms_reader->missingData()) {
        os << "      column missing in  " << reader.name << '\n';
      }
    } else {
      os << "      MS missing         " << reader.name << '\n';
    }
  }
  os << "  WEIGHT column:  " << itsWeightColName << '\n';
  os << "  autoweight:     " << std::boolalpha << itsAutoWeight << '\n';
}

void MultiMSReader::showCounts(std::ostream& os) const {
  for (const Reader& reader : readers_) {
    if (reader.ms_reader) {
      reader.ms_reader->showCounts(os);
    }
  }
}

void MultiMSReader::showTimings(std::ostream& os, double duration) const {
  for (const Reader& reader : readers_) {
    if (reader.ms_reader) {
      reader.ms_reader->showTimings(os, duration);
    }
  }
}

void MultiMSReader::getWeights(std::unique_ptr<base::DPBuffer>& buffer) {
  buffer->GetWeights().resize({itsNrBl, itsNrChan, itsNrCorr});
  int first_channel = 0;
  int last_channel = 0;
  for (Reader& reader : readers_) {
    const std::shared_ptr<ResultStep>& result = reader.result;
    if (result) {
      last_channel = first_channel + result->get().GetWeights().shape(1);
      xt::view(buffer->GetWeights(), xt::all(),
               xt::range(first_channel, last_channel), xt::all()) =
          result->get().GetWeights();
    } else {
      last_channel = first_channel + itsFillNChan;
      xt::view(buffer->GetWeights(), xt::all(),
               xt::range(first_channel, last_channel), xt::all())
          .fill(0.0f);
    }
    first_channel = last_channel;
  }
}

}  // namespace steps
}  // namespace dp3
