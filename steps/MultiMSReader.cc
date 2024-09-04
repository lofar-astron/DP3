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

  if (parset.getString(prefix + "startchan", "0") != "0")
    throw std::runtime_error(
        "startchan is not supported when reading multiple MeasurementSets.");

  if (parset.getString(prefix + "nchan", "0") != "0")
    throw std::runtime_error(
        "nchan is not supported when reading multiple MeasurementSets.");

  // Inherited members from MSReader must be initialized here.
  const bool allow_missing_data = parset.getBool(prefix + "missingdata", false);

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
          std::make_shared<MSReader>(ms, parset, prefix, allow_missing_data);
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

void MultiMSReader::ValidateBands() {
  for (const Reader& reader : readers_) {
    const std::shared_ptr<MSReader>& ms_reader = reader.ms_reader;
    if (ms_reader) {
      const DPInfo& reader_info = ms_reader->getInfo();
      const std::string& name = reader.name;
      const std::string& first_name = readers_.front().name;
      if (!casacore::near(getInfo().firstTime(), reader_info.firstTime()))
        throw std::runtime_error("First time of MS " + name + " differs from " +
                                 first_name);
      if (!casacore::near(getInfo().lastTime(), reader_info.lastTime()))
        throw std::runtime_error("Last time of MS " + name + " differs from " +
                                 first_name);
      if (!casacore::near(getInfo().timeInterval(), reader_info.timeInterval()))
        throw std::runtime_error("Time interval of MS " + name +
                                 " differs from " + first_name);
      if (getInfoOut().ncorr() != reader_info.ncorr())
        throw std::runtime_error("Number of correlations of MS " + name +
                                 " differs from " + first_name);
      if (getInfoOut().nbaselines() != reader_info.nbaselines())
        throw std::runtime_error("Number of baselines of MS " + name +
                                 " differs from " + first_name);
      if (getInfoOut().antennaSet() != reader_info.antennaSet())
        throw std::runtime_error("Antenna set of MS " + name +
                                 " differs from " + first_name);
      if (getInfoOut().getAnt1() != reader_info.getAnt1())
        throw std::runtime_error("Baseline order (ant1) of MS " + name +
                                 " differs from " + first_name);
      if (getInfoOut().getAnt2() != reader_info.getAnt2())
        throw std::runtime_error("Baseline order (ant2) of MS " + name +
                                 " differs from " + first_name);
    }
  }
}

void MultiMSReader::HandleBands() {
  // Collect the channel info of all MSs.
  unsigned int n_channels = 0;
  for (const Reader& reader : readers_) {
    n_channels += reader.ms_reader->getInfoOut().nchan();
  }

  std::vector<double> frequencies(n_channels);
  std::vector<double> widths(n_channels);
  std::vector<double> resolutions(n_channels);
  std::vector<double> effectiveBW(n_channels);
  unsigned int index = 0;
  for (const Reader& reader : readers_) {
    const DPInfo& reader_info = reader.ms_reader->getInfoOut();
    std::copy_n(reader_info.chanFreqs().data(), reader_info.nchan(),
                frequencies.data() + index);
    std::copy_n(reader_info.chanWidths().data(), reader_info.nchan(),
                widths.data() + index);
    std::copy_n(reader_info.resolutions().data(), reader_info.nchan(),
                resolutions.data() + index);
    std::copy_n(reader_info.effectiveBW().data(), reader_info.nchan(),
                effectiveBW.data() + index);
    index += reader_info.nchan();
  }

  const DPInfo first_reader_info = readers_.front().ms_reader->getInfoOut();
  info().setChannels(std::move(frequencies), std::move(widths),
                     std::move(resolutions), std::move(effectiveBW),
                     first_reader_info.refFreq(),
                     first_reader_info.spectralWindow());
}

void MultiMSReader::SortBands() {
  if (itsNMissing > 0) {
    throw std::runtime_error("Cannot order the MSs if some are missing");
  }

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

void MultiMSReader::FillBands() {
  unsigned int n_channels = 0;
  for (const Reader& reader : readers_) {
    if (reader.ms_reader) {
      n_channels += reader.ms_reader->getInfoOut().nchan();
    }
  }

  // Get channel width (which should be the same for all bands).
  double chanw = readers_[itsFirst].ms_reader->getInfo().chanWidths().data()[0];
  // Get frequency for first subband.
  double freq = readers_[itsFirst].ms_reader->getInfo().chanFreqs().data()[0];
  freq -= itsFirst * itsFillNChan * chanw;
  // Add missing channels to the total nr.
  n_channels += itsNMissing * itsFillNChan;
  // Collect the channel info of all MSs.
  std::vector<double> chanFreqs(n_channels);
  std::vector<double> chanWidths(n_channels);
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
  MSReader& first_ms_reader = *readers_[itsFirst].ms_reader;
  if (!first_ms_reader.process(std::move(recycled_buffer))) {
    return false;  // end of input
  }

  const std::vector<std::string>& extra_data_column_names =
      first_ms_reader.ExtraDataColumnNames();

  const DPBuffer& buf1 = readers_[itsFirst].result->get();
  buffer->SetTime(buf1.GetTime());
  buffer->SetExposure(buf1.GetExposure());
  buffer->SetRowNumbers(buf1.GetRowNumbers());
  // Size the buffers if they should be read.
  const std::array<size_t, 3> shape{getInfoOut().nbaselines(),
                                    getInfoOut().nchan(), getInfoOut().ncorr()};
  if (getFieldsToRead().Data()) {
    buffer->GetData().resize(shape);
    for (const std::string& column_name : extra_data_column_names) {
      buffer->AddData(column_name);
    }
  }
  if (getFieldsToRead().Flags()) {
    buffer->GetFlags().resize(shape);
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
        for (const std::string& column_name : extra_data_column_names) {
          xt::view(buffer->GetData(column_name), xt::all(), channel_range,
                   xt::all()) = msBuf.GetData(column_name);
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
        for (const std::string& column_name : extra_data_column_names) {
          xt::view(buffer->GetData(column_name), xt::all(), channel_range,
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

  const MSReader& first_reader = *readers_[itsFirst].ms_reader;
  Step::updateInfo(first_reader.getInfoOut());

  // Use the first valid MS as the standard MS (for meta data)
  // Get meta data and check they are equal for all MSs.
  itsFillNChan = getInfoOut().nchan();

  ValidateBands();

  if (itsOrderMS) {
    SortBands();
  }

  if (itsNMissing > 0) {
    FillBands();
  } else {
    HandleBands();
  }
}

void MultiMSReader::show(std::ostream& os) const {
  const MSReader& first_reader = *readers_[itsFirst].ms_reader;
  os << "MultiMSReader" << '\n';
  os << "  input MSs:      " << readers_.front().name << '\n';
  for (std::size_t i = 1; i < readers_.size(); ++i) {
    os << "                  " << readers_[i].name << '\n';
  }
  if (!first_reader.baselineSelection().empty()) {
    os << "  baseline:       " << first_reader.baselineSelection() << '\n';
  }
  os << "  band            " << getInfoOut().spectralWindow() << '\n';
  os << "  nchan:          " << getInfoOut().nchan();
  if (getInfoOut().channelsAreRegular()) {
    os << " (regularly spaced)" << '\n';
  } else {
    os << " (NOT regularly spaced)" << '\n';
  }
  os << "  ncorrelations:  " << getInfoOut().ncorr() << '\n';
  os << "  nbaselines:     " << getInfoOut().nbaselines() << '\n';
  os << "  ntimes:         " << getInfoOut().ntime() << '\n';
  os << "  time interval:  " << getInfoOut().timeInterval() << '\n';
  os << "  DATA column:    " << first_reader.DataColumnName() << '\n';
  for (const Reader& reader : readers_) {
    if (reader.ms_reader) {
      if (reader.ms_reader->MissingData()) {
        os << "      column missing in  " << reader.name << '\n';
      }
    } else {
      os << "      MS missing         " << reader.name << '\n';
    }
  }
  os << "  WEIGHT column:  " << first_reader.WeightColumnName() << '\n';
  os << "  autoweight:     " << std::boolalpha << first_reader.AutoWeight()
     << '\n';
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
  buffer->GetWeights().resize(
      {getInfoOut().nbaselines(), getInfoOut().nchan(), getInfoOut().ncorr()});
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
