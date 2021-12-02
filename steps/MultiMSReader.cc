// MultiMSReader.cc: DPPP step reading from multiple MSs
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "MultiMSReader.h"

#include "../base/DPBuffer.h"
#include "../base/DPLogger.h"
#include "../base/DPInfo.h"
#include "../base/Exceptions.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

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

#include <iostream>

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
    : itsFirst(-1),
      itsNMissing(0),
      itsMSNames(msNames),
      itsRegularChannels(true) {
  if (msNames.empty()) throw Exception("No names of MeasurementSets given");
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
  auto nullStep = std::make_shared<NullStep>();
  itsReaders.reserve(msNames.size());
  for (const std::string& name : msNames) {
    if (!casacore::Table::isReadable(name)) {
      // Ignore if the MS is missing.
      itsReaders.push_back(nullptr);
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
      // Add a null step for the reader.
      reader->setNextStep(nullStep);
      itsReaders.push_back(std::move(reader));
      if (itsFirst < 0) {
        itsFirst = itsReaders.size() - 1;
      }
    }
  }

  // TODO: check if frequencies are regular, insert some empy readers
  // if necessary

  if (itsFirst < 0) throw Exception("All input MeasurementSets do not exist");
  itsBuffers.resize(itsReaders.size());
}

MultiMSReader::~MultiMSReader() {}

void MultiMSReader::setReadVisData(bool readVisData) {
  itsReadVisData = readVisData || itsReadVisData;
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      itsReaders[i]->setReadVisData(readVisData);
    }
  }
}

std::string MultiMSReader::msName() const { return itsMSNames.front(); }

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
  info().set(std::move(chanFreqs), std::move(chanWidths),
             std::move(resolutions), std::move(effectiveBW));
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
  if (itsOrderMS) throw Exception("Cannot order the MSs if some are missing");
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
        throw Exception("An MS is missing; the others should have equal nchan");
      // Check if all channels have the same width and are consecutive.
      const casacore::Vector<double>& freqs =
          itsReaders[i]->getInfo().chanFreqs();
      const casacore::Vector<double>& width =
          itsReaders[i]->getInfo().chanWidths();
      if (freqs[0] < freq && !casacore::near(freqs[0], freq, 1e-5))
        throw Exception(
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

  info().set(std::move(chanFreqs), std::move(chanWidths));
}

bool MultiMSReader::process(const DPBuffer& buf) {
  // Stop if at end.
  if (!itsReaders[itsFirst]->process(buf)) {
    return false;  // end of input
  }
  const DPBuffer& buf1 = itsReaders[itsFirst]->getBuffer();
  itsBuffer.setTime(buf1.getTime());
  itsBuffer.setExposure(buf1.getExposure());
  itsBuffer.setRowNrs(buf1.getRowNrs());
  // Size the buffers.
  if (itsBuffer.getFlags().empty()) {
    if (itsReadVisData) {
      itsBuffer.getData().resize(IPosition(3, itsNrCorr, itsNrChan, itsNrBl));
    }
    itsBuffer.getFlags().resize(IPosition(3, itsNrCorr, itsNrChan, itsNrBl));
  }
  // Loop through all readers and get data and flags.
  IPosition s(3, 0, 0, 0);
  IPosition e(3, itsNrCorr - 1, 0, itsNrBl - 1);
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      if (int(i) != itsFirst) {
        itsReaders[i]->process(buf);
      }
      const DPBuffer& msBuf = itsReaders[i]->getBuffer();
      if (msBuf.getRowNrs().empty())
        throw Exception(
            "When using multiple MSs, the times in all MSs have to be "
            "consecutive; this is not the case for MS " +
            std::to_string(i));
      // Copy data and flags.
      e[1] = s[1] + itsReaders[i]->getInfo().nchan() - 1;
      if (itsReadVisData) {
        itsBuffer.getData()(s, e) = msBuf.getData();
      }
      itsBuffer.getFlags()(s, e) = msBuf.getFlags();
    } else {
      e[1] = s[1] + itsFillNChan - 1;
      if (itsReadVisData) {
        itsBuffer.getData()(s, e) = casacore::Complex();
      }
      itsBuffer.getFlags()(s, e) = true;
    }
    s[1] = e[1] + 1;
  }
  getNextStep()->process(itsBuffer);
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
  itsStartTime = getInfo().startTime();
  itsFirstTime = itsReaders[itsFirst]->firstTime();
  itsLastTime = itsReaders[itsFirst]->lastTime();
  itsTimeInterval = getInfo().timeInterval();
  itsSelBL = itsReaders[itsFirst]->baselineSelection();
  itsSpw = itsReaders[itsFirst]->spectralWindow();
  itsNrCorr = getInfo().ncorr();
  itsNrBl = getInfo().nbaselines();
  itsNrChan = 0;
  itsFillNChan = getInfo().nchan();
  itsStartChan = itsReaders[itsFirst]->startChan();
  itsFullResNChanAvg = itsReaders[itsFirst]->nchanAvgFullRes();
  itsFullResNTimeAvg = itsReaders[itsFirst]->ntimeAvgFullRes();
  itsHasFullResFlags = itsReaders[itsFirst]->hasFullResFlags();
  itsBaseRowNrs = itsReaders[itsFirst]->getBaseRowNrs();
  for (std::size_t i = 0; i < itsMSNames.size(); ++i) {
    if (itsReaders[i]) {
      const DPInfo& rdinfo = itsReaders[i]->getInfo();
      if (!casacore::near(itsStartTime, rdinfo.startTime()))
        throw Exception("Start time of MS " + itsMSNames[i] + " differs from " +
                        itsMSNames[itsFirst]);
      if (!casacore::near(itsLastTime, itsReaders[i]->lastTime()))
        throw Exception("Last time of MS " + itsMSNames[i] + " differs from " +
                        itsMSNames[itsFirst]);
      if (!casacore::near(itsTimeInterval, rdinfo.timeInterval()))
        throw Exception("Time interval of MS " + itsMSNames[i] +
                        " differs from " + itsMSNames[itsFirst]);
      if (itsNrCorr != rdinfo.ncorr())
        throw Exception("Number of correlations of MS " + itsMSNames[i] +
                        " differs from " + itsMSNames[itsFirst]);
      if (itsNrBl != rdinfo.nbaselines())
        throw Exception("Number of baselines of MS " + itsMSNames[i] +
                        " differs from " + itsMSNames[itsFirst]);
      if (itsFullResNChanAvg != itsReaders[i]->nchanAvgFullRes())
        throw Exception("FullResNChanAvg of MS " + itsMSNames[i] +
                        " differs from " + itsMSNames[itsFirst]);
      if (itsFullResNTimeAvg != itsReaders[i]->ntimeAvgFullRes())
        throw Exception("FullResNTimeAvg of MS " + itsMSNames[i] +
                        " differs from " + itsMSNames[itsFirst]);
      if (getInfo().antennaSet() != rdinfo.antennaSet())
        throw Exception("Antenna set of MS " + itsMSNames[i] +
                        " differs from " + itsMSNames[itsFirst]);
      if (getInfo().getAnt1() != rdinfo.getAnt1())
        throw Exception("Baseline order (ant1) of MS " + itsMSNames[i] +
                        " differs from " + itsMSNames[itsFirst]);
      if (getInfo().getAnt2() != rdinfo.getAnt2())
        throw Exception("Baseline order (ant2) of MS " + itsMSNames[i] +
                        " differs from " + itsMSNames[itsFirst]);
      itsNrChan += rdinfo.nchan();
      itsHasFullResFlags =
          (itsHasFullResFlags && itsReaders[i]->hasFullResFlags());
    }
  }
  // Handle the bands and take care of missing MSs.
  // Sort them if needed.
  handleBands();

  // check that channels are regularly spaced, give warning otherwise
  itsRegularChannels = info().channelsAreRegular();

  // Set correct nr of channels.
  info().setNChan(itsNrChan);
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
  os << "  band            " << itsSpw << '\n';
  os << "  startchan:      " << itsStartChan << "  (" << itsStartChanStr << ')'
     << '\n';
  os << "  nchan:          " << itsNrChan << "  (" << itsNrChanStr << ')';
  if (itsRegularChannels) {
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

void MultiMSReader::getUVW(const RefRows& rowNrs, double time, DPBuffer& buf) {
  // All MSs have the same UVWs, so use first one.
  itsReaders[itsFirst]->getUVW(rowNrs, time, buf);
}

void MultiMSReader::getWeights(const RefRows& rowNrs, DPBuffer& buf) {
  Cube<float>& weights = buf.getWeights();
  // Resize if needed (probably when called for first time).
  if (weights.empty()) {
    weights.resize(itsNrCorr, itsNrChan, itsNrBl);
  }
  IPosition s(3, 0, 0, 0);
  IPosition e(3, itsNrCorr - 1, 0, itsNrBl - 1);
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      unsigned int nchan = itsReaders[i]->getInfo().nchan();
      e[1] = s[1] + nchan - 1;
      itsReaders[i]->getWeights(rowNrs, itsBuffers[i]);
      weights(s, e) = itsBuffers[i].getWeights();
    } else {
      e[1] = s[1] + itsFillNChan - 1;
      weights(s, e) = float(0);
    }
    s[1] = e[1] + 1;
  }
}

bool MultiMSReader::getFullResFlags(const RefRows& rowNrs, DPBuffer& buf) {
  Cube<bool>& flags = buf.getFullResFlags();
  // Resize if needed (probably when called for first time).
  if (flags.empty()) {
    int norigchan = itsNrChan * itsFullResNChanAvg;
    flags.resize(norigchan, itsFullResNTimeAvg, itsNrBl);
  }
  // Return false if no fullRes flags available.
  if (!itsHasFullResFlags) {
    flags = false;
    return false;
  }
  // Flag everything if data rows are missing.
  if (rowNrs.rowVector().empty()) {
    flags = true;
    return true;
  }
  // Get the flags from all MSs and combine them.
  IPosition s(3, 0);
  IPosition e(flags.shape() - 1);
  for (unsigned int i = 0; i < itsReaders.size(); ++i) {
    if (itsReaders[i]) {
      itsReaders[i]->getFullResFlags(rowNrs, itsBuffers[i]);
      e[0] = s[0] + itsBuffers[i].getFullResFlags().shape()[0] - 1;
      flags(s, e) = itsBuffers[i].getFullResFlags();
    } else {
      e[0] = s[0] + itsFillNChan - 1;
      flags(s, e) = true;
    }
    s[0] = e[0] + 1;
  }
  return true;
}

}  // namespace steps
}  // namespace dp3
