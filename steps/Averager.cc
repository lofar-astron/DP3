// Averager.cc: DP3 step class to average in time and/or freq
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Averager.h"

#include <iomanip>

#include <aocommon/parallelfor.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Utilities/Regex.h>
#include <boost/algorithm/string/trim.hpp>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"
#include "../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using casacore::Cube;
using casacore::IPosition;

namespace dp3 {
namespace steps {

const common::Fields Averager::kRequiredFields =
    kDataField | kFlagsField | kWeightsField | kFullResFlagsField | kUvwField;
const common::Fields Averager::kProvidedFields = kRequiredFields;

Averager::Averager(const common::ParameterSet& parset, const string& prefix)
    : itsName(prefix),
      itsMinNPoint(parset.getUint(prefix + "minpoints", 1)),
      itsMinPerc(parset.getFloat(prefix + "minperc", 0.) / 100.),
      itsNTimes(0),
      itsOriginalTimeInterval(0),
      itsNoAvg(true) {
  string freqResolutionStr = parset.getString(prefix + "freqresolution", "0");
  itsFreqResolution = getFreqHz(freqResolutionStr);

  if (itsFreqResolution > 0) {
    itsNChanAvg = 0;  // Will be set later in updateinfo
  } else {
    itsNChanAvg = parset.getUint(prefix + "freqstep", 1);
  }

  itsTimeResolution = parset.getFloat(prefix + "timeresolution", 0.);
  if (itsTimeResolution > 0) {
    itsNTimeAvg = 0;  // Will be set later in updateInfo
  } else {
    itsNTimeAvg = parset.getUint(prefix + "timestep", 1);
  }
}

Averager::Averager(const string& stepName, unsigned int nchanAvg,
                   unsigned int ntimeAvg)
    : itsName(stepName),
      itsFreqResolution(0),
      itsTimeResolution(0),
      itsNChanAvg(nchanAvg == 0 ? 1 : nchanAvg),
      itsNTimeAvg(ntimeAvg == 0 ? 1 : ntimeAvg),
      itsMinNPoint(1),
      itsMinPerc(0),
      itsNTimes(0),
      itsOriginalTimeInterval(0),
      itsNoAvg(itsNChanAvg == 1 && itsNTimeAvg == 1) {}

Averager::Averager(const string& stepName, double freq_resolution,
                   double time_resolution)
    : itsName(stepName),
      itsFreqResolution(freq_resolution),
      itsTimeResolution(time_resolution),
      itsNChanAvg(0),
      itsNTimeAvg(0),
      itsMinNPoint(1),
      itsMinPerc(0),
      itsNTimes(0),
      itsOriginalTimeInterval(0),
      itsNoAvg(itsNChanAvg == 1 && itsNTimeAvg == 1) {}

Averager::~Averager() {}

void Averager::updateInfo(const base::DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  info().setMetaChanged();

  if (itsNChanAvg <= 0) {
    if (itsFreqResolution > 0) {
      double chanwidth = infoIn.chanWidths()[0];
      itsNChanAvg = std::max(1, int(itsFreqResolution / chanwidth + 0.5));
    } else {
      itsNChanAvg = 1;
    }
  }

  itsOriginalTimeInterval = infoIn.timeInterval();
  if (itsNTimeAvg <= 0) {
    if (itsTimeResolution > 0) {
      itsNTimeAvg =
          std::max(1, int(itsTimeResolution / itsOriginalTimeInterval + 0.5));
    } else {
      itsNTimeAvg = 1;
    }
  }

  itsNoAvg = (itsNChanAvg == 1 && itsNTimeAvg == 1);

  // Adapt averaging to available nr of channels and times.
  itsNTimeAvg = std::min(itsNTimeAvg, infoIn.ntime());
  itsNChanAvg = info().update(itsNChanAvg, itsNTimeAvg);

  if (!itsNoAvg) loop_.SetNThreads(getInfo().nThreads());
}

void Averager::show(std::ostream& os) const {
  os << "Averager " << itsName << '\n';
  os << "  freqstep:       " << itsNChanAvg;
  if (itsFreqResolution > 0) {
    os << " (set by freqresolution: " << itsFreqResolution << " Hz)" << '\n';
  }
  os << "  timestep:       " << itsNTimeAvg;
  if (itsTimeResolution > 0) {
    os << " (set by timeresolution: " << itsTimeResolution << ")";
  }
  os << '\n';
  os << "  minpoints:      " << itsMinNPoint << '\n';
  os << "  minperc:        " << 100 * itsMinPerc << '\n';
}

void Averager::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " Averager " << itsName << '\n';
}

bool Averager::process(std::unique_ptr<base::DPBuffer> buffer) {
  assert(buffer->GetData().shape(1) == buffer->GetWeights().shape(1));

  // Nothing needs to be done if no averaging.
  if (itsNoAvg) {
    getNextStep()->process(std::move(buffer));
    return true;
  }
  itsTimer.start();
  // Sum the data in time applying the weights.
  // The summing in channel and the averaging is done in function average.
  if (itsNTimes == 0) {
    // The first time we move because that is faster than first clearing
    // and adding thereafter.
    assert(!itsBuf);
    itsBuf = std::move(buffer);
    // Ensure `itsBuf` does not refer to data in other DPBuffers.
    // Full res flags are not needed here, since they are resized below.
    itsBuf->MakeIndependent(kDataField | kWeightsField | kFlagsField |
                            kUvwField);
    itsNPoints.resize(itsBuf->GetData().shape());
    assert(itsBuf->GetData().shape(1) == itsBuf->GetWeights().shape(1));
    itsAvgAll = itsBuf->GetData() * itsBuf->GetWeights();
    itsWeightAll = itsBuf->GetWeights();
    // Take care of the fullRes flags.
    // We have to shape the output array and copy to a part of it.

    // Make sure the current fullResFlags are up to date with the flags
    itsBuf->MergeFullResFlags();

    // Extract fullResFlags before resizing the field in itsBuffer.
    // Until DPBuffer stores data in XTensor objects, copy the Casacore object
    // that contains the full resolution flags. After resizing them in itsBuf,
    // the Casacore object is the only owner of the old flags. 'full_res_flags'
    // is therefore valid as long as 'casa_full_res_flags' is in scope.
    casacore::Cube<bool> casa_full_res_flags =
        itsBuf->GetCasacoreFullResFlags();
    aocommon::xt::Span<bool, 3> full_res_flags = itsBuf->GetFullResFlags();

    const std::array<size_t, 3> shape = full_res_flags.shape();
    // More time entries, same chan and bl
    itsBuf->ResizeFullResFlags(shape[0], shape[1] * itsNTimeAvg, shape[2]);
    itsBuf->MakeIndependent(kFullResFlagsField);
    itsBuf->GetFullResFlags().fill(
        true);  // initialize for times missing at end
    copyFullResFlags(full_res_flags, itsBuf->GetFlags(), 0);

    // Set middle of new interval.
    const double time = itsBuf->getTime() + 0.5 * (getInfo().timeInterval() -
                                                   itsOriginalTimeInterval);
    itsBuf->setTime(time);
    itsBuf->setExposure(getInfo().timeInterval());
    // Only set.
    itsNPoints.fill(1);
    // Set flagged points to zero.
    auto flag_iterator = itsBuf->GetFlags().cbegin();
    auto data_iterator = itsBuf->GetData().begin();
    auto weight_iterator = itsBuf->GetWeights().begin();

    for (int& n : itsNPoints) {
      if (*flag_iterator) {
        // Flagged data point
        n = 0;
        *data_iterator = casacore::Complex();
        *weight_iterator = 0;
      } else {
        // Weigh the data point
        *data_iterator *= *weight_iterator;
      }
      ++flag_iterator;
      ++data_iterator;
      ++weight_iterator;
    }
  } else {
    // Not the first time.
    // For now we assume that all timeslots have the same nr of baselines,
    // so check if the buffer sizes are the same.
    assert(itsBuf);
    if (itsBuf->GetData().shape() != buffer->GetData().shape())
      throw std::runtime_error(
          "Inconsistent buffer sizes in Averager, possibly because of "
          "inconsistent nr of baselines in timeslots");
    itsBuf->GetUvw() += buffer->GetUvw();

    // Make sure the current fullResFlags are up to date with the flags
    buffer->MergeFullResFlags();
    copyFullResFlags(buffer->GetFullResFlags(), buffer->GetFlags(), itsNTimes);

    // Ignore flagged points.
    auto data_in_iterator = buffer->GetData().cbegin();
    auto weights_in_iterator = buffer->GetWeights().cbegin();
    auto flags_in_iterator = buffer->GetFlags().cbegin();
    auto data_out_iterator = itsBuf->GetData().begin();
    auto weights_out_iterator = itsBuf->GetWeights().begin();
    auto avg_all_iterator = itsAvgAll.begin();
    auto weight_all_iterator = itsWeightAll.begin();
    for (int& n : itsNPoints) {
      *avg_all_iterator += *data_in_iterator * *weights_in_iterator;
      *weight_all_iterator += *weights_in_iterator;
      if (!*flags_in_iterator) {
        *data_out_iterator += *data_in_iterator * *weights_in_iterator;
        *weights_out_iterator += *weights_in_iterator;
        ++n;
      }
      ++data_in_iterator;
      ++weights_in_iterator;
      ++flags_in_iterator;
      ++data_out_iterator;
      ++avg_all_iterator;
      ++weights_out_iterator;
      ++weight_all_iterator;
    }
  }
  // Do the averaging if enough time steps have been processed.
  itsNTimes += 1;
  if (itsNTimes >= itsNTimeAvg) {
    average();
    itsTimer.stop();
    getNextStep()->process(std::move(itsBuf));
    itsBuf.reset();
    itsNTimes = 0;
  } else {
    itsTimer.stop();
  }
  return true;
}

void Averager::finish() {
  // Average remaining entries.
  if (itsNTimes > 0) {
    itsTimer.start();
    average();
    itsTimer.stop();
    getNextStep()->process(std::move(itsBuf));
    itsNTimes = 0;
  }
  // Let the next steps finish.
  getNextStep()->finish();
}

void Averager::average() {
  // Extract the input data and weights before resizing itsBuf.
  // Until DPBuffer stores data in XTensor objects, Casacore owns the data.
  // Copying the Casacore objects ensures that the data remains valid, also
  // in 'data_in' and 'weights_in', which are views on that data.
  casacore::Cube<casacore::Complex> casa_data_in = itsBuf->GetCasacoreData();
  casacore::Cube<float> casa_weights_in = itsBuf->GetCasacoreWeights();
  aocommon::xt::Span<std::complex<float>, 3> data_in = itsBuf->GetData();
  aocommon::xt::Span<float, 3> weights_in = itsBuf->GetWeights();

  const unsigned int n_bl = data_in.shape(0);
  const unsigned int n_chan_in = data_in.shape(1);
  const unsigned int n_chan_out = (n_chan_in + itsNChanAvg - 1) / itsNChanAvg;
  const unsigned int n_corr = data_in.shape(2);
  itsBuf->ResizeData(n_bl, n_chan_out, n_corr);
  itsBuf->ResizeWeights(n_bl, n_chan_out, n_corr);
  itsBuf->ResizeFlags(n_bl, n_chan_out, n_corr);
  itsBuf->MakeIndependent(kDataField | kWeightsField);
  aocommon::xt::Span<std::complex<float>, 3>& data_out = itsBuf->GetData();
  aocommon::xt::Span<float, 3>& weights_out = itsBuf->GetWeights();
  aocommon::xt::Span<bool, 3>& flags_out = itsBuf->GetFlags();
  assert(data_in.data() != data_out.data());
  assert(weights_in.data() != weights_out.data());

  loop_.Run(0, n_bl, [&](size_t bl_begin, size_t bl_end) {
    for (unsigned int bl = bl_begin; bl < bl_end; ++bl) {
      for (unsigned int corr = 0; corr < n_corr; ++corr) {
        for (unsigned int chan_out = 0; chan_out < n_chan_out; ++chan_out) {
          const unsigned int chan_in_begin = chan_out * itsNChanAvg;
          const unsigned int n_averaged_chan =
              std::min(itsNChanAvg, n_chan_in - chan_in_begin);
          const unsigned int chan_in_end = chan_in_begin + n_averaged_chan;
          const unsigned int n_averaged_all = n_averaged_chan * itsNTimes;

          std::complex<float> sum_data;
          std::complex<float> sum_all_data;
          float sum_weights = 0;
          float sum_all_weights = 0;
          unsigned int sum_n_points = 0;

          for (unsigned int chan_in = chan_in_begin; chan_in < chan_in_end;
               ++chan_in) {
            // Note: weight is accounted for in process().
            sum_data += data_in(bl, chan_in, corr);
            sum_all_data += itsAvgAll(bl, chan_in, corr);
            sum_weights += weights_in(bl, chan_in, corr);
            sum_all_weights += itsWeightAll(bl, chan_in, corr);
            sum_n_points += itsNPoints(bl, chan_in, corr);
          }

          // Flag the point if insufficient unflagged data.
          if (sum_weights == 0 || sum_n_points < itsMinNPoint ||
              sum_n_points < n_averaged_all * itsMinPerc) {
            data_out(bl, chan_out, corr) =
                (sum_all_weights == 0 ? std::complex<float>()
                                      : sum_all_data / sum_all_weights);
            flags_out(bl, chan_out, corr) = true;
            weights_out(bl, chan_out, corr) = sum_all_weights;
          } else {
            data_out(bl, chan_out, corr) = sum_data / sum_weights;
            flags_out(bl, chan_out, corr) = false;
            weights_out(bl, chan_out, corr) = sum_weights;
          }
        }
      }
    }
  });
  // The result UVWs are the average of the input.
  // If ever needed, UVWCalculator can be used to calculate the UVWs.
  itsBuf->GetUvw() /= double(itsNTimes);
}

void Averager::copyFullResFlags(
    const aocommon::xt::Span<bool, 3>& full_res_flags,
    const aocommon::xt::Span<bool, 3>& flags, int time_index) {
  // Copy the fullRes flags to the given index.
  // Furthermore the appropriate FullRes flags are set for a
  // flagged data point. It can be the case that an input data point
  // has been averaged before, thus has fewer channels than FullResFlags.
  // nchan and nbl are the same for in and out.
  // ntimout is a multiple of ntimavg.
  const std::array<size_t, 3> shape_in = full_res_flags.shape();
  const std::array<size_t, 3> shape_out = itsBuf->GetFullResFlags().shape();
  const std::array<size_t, 3> shape_flags = flags.shape();
  // original nr of channels
  const unsigned int nchan = shape_in[2];
  // nr of averaged times in input data
  const unsigned int ntimavg = shape_in[1];
  // nr of averaged channels in input data
  const unsigned int nchanavg = nchan / shape_flags[1];
  // nr of averaged times in output data
  const unsigned int ntimout = shape_out[1];
  const unsigned int nbl = shape_in[0];       // nr of baselines
  const unsigned int ncorr = shape_flags[2];  // nr of correlations (in FLAG)
  // in has to be copied to the correct time index in out.
  bool* outBase =
      itsBuf->GetFullResFlags().data() + nchan * ntimavg * time_index;
  for (unsigned int k = 0; k < nbl; ++k) {
    const bool* inPtr = full_res_flags.data() + k * nchan * ntimavg;
    const bool* flagPtr = flags.data() + k * ncorr * shape_flags[1];
    bool* outPtr = outBase + k * nchan * ntimout;
    memcpy(outPtr, inPtr, nchan * ntimavg * sizeof(bool));
    // Applying the flags only needs to be done if the input data
    // was already averaged before.
    if (ntimavg > 1 || nchanavg > 1) {
      for (unsigned int j = 0; j < shape_flags[1]; ++j) {
        // If a data point is flagged, the flags in the corresponding
        // FullRes window have to be set.
        // Only look at the flags of the first correlation.
        if (*flagPtr) {
          bool* avgPtr = outPtr + j * nchanavg;
          for (unsigned int i = 0; i < ntimavg; ++i) {
            std::fill(avgPtr, avgPtr + nchanavg, true);
            avgPtr += nchan;
          }
        }
        flagPtr += ncorr;
      }
    }
  }
}

double Averager::getFreqHz(const string& freqstr) {
  casacore::String unit;
  // See if a unit is given at the end.
  casacore::String v(freqstr);
  // Remove possible trailing blanks.
  boost::algorithm::trim_right(v);
  casacore::Regex regex("[a-zA-Z]+$");
  casacore::String::size_type pos = v.index(regex);
  if (pos != casacore::String::npos) {
    unit = v.from(pos);
    v = v.before(pos);
  }
  // Set value and unit.

  double value = common::strToDouble(v);
  if (unit.empty()) {
    return value;
  } else {
    casacore::Quantity q(value, unit);
    return q.getValue("Hz", true);
  }
}

}  // namespace steps
}  // namespace dp3
