// Copyright (C) 2022 ASRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ImageStep.h"

#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/views/xview.hpp>
#include <xtensor/views/xindex_view.hpp>

#include <aocommon/banddata.h>
#include <wsclean/wsclean.h>

#include "base/FlagCounter.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3::steps {

ImageStep::ImageStep(const common::ParameterSet& parset,
                     const std::string& prefix)
    : name_(prefix),
      cadence_(parset.getDouble(prefix + "cadence")),
      wsclean_options_(parset.getString(prefix + "options")),
      fits_prefix_(parset.getString(prefix + "fitsprefix", "")) {}

void ImageStep::updateInfo(const DPInfo& info_in) {
  Step::updateInfo(info_in);

  const double time_interval = getInfoOut().timeInterval();
  if (cadence_ < time_interval) {
    throw std::runtime_error(
        "The imaging cadence (" + std::to_string(cadence_) +
        ") should be larger than the time resolution of the input data (" +
        std::to_string(time_interval) + ")");
  }
  n_timesteps_per_image_ = std::round(cadence_ / time_interval);

  has_frequency_bda_ = getInfoOut().hasBDAChannels();
  has_time_bda_ = getInfoOut().ntimeAvgs().size() > 1;
}

void ImageStep::show(std::ostream& os) const {
  os << "ImageStep " << name_ << '\n';
  os << "  cadence:                 " << cadence_ << '\n';
  os << "  timesteps per image:     " << n_timesteps_per_image_ << '\n';
  if (!fits_prefix_.empty()) {
    os << "  FITS prefix:             " << fits_prefix_ << '\n';
  }
  os << "  wsclean options:         " << wsclean_options_ << '\n';
  os << "  in-memory MS:            " << '\n';
  os << "    has time BDA:          " << std::boolalpha << has_time_bda_
     << '\n';
  os << "    has frequency BDA:     " << std::boolalpha << has_frequency_bda_
     << '\n';
  os << "    observer:              " << getInfoOut().GetObserver() << '\n';
  os << "    telescope:             " << getInfoOut().GetTelescopeName()
     << '\n';
  os << "    field:                 " << getInfoOut().GetFieldName() << '\n';
}

void ImageStep::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " ImageStep " << name_ << '\n';
}

bool ImageStep::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  if (timestep_counter_ <= 0) {
    if (in_memory_ms_) {
      in_memory_ms_.reset();
    }

    InitializeInMemoryMs();
  }

  DPBuffer::WeightsType weights = buffer->GetWeights();
  xt::filtration(weights, buffer->GetFlags()) = 0.0f;

  const double current_time = buffer->GetTime();
  AddCurrentBufferToInMemoryMs(buffer->GetData(), weights, buffer->GetUvw(),
                               current_time);

  ++timestep_counter_;

  if (timestep_counter_ >= n_timesteps_per_image_) {
    std::vector<wsclean::InMemoryMs> in_memory_ms_list;
    in_memory_ms_list.push_back(std::move(*in_memory_ms_));

    const std::string prefix =
        fits_prefix_.empty() ? "" : FormatImageNamePrefix();
    wsclean::Image(wsclean_options_, std::move(in_memory_ms_list), prefix);

    timestep_counter_ = 0;

    ++image_counter_;
  }

  timer_.stop();
  getNextStep()->process(std::move(buffer));
  return false;
}

bool ImageStep::process(std::unique_ptr<base::BdaBuffer> buffer) {
  throw std::runtime_error(
      "Processing BDA data has not yet been implemented for the ImageStep");
  return getNextStep()->process(std::move(buffer));
}

void ImageStep::finish() { getNextStep()->finish(); }

void ImageStep::InitializeInMemoryMs() {
  assert(!has_frequency_bda_);

  in_memory_ms_ = std::make_unique<wsclean::InMemoryMs>();

  in_memory_ms_->has_frequency_bda = has_frequency_bda_;
  in_memory_ms_->has_time_bda = has_time_bda_;

  in_memory_ms_->observer = getInfoOut().GetObserver();
  in_memory_ms_->telescope_name = getInfoOut().GetTelescopeName();
  in_memory_ms_->field_name = getInfoOut().GetFieldName();

  in_memory_ms_->start_time = getInfoOut().startTime();
  in_memory_ms_->interval = getInfoOut().timeInterval();

  const base::Direction& phase_centre = getInfoOut().phaseCenterDirection();
  in_memory_ms_->phase_centre_ra = phase_centre.ra;
  in_memory_ms_->phase_centre_dec = phase_centre.dec;

  in_memory_ms_->antenna_names = getInfoOut().GetUsedAntennaNames();

  const std::set<aocommon::PolarizationEnum>& polarizations =
      getInfoOut().polarizations();
  in_memory_ms_->polarizations = std::vector<aocommon::PolarizationEnum>(
      polarizations.begin(), polarizations.end());

  const size_t n_channels = getInfoOut().nchan();
  const std::vector<double>& channel_frequencies = getInfoOut().chanFreqs();
  const std::vector<double>& channel_widths = getInfoOut().chanWidths();

  std::vector<aocommon::ChannelInfo> channels;
  channels.reserve(n_channels);
  for (size_t channel_index = 0; channel_index < n_channels; ++channel_index) {
    channels.emplace_back(channel_frequencies[channel_index],
                          channel_widths[channel_index]);
  }

  const double reference_frequency = getInfoOut().refFreq();
  aocommon::BandData band(channels, reference_frequency);
  in_memory_ms_->bands.AddBand(band);
}

void ImageStep::AddCurrentBufferToInMemoryMs(
    const DPBuffer::DataType& data, const DPBuffer::WeightsType& weights,
    const DPBuffer::UvwType& uvws, const double time) {
  const size_t n_baselines = data.shape(0);
  const std::vector<int>& antenna1 = getInfoOut().getAnt1();
  const std::vector<int>& antenna2 = getInfoOut().getAnt2();
  for (size_t baseline_index = 0; baseline_index < n_baselines;
       ++baseline_index) {
    wsclean::InMemoryRow& row = in_memory_ms_->rows.emplace_back();
    row.antenna1 = antenna1[baseline_index];
    row.antenna2 = antenna2[baseline_index];
    row.time = time;

    row.data_desc_id = 0;  // Always zero when not using BDA
    row.field_id = 0;      // Always zero for LOFAR

    const size_t n_frequencies = data.shape(1);
    const size_t n_polarizations = data.shape(2);
    const size_t row_data_size = n_frequencies * n_polarizations;
    row.data.resize(row_data_size);
    row.weights.resize(row_data_size);

    std::copy_n(&data(baseline_index, 0, 0), row_data_size, row.data.data());
    std::copy_n(&weights(baseline_index, 0, 0), row_data_size,
                row.weights.data());

    xt::adapt(row.uvw) = xt::view(uvws, baseline_index, xt::all());
  }
}

std::string ImageStep::FormatImageNamePrefix() const {
  std::ostringstream counter;

  // 5 digits are used, to be able to image 8 hour with 1 second chunks
  // (which would give 28800 images) without breaking the ordering of
  // the files.
  counter << std::setw(5) << std::setfill('0') << image_counter_;

  return fits_prefix_ + "-" + counter.str();
}

}  // namespace dp3::steps
