// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Mick Veldhuis

#include "Transfer.h"

#include <xtensor/xview.hpp>

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Cube.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "../common/ParameterSet.h"
#include "../base/FlagCounter.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

Transfer::Transfer(const common::ParameterSet& parset,
                   const std::string& prefix)
    : name_(prefix),
      source_ms_path_(parset.getString(prefix + "source_ms")),
      source_data_column_(parset.getString(
          prefix + "datacolumn", casacore::MS::columnName(casacore::MS::DATA))),
      transfer_data_(parset.getBool(prefix + "data", false)),
      transfer_flags_(parset.getBool(prefix + "flags", false)) {
  if (!(transfer_data_ || transfer_flags_)) {
    throw std::runtime_error(
        "The type of data (visibilties and/or flags) to be transfered has not "
        "specified.");
  }

  // Initialise source MeasurementSet and its respective iterator
  ms_ = casacore::MeasurementSet(source_ms_path_,
                                 casacore::TableLock::AutoNoReadLocking);

  if (transfer_data_ && !ms_.tableDesc().isColumn(source_data_column_)) {
    throw std::runtime_error("The requested column (" + source_data_column_ +
                             ") does not exist in the source MS.");
  }

  ms_iterator_ = casacore::TableIterator(
      ms_, casacore::MS::columnName(casacore::MS::TIME));

  // Read time interval, channel width, baseline, channel, and correlation count
  // from source_ms
  time_interval_ = casacore::ScalarColumn<double>(ms_, "INTERVAL")(0);

  casacore::Table spw_table(ms_.keywordSet().asTable("SPECTRAL_WINDOW"));
  casacore::ArrayColumn<double> width_column(spw_table, "CHAN_WIDTH");
  const std::vector<double> channel_widths = width_column(0).tovector();

  casacore::Table table = ms_iterator_.table();
  const casacore::IPosition shape(
      casacore::ArrayColumn<casacore::Complex>(table, source_data_column_)
          .shape(0));
  const std::size_t n_correlations = shape[0];
  const std::size_t n_channels = shape[1];
  const std::size_t n_baselines = table.nrow();

  // Initialise time averaging factor, assigned in updateInfo()
  time_averaging_factor_ = 0;

  // Initialise channel bins using the central frequencies from CHAN_FREQ
  casacore::ArrayColumn<double> frequency_column(spw_table, "CHAN_FREQ");
  source_channel_upper_edges_ = frequency_column(0).tovector();
  for (std::size_t channel_index = 0;
       channel_index < source_channel_upper_edges_.size(); ++channel_index) {
    source_channel_upper_edges_[channel_index] +=
        0.5 * channel_widths[channel_index];
  }

  // Initialise data_ and flags_ to match the shape of what's stored in DPBuffer
  const std::array<std::size_t, 3> buffer_shape{n_baselines, n_channels,
                                                n_correlations};
  data_ = base::DPBuffer::DataType(buffer_shape);
  flags_ = base::DPBuffer::FlagsType(buffer_shape);

  // Create a filter substep and connect it to a result step
  filter_step_ = std::make_shared<Filter>(parset, prefix);
  result_step_ = std::make_shared<ResultStep>();
  filter_step_->setNextStep(result_step_);
}

void Transfer::ReadSourceMsVisibilities() {
  casacore::Table table = ms_iterator_.table();
  casacore::ArrayColumn<casacore::Complex> data_col(table, source_data_column_);

  const casacore::IPosition casa_shape{data_.shape()[2], data_.shape()[1],
                                       data_.shape()[0]};
  casacore::Cube<casacore::Complex> casa_data(casa_shape, data_.data(),
                                              casacore::SHARE);

  casa_data = data_col.getColumn();
}

void Transfer::ReadSourceMsFlags() {
  casacore::Table table = ms_iterator_.table();
  casacore::ArrayColumn<bool> flag_col(
      table, casacore::MS::columnName(casacore::MS::FLAG));

  const casacore::IPosition casa_shape{flags_.shape()[2], flags_.shape()[1],
                                       flags_.shape()[0]};
  casacore::Cube<bool> casa_flags(casa_shape, flags_.data(), casacore::SHARE);

  casa_flags = flag_col.getColumn();
}

template <typename T>
void Transfer::TransferSingleTimeSlot(
    T source, T& target, const std::vector<common::rownr_t>& target_row_numbers,
    const std::vector<double>& target_frequencies) const {
  const std::size_t n_source_channels = source.shape()[1];
  const std::size_t n_source_baselines = source.shape()[0];
  const std::size_t n_target_channels = target.shape()[1];
  const std::size_t n_target_baselines = target.shape()[0];
  if (n_source_channels == n_target_channels &&
      n_source_baselines == n_target_baselines) {
    target = source;
  } else if (n_source_channels == n_target_channels) {
    for (std::size_t source_row = 0; source_row < n_source_baselines;
         ++source_row) {
      const common::rownr_t target_row = target_row_numbers[source_row];
      xt::view(target, target_row, xt::all(), xt::all()) =
          xt::view(source, source_row, xt::all(), xt::all());
    }
  } else {
    std::size_t target_channel = 0;
    for (std::size_t channel_block = 0; channel_block < n_source_channels;
         ++channel_block) {
      const double source_channel_edge =
          source_channel_upper_edges_[channel_block];
      while (target_frequencies[target_channel] < source_channel_edge &&
             target_channel < n_target_channels) {
        for (std::size_t source_row = 0; source_row < n_source_baselines;
             ++source_row) {
          const common::rownr_t target_row = target_row_numbers[source_row];
          xt::view(target, target_row, target_channel, xt::all()) =
              xt::view(source, source_row, channel_block, xt::all());
        }

        ++target_channel;
      }
    }
  }
}

void Transfer::show(std::ostream& os) const {
  os << "Transfer " << name_ << '\n'
     << "  transfering: " << '\n'
     << "    data:  " << std::boolalpha << transfer_data_ << '\n'
     << "    flags: " << std::boolalpha << transfer_flags_ << '\n'
     << "  Source MS: " << source_ms_path_ << '\n'
     << "  time avg factor: " << time_averaging_factor_ << '\n'
     << "  source interval: " << time_interval_ << '\n'
     << "  target interval: " << getInfoOut().timeInterval() << '\n';
  filter_step_->show(os);
}

void Transfer::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " Transfer " << name_ << '\n';
}

bool Transfer::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  if (timestep_counter_ <= 0) {
    timestep_counter_ = time_averaging_factor_;
    if (!ms_iterator_.pastEnd()) {
      if (transfer_data_) ReadSourceMsVisibilities();
      if (transfer_flags_) ReadSourceMsFlags();
      ms_iterator_.next();
    }
  }

  common::Fields filter_fields = filter_step_->getRequiredFields();
  std::unique_ptr<DPBuffer> substep_buffer =
      std::make_unique<DPBuffer>(*buffer, filter_fields);
  filter_step_->process(std::move(substep_buffer));
  std::unique_ptr<DPBuffer> result_buffer = result_step_->take();

  std::vector<common::rownr_t> filtered_row_numbers =
      result_buffer->GetRowNumbers().tovector();

  const common::rownr_t first_row_number = buffer->GetRowNumbers()[0];
  for (auto& row_number : filtered_row_numbers) {
    row_number -= first_row_number;
  }

  if (filtered_row_numbers.size() != data_.shape()[0]) {
    throw std::runtime_error(
        "Transfer requires that the source and target MS have an equal "
        "number of baselines. Note: baselines might not have been properly "
        "filtered!");
  }

  const std::vector<double>& target_frequencies = getInfoOut().chanFreqs(0);

  // Fill data if DPBuffer is empty
  base::DPBuffer::DataType& data = buffer->GetData();
  if (data.size() == 0) {
    data.resize({getInfoOut().nbaselines(), getInfoOut().nchan(),
                 getInfoOut().ncorr()});
    data.fill(std::complex<float>(0, 0));
  }

  if (transfer_data_) {
    TransferSingleTimeSlot<base::DPBuffer::DataType>(
        data_, data, filtered_row_numbers, target_frequencies);
  }

  // Fill flags if DPBuffer is empty
  base::DPBuffer::FlagsType& flags = buffer->GetFlags();
  if (flags.size() == 0) {
    flags.resize({getInfoOut().nbaselines(), getInfoOut().nchan(),
                  getInfoOut().ncorr()});
    flags.fill(false);
  }

  if (transfer_flags_) {
    TransferSingleTimeSlot<base::DPBuffer::FlagsType>(
        flags_, flags, filtered_row_numbers, target_frequencies);
  }

  --timestep_counter_;

  timer_.stop();
  getNextStep()->process(std::move(buffer));
  return true;
}

void Transfer::updateInfo(const base::DPInfo& info_in) {
  Step::updateInfo(info_in);

  filter_step_->setInfo(info_in);

  if (getInfoOut().ncorr() != data_.shape()[2]) {
    throw std::runtime_error(
        "The transfer step requires that the source and target MS have an "
        "equal "
        "number of correlations");
  }

  // Check whether stations in source_ms exist in the the high-res target.
  casacore::Table source_antenna_table(source_ms_path_ + "/ANTENNA");
  casacore::ScalarColumn<casacore::String> source_antenna_name_column(
      source_antenna_table, "NAME");
  const std::size_t n_source_antennas = source_antenna_table.nrow();
  if (n_source_antennas != getInfoOut().nantenna()) {
    const std::vector<std::string>& antenna_names = getInfoOut().antennaNames();
    for (std::size_t row_number = 0; row_number < n_source_antennas;
         ++row_number) {
      const std::string source_antenna_name =
          source_antenna_name_column.get(row_number);
      const std::string antenna_name = antenna_names[row_number];
      if (source_antenna_name != antenna_name) {
        throw std::runtime_error(
            "Source MS antennas do not match with the target MS! Found '" +
            source_antenna_name + "' in source MS, but '" + antenna_name +
            "' in the target.");
      }
    }
  }

  // Determine the ratio of averaged time slots.
  const double target_time_interval = getInfoOut().timeInterval();
  const double averaging_factor = time_interval_ / target_time_interval;
  time_averaging_factor_ = std::round(averaging_factor);

  // Currently, Transfer only supports integer averaging factors.
  const double tolerance = 1.0e-5;
  if (std::abs(averaging_factor - time_averaging_factor_) > tolerance) {
    throw std::runtime_error("The time averaging factor is not integer in " +
                             name_);
  }

  const double start_time =
      casacore::ScalarColumn<double>(ms_, "TIME")(0) - 0.5 * time_interval_;
  const double difference = std::abs(getInfoOut().startTime() - start_time);
  if (difference > tolerance) {
    throw std::runtime_error(
        "Observation start times do not match! Transfer expects that the "
        "source MS is an averaged version of the target.");
  }
}

void Transfer::finish() { getNextStep()->finish(); }

}  // namespace steps
}  // namespace dp3
