// FlagTransfer.cc: DP3 step class to transfer flags
// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Mick Veldhuis

#include "FlagTransfer.h"

#include <xtensor/xview.hpp>

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Cube.h>

#include "base/DPBuffer.h"
#include "base/DPInfo.h"

#include <aocommon/logger.h>

#include "../common/ParameterSet.h"
#include "../base/FlagCounter.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

FlagTransfer::FlagTransfer(const common::ParameterSet& parset,
                           const std::string& prefix)
    : name_(prefix),
      source_ms_path_(parset.getString(prefix + "source_ms")),
      timestep_counter_(0) {
  aocommon::Logger::Warn
      << "FlagTransfer deprecation warning: please use the Transfer step with "
         "'transfer.flags=True' instead!\n";

  // Initialise source MeasurementSet and its respective iterator
  ms_ = casacore::MeasurementSet(source_ms_path_,
                                 casacore::TableLock::AutoNoReadLocking);

  if (!ms_.isColumn(casacore::MSMainEnums::PredefinedColumns::FLAG)) {
    throw std::runtime_error(
        "The FLAG column does not exist in the source MS.");
  }

  ms_iterator_ = casacore::TableIterator(
      ms_, casacore::MS::columnName(casacore::MS::TIME));

  // Read time interval, channel width, baseline, channel, and correlation count
  // from source_ms
  time_interval_ = casacore::ScalarColumn<double>(ms_, "INTERVAL")(0);

  casacore::Table spw_table(ms_.keywordSet().asTable("SPECTRAL_WINDOW"));
  casacore::ArrayColumn<double> width_column(spw_table, "CHAN_WIDTH");
  double channel_width = width_column(0).tovector().at(0);

  casacore::Table table = ms_iterator_.table();
  const casacore::IPosition shape(
      casacore::ArrayColumn<bool>(table,
                                  casacore::MS::columnName(casacore::MS::FLAG))
          .shape(0));
  const std::size_t n_correlations = shape[0];
  const std::size_t n_channels = shape[1];
  const std::size_t n_baselines = table.nrow();

  // Initialise time averaging factor, assigned in updateInfo()
  time_averaging_factor_ = 0;

  // Initialise channel bins using the central frequencies from CHAN_FREQ
  casacore::ArrayColumn<double> frequency_column(spw_table, "CHAN_FREQ");
  source_channel_upper_edges_ = frequency_column(0).tovector();
  for (double& bin : source_channel_upper_edges_) {
    bin += 0.5 * channel_width;
  }

  // Initialise flags_ to match the shape of the flags stored in DPBuffer
  flags_ = xt::xtensor<bool, 3>({n_baselines, n_channels, n_correlations});

  // Create a filter substep and connect it to a result step
  filter_step_ = std::make_shared<Filter>(parset, prefix);
  result_step_ = std::make_shared<ResultStep>();
  filter_step_->setNextStep(result_step_);
}

void FlagTransfer::ReadSourceMsFlags() {
  casacore::Table table = ms_iterator_.table();
  casacore::ArrayColumn<bool> flag_col(
      table, casacore::MS::columnName(casacore::MS::FLAG));

  const casacore::IPosition casa_shape(3, flags_.shape()[2], flags_.shape()[1],
                                       flags_.shape()[0]);
  casacore::Cube<bool> casa_flags(casa_shape, flags_.data(), casacore::SHARE);

  casa_flags = flag_col.getColumn();
}

void FlagTransfer::show(std::ostream& os) const {
  os << "FlagTransfer " << name_ << '\n'
     << " Source MS: " << source_ms_path_ << '\n'
     << " Time avg factor: " << time_averaging_factor_ << '\n'
     << " Source interval: " << time_interval_ << '\n'
     << " Target interval: " << getInfoOut().timeInterval() << '\n';
  filter_step_->show(os);
}

void FlagTransfer::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " FlagTransfer " << name_ << '\n';
}

bool FlagTransfer::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  if (timestep_counter_ <= 0) {
    timestep_counter_ = time_averaging_factor_;
    if (!ms_iterator_.pastEnd()) {
      ReadSourceMsFlags();
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

  if (filtered_row_numbers.size() != flags_.shape()[0]) {
    throw std::runtime_error(
        "FlagTransfer requires that the source and target MS have an equal "
        "number of baselines. Note: baselines might not have been properly "
        "filtered!");
  }

  // Fill flags if DPBuffer is empty
  base::DPBuffer::FlagsType& flags = buffer->GetFlags();
  if (flags.size() == 0) {
    flags.resize({getInfoOut().nbaselines(), getInfoOut().nchan(),
                  getInfoOut().ncorr()});
    flags.fill(false);
  }

  const std::size_t n_source_channels = flags_.shape()[1];
  const std::size_t n_source_baselines = flags_.shape()[0];
  if (n_source_channels == getInfoOut().nchan() &&
      n_source_baselines == getInfoOut().nbaselines()) {
    flags = flags_;
  } else if (n_source_channels == getInfoOut().nchan()) {
    for (std::size_t source_row = 0; source_row < n_source_baselines;
         ++source_row) {
      common::rownr_t target_row = filtered_row_numbers[source_row];
      xt::view(flags, target_row, xt::all(), xt::all()) =
          xt::view(flags_, source_row, xt::all(), xt::all());
    }
  } else {
    std::size_t target_channel = 0;
    const std::vector<double> target_frequencies = getInfoOut().chanFreqs(0);
    for (std::size_t channel_block = 0; channel_block < n_source_channels;
         ++channel_block) {
      const double source_channel_edge =
          source_channel_upper_edges_[channel_block];
      while (target_frequencies[target_channel] < source_channel_edge) {
        if (target_channel >= getInfoOut().nchan()) {
          break;
        }

        for (std::size_t source_row = 0; source_row < n_source_baselines;
             ++source_row) {
          common::rownr_t target_row = filtered_row_numbers[source_row];
          xt::view(flags, target_row, target_channel, xt::all()) =
              xt::view(flags_, source_row, channel_block, xt::all());
        }

        ++target_channel;
      }
    }
  }

  --timestep_counter_;

  timer_.stop();
  getNextStep()->process(std::move(buffer));
  return true;
}

void FlagTransfer::updateInfo(const base::DPInfo& info_in) {
  Step::updateInfo(info_in);

  filter_step_->setInfo(info_in);

  if (getInfoOut().ncorr() != flags_.shape()[2]) {
    throw std::runtime_error(
        "FlagTransfer requires that the source and target MS have an equal "
        "number of correlations");
  }

  // Determine the ratio of averaged time slots
  const double target_time_interval = getInfoOut().timeInterval();
  const double averaging_factor = time_interval_ / target_time_interval;
  time_averaging_factor_ = std::round(averaging_factor);

  // Currently, FlagTransfer only supports integer averaging factors.
  const double tolerance = 1.0e-5;
  if (std::abs(averaging_factor - time_averaging_factor_) > tolerance) {
    throw std::runtime_error("The time averaging factor is not integer in " +
                             name_);
  }

  double start_time =
      casacore::ScalarColumn<double>(ms_, "TIME")(0) - 0.5 * time_interval_;
  double difference = std::abs(getInfoOut().startTime() - start_time);
  if (difference > tolerance) {
    throw std::runtime_error(
        "Observation start times do not match! FlagTransfer expects that the "
        "source MS is an averaged version of the target.");
  }
}

void FlagTransfer::finish() { getNextStep()->finish(); }

}  // namespace steps
}  // namespace dp3