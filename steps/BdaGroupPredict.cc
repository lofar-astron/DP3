// BdaGroupPredict.cc: DP3 step class that predicts BDA'd visibilities.
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Sebastiaan van der Tol

#include "BdaGroupPredict.h"
#include "Predict.h"

#include <iostream>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using dp3::base::BDABuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

/// class representing a group of baselines that have the same averaging / data
/// shape
class BdaGroupPredict::BaselineGroup {
 public:
  /// Add a baseline given by its index bl to this group
  void AddBaseline(std::size_t bl) { baselines_.push_back(bl); }

  /// Get the number of baselines in this group
  std::size_t GetSize() const { return baselines_.size(); }

  /// Create the predict and result step for this group
  /// to be called afer all baselines have been added
  void MakeSteps(InputStep &input, base::DPInfo &info_in,
                 const common::ParameterSet &parset, std::string &prefix,
                 std::vector<std::string> source_patterns) {
    predict_step_ =
        std::make_shared<Predict>(input, parset, prefix, source_patterns);
    result_step_ = std::make_shared<ResultStep>();

    predict_step_->setNextStep(result_step_);

    DPInfo info(info_in);
    std::size_t nr_baselines = baselines_.size();
    casacore::Vector<casacore::Int> ant1(nr_baselines), ant2(nr_baselines);

    for (std::size_t bl_idx = 0; bl_idx < nr_baselines; ++bl_idx) {
      ant1(bl_idx) = info_in.getAnt1()[baselines_[bl_idx]];
      ant2(bl_idx) = info_in.getAnt2()[baselines_[bl_idx]];
    }

    write_back_info_.resize(nr_baselines);

    info.set(info_in.antennaNames(), info_in.antennaDiam(),
             info_in.antennaPos(), ant1, ant2);

    std::vector<double> chanFreqs(info_in.chanFreqs(baselines_[0]));
    std::vector<double> chanWidths(info_in.chanWidths(baselines_[0]));
    std::size_t nr_chan = chanFreqs.size();
    info.set(std::move(chanFreqs),
             std::move(chanWidths));  // This does not update info.nchan() !!
    info.setNChan(nr_chan);           // So we need to set it
    predict_step_->setInfo(info);
    dpbuffer_.getData().resize(info_in.ncorr(), nr_chan, nr_baselines);
    dpbuffer_.getWeights().resize(info_in.ncorr(), nr_chan, nr_baselines);
    dpbuffer_.getFlags().resize(info_in.ncorr(), nr_chan, nr_baselines);
    dpbuffer_.getUVW().resize(3, nr_baselines);
  }

  /// Write information on the predict step for this baseline group to the
  /// output stream os
  void Show(std::ostream &os) const {
    Step::ShPtr step = predict_step_;
    do {
      step->show(os);
    } while (step = step->getNextStep());
  }

  void ShowTimings(std::ostream &os, const double duration) const {
    for (Step::ShPtr step = predict_step_; step; step = step->getNextStep()) {
      step->showTimings(os, duration);
    }
  }

  /// Process one row if BDA data
  /// requests are buffered until the baseline group is complete
  void ProcessRow(const base::BDABuffer::Row &row, std::size_t &row_counter,
                  int bl_idx) {
    double time = row.time + row.interval / 2;

    // Check whether the time for this baseline is the same as the time for the
    // group
    if (nr_baselines_requested_ && (abs(dpbuffer_.getTime() - time) > 1e-3)) {
      throw std::runtime_error(
          "Incomplete data: missing baselines in BDA buffer.");
    }

    // Copy data from BDA buffer row into the (regular) buffer for this baseline
    // group
    std::copy(row.uvw, row.uvw + 3, dpbuffer_.getUVW()[bl_idx].begin());

    write_back_info_[bl_idx] = {row.data, &row_counter};
    dpbuffer_.setTime(time);
    std::size_t nr_baselines = baselines_.size();

    // Flush if the baseline group is complete
    if (++nr_baselines_requested_ == nr_baselines) {
      Flush();
    }
  }

  // Execute the predict step, processing all baselines in group at once
  // Results are copied back the the BDA Buffer
  void Flush() {
    // Send data to the OnePredict step
    predict_step_->process(dpbuffer_);

    // Get the result out of the Result step
    base::DPBuffer &buf_out = result_step_->get();

    // Loop over all baselines in baselinegroup
    std::size_t nr_baselines = baselines_.size();
    for (std::size_t bl = 0; bl < nr_baselines; ++bl) {
      // Copy result from regular predict into corresponding BDABuffer in queue
      std::copy(buf_out.getData()[bl].begin(), buf_out.getData()[bl].end(),
                write_back_info_[bl].data);
      // Increment counter of rows filled in BDABuffer
      (*write_back_info_[bl].row_counter)++;
    }
    // Reset the number of baselines requested
    nr_baselines_requested_ = 0;
  }

  base::Direction GetFirstDirection() const {
    return predict_step_->GetFirstDirection();
  }

 private:
  std::vector<std::size_t> baselines_;
  std::shared_ptr<Predict> predict_step_;
  std::shared_ptr<ResultStep> result_step_;
  base::DPBuffer dpbuffer_;
  struct WriteBackInfo {
    std::complex<float> *data;
    std::size_t *row_counter;
  };
  std::vector<WriteBackInfo> write_back_info_;
  std::size_t nr_baselines_requested_;
};

BdaGroupPredict::BdaGroupPredict(InputStep &input,
                                 const common::ParameterSet &parset,
                                 const string &prefix)
    : input_(input), parset_(parset), name_(prefix) {}

BdaGroupPredict::BdaGroupPredict(
    InputStep &input, const common::ParameterSet &parset, const string &prefix,
    const std::vector<std::string> &source_patterns)
    : input_(input),
      parset_(parset),
      name_(prefix),
      source_patterns_(source_patterns) {}

BdaGroupPredict::~BdaGroupPredict() {}

void BdaGroupPredict::updateInfo(const DPInfo &infoIn) {
  Step::updateInfo(infoIn);
  info().setNeedVisData();
  info().setWriteData();

  // Loop over all baselines, grouping them by averaging parameters
  for (std::size_t bl = 0; bl < info().nbaselines(); ++bl) {
    // Create a key describing the averaging in time and frequency
    auto averaging_key =
        std::make_pair(info().ntimeAvg(bl), info().chanFreqs(bl).size());
    // Get the baselinegroup for this averaging, if needed a new baselinegroup
    // will be created
    BaselineGroup &blg = averaging_to_baseline_group_map_[averaging_key];
    // Baseline will be added to the group, its index in the group will be the
    // current size of the baseline group
    std::size_t idx_in_blg = blg.GetSize();
    blg.AddBaseline(bl);
    index_to_baseline_group_map_.push_back(std::make_pair(&blg, idx_in_blg));
  }

  for (auto &entry : averaging_to_baseline_group_map_) {
    BaselineGroup &blg = entry.second;
    blg.MakeSteps(input_, info(), parset_, name_, source_patterns_);
  }
}

base::Direction BdaGroupPredict::GetFirstDirection() const {
  if (index_to_baseline_group_map_.empty()) {
    throw std::runtime_error("BdaGroupPredict is not initialized");
  }
  return index_to_baseline_group_map_.front().first->GetFirstDirection();
}

void BdaGroupPredict::show(std::ostream &os) const {
  os << "BdaGroupPredict " << name_ << '\n';
  os << "Using a regular predict per baseline group. Baseline groups total: "
     << averaging_to_baseline_group_map_.size() << "\n";
  if (!averaging_to_baseline_group_map_.empty()) {
    os << "Predict for first baseline group\n";
    averaging_to_baseline_group_map_.begin()->second.Show(os);
  }
}

void BdaGroupPredict::showTimings(std::ostream &os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " BdaGroupPredict " << name_ << '\n';
  os << " Predict for first baseline group\n";
  averaging_to_baseline_group_map_.begin()->second.ShowTimings(os, duration);
}

bool BdaGroupPredict::process(std::unique_ptr<base::BDABuffer> buffer) {
  timer_.start();

  buffers_.push({std::move(buffer), 0});

  std::size_t &row_counter = buffers_.back().nr_rows_filled;

  const std::vector<base::BDABuffer::Row> &rows =
      buffers_.back().buffer->GetRows();
  for (auto const &row : rows) {
    BaselineGroup &blg = *index_to_baseline_group_map_[row.baseline_nr].first;
    int bl_in_group_idx = index_to_baseline_group_map_[row.baseline_nr].second;
    blg.ProcessRow(row, row_counter, bl_in_group_idx);
  }

  timer_.stop();
  while (!buffers_.empty() && buffers_.front().buffer->GetRows().size() ==
                                  buffers_.front().nr_rows_filled) {
    getNextStep()->process(std::move(buffers_.front().buffer));
    buffers_.pop();
  }
  return false;
}

void BdaGroupPredict::finish() {
  // Let the next steps finish.
  if (!buffers_.empty()) {
    throw std::runtime_error(
        "Incomplete data: missing baselines in BdaGroupPredict::finish().");
  }
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
