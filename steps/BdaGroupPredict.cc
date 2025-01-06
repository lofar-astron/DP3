// BdaGroupPredict.cc: DP3 step class that predicts BDA'd visibilities.
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Sebastiaan van der Tol

#include "BdaGroupPredict.h"
#include "Predict.h"
#include "ResultStep.h"

#include <iostream>

#include <dp3/base/DP3.h>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using dp3::base::BdaBuffer;
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

  /// @return The required fields for the sub-steps in this group.
  common::Fields GetRequiredFields() const {
    return base::GetChainRequiredFields(predict_step_);
  }

  /// Create the predict and result step for this group
  /// to be called after all baselines have been added
  void MakeSteps(base::DPInfo& info_in, const common::ParameterSet& parset,
                 std::string& prefix,
                 std::vector<std::string> source_patterns) {
    predict_step_ = std::make_shared<Predict>(parset, prefix, source_patterns);
    result_step_ = std::make_shared<ResultStep>();

    predict_step_->setNextStep(result_step_);

    DPInfo info(info_in);
    std::size_t nr_baselines = baselines_.size();
    std::vector<int> ant1(nr_baselines), ant2(nr_baselines);

    for (std::size_t bl_idx = 0; bl_idx < nr_baselines; ++bl_idx) {
      ant1[bl_idx] = info_in.getAnt1()[baselines_[bl_idx]];
      ant2[bl_idx] = info_in.getAnt2()[baselines_[bl_idx]];
    }

    write_back_info_.resize(nr_baselines);

    info.setAntennas(info_in.antennaNames(), info_in.antennaDiam(),
                     info_in.antennaPos(), ant1, ant2);

    std::vector<double> chanFreqs(info_in.chanFreqs(baselines_[0]));
    std::vector<double> chanWidths(info_in.chanWidths(baselines_[0]));
    std::size_t nr_chan = chanFreqs.size();
    info.setChannels(std::move(chanFreqs), std::move(chanWidths));
    predict_step_->setInfo(info);
    const std::array<std::size_t, 3> shape{nr_baselines, nr_chan,
                                           info_in.ncorr()};
    dpbuffer_ = std::make_unique<base::DPBuffer>();
    dpbuffer_->GetData().resize(shape);
    dpbuffer_->GetWeights().resize(shape);
    dpbuffer_->GetFlags().resize(shape);
    dpbuffer_->GetUvw().resize({nr_baselines, 3});
  }

  /// Write information on the predict step for this baseline group to the
  /// output stream os
  void Show(std::ostream& os) const {
    for (std::shared_ptr<Step> step = predict_step_; step;
         step = step->getNextStep()) {
      step->show(os);
    };
  }

  void ShowTimings(std::ostream& os, const double duration) const {
    for (std::shared_ptr<Step> step = predict_step_; step;
         step = step->getNextStep()) {
      step->showTimings(os, duration);
    }
  }

  /// Process one row if BDA data
  /// requests are buffered until the baseline group is complete
  void ProcessRow(const base::BdaBuffer::Row& row,
                  std::complex<float>* row_data, std::size_t& row_counter,
                  int bl_idx) {
    double time = row.time + row.interval / 2;

    // Check whether the time for this baseline is the same as the time for the
    // group
    if (nr_baselines_requested_ && (abs(dpbuffer_->GetTime() - time) > 1e-3)) {
      throw std::runtime_error(
          "Incomplete data: missing baselines in BDA buffer.");
    }

    // Copy data from BDA buffer row into the (regular) buffer for this baseline
    // group
    std::copy_n(row.uvw, 3, &dpbuffer_->GetUvw()(bl_idx, 0));

    write_back_info_[bl_idx] = {row_data, &row_counter};
    dpbuffer_->SetTime(time);
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
    predict_step_->process(std::move(dpbuffer_));

    // Get the result out of the Result step and reuse the DPBuffer.
    dpbuffer_ = result_step_->take();

    // Loop over all baselines in baselinegroup
    const std::size_t nr_baselines = baselines_.size();
    assert(dpbuffer_->GetData().shape(0) == nr_baselines);
    const std::size_t baseline_size =
        dpbuffer_->GetData().shape(1) * dpbuffer_->GetData().shape(2);
    for (std::size_t bl = 0; bl < nr_baselines; ++bl) {
      // Copy result from regular predict into corresponding BdaBuffer in queue
      std::copy_n(&dpbuffer_->GetData()(bl, 0, 0), baseline_size,
                  write_back_info_[bl].data);
      // Increment counter of rows filled in BdaBuffer
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
  std::unique_ptr<base::DPBuffer> dpbuffer_;
  struct WriteBackInfo {
    std::complex<float>* data;
    std::size_t* row_counter;
  };
  std::vector<WriteBackInfo> write_back_info_;
  std::size_t nr_baselines_requested_;
};

BdaGroupPredict::BdaGroupPredict(const common::ParameterSet& parset,
                                 const std::string& prefix)
    : parset_(parset), name_(prefix) {}

BdaGroupPredict::BdaGroupPredict(
    const common::ParameterSet& parset, const std::string& prefix,
    const std::vector<std::string>& source_patterns)
    : parset_(parset), name_(prefix), source_patterns_(source_patterns) {}

BdaGroupPredict::~BdaGroupPredict() {}

common::Fields BdaGroupPredict::getRequiredFields() const {
  common::Fields fields;
  for (const auto& entry : averaging_to_baseline_group_map_) {
    const BaselineGroup& blg = entry.second;
    fields |= blg.GetRequiredFields();
  }
  return fields;
}

void BdaGroupPredict::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);

  // Loop over all baselines, grouping them by averaging parameters
  for (std::size_t bl = 0; bl < info().nbaselines(); ++bl) {
    // Create a key describing the averaging in time and frequency
    auto averaging_key =
        std::make_pair(info().ntimeAvg(bl), info().chanFreqs(bl).size());
    // Get the baselinegroup for this averaging, if needed a new baselinegroup
    // will be created
    BaselineGroup& blg = averaging_to_baseline_group_map_[averaging_key];
    // Baseline will be added to the group, its index in the group will be the
    // current size of the baseline group
    std::size_t idx_in_blg = blg.GetSize();
    blg.AddBaseline(bl);
    index_to_baseline_group_map_.push_back(std::make_pair(&blg, idx_in_blg));
  }

  for (auto& entry : averaging_to_baseline_group_map_) {
    BaselineGroup& blg = entry.second;
    blg.MakeSteps(info(), parset_, name_, source_patterns_);
  }
}

base::Direction BdaGroupPredict::GetFirstDirection() const {
  if (index_to_baseline_group_map_.empty()) {
    throw std::runtime_error("BdaGroupPredict is not initialized");
  }
  return index_to_baseline_group_map_.front().first->GetFirstDirection();
}

void BdaGroupPredict::show(std::ostream& os) const {
  os << "BdaGroupPredict " << name_ << '\n';
  os << "Using a regular predict per baseline group. Baseline groups total: "
     << averaging_to_baseline_group_map_.size() << "\n";
  if (!averaging_to_baseline_group_map_.empty()) {
    os << "Predict for first baseline group\n";
    averaging_to_baseline_group_map_.begin()->second.Show(os);
  }
}

void BdaGroupPredict::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " BdaGroupPredict " << name_ << '\n';
  os << " Predict for first baseline group\n";
  averaging_to_baseline_group_map_.begin()->second.ShowTimings(os, duration);
}

bool BdaGroupPredict::process(std::unique_ptr<base::BdaBuffer> buffer) {
  timer_.start();

  buffers_.push({std::move(buffer), 0});

  std::size_t& row_counter = buffers_.back().nr_rows_filled;

  buffers_.back().buffer->AddData();  // Input buffer may not contain data.

  const std::vector<base::BdaBuffer::Row>& rows =
      buffers_.back().buffer->GetRows();
  for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
    const base::BdaBuffer::Row& row = rows[row_index];
    std::complex<float>* row_data = buffers_.back().buffer->GetData(row_index);

    BaselineGroup& blg = *index_to_baseline_group_map_[row.baseline_nr].first;
    int bl_in_group_idx = index_to_baseline_group_map_[row.baseline_nr].second;
    blg.ProcessRow(row, row_data, row_counter, bl_in_group_idx);
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
