// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BdaDdeCal.h"

#include "../base/DPInfo.h"
#include "BDAResultStep.h"
#include "BdaPredict.h"

#include "../common/ParameterValue.h"
#include "../common/StreamUtil.h"

#include <boost/make_unique.hpp>

using dp3::base::BDABuffer;
using dp3::ddecal::BDASolverBuffer;
using dp3::steps::BdaPredict;
using dp3::common::operator<<;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

BdaDdeCal::BdaDdeCal(InputStep* input, const common::ParameterSet& parset,
                     const std::string& prefix)
    : settings_(parset, prefix),
      steps_(),
      result_steps_(),
      solver_buffer_(settings_.directions.size(), info().startTime(),
                     info().timeInterval()) {
  if (settings_.directions.empty()) {
    throw std::invalid_argument(
        "Invalid info in input parset: direction(s) must be specified");
  }
  InitializePredictSteps(input, parset, prefix);
}

void BdaDdeCal::InitializePredictSteps(InputStep* input,
                                       const common::ParameterSet& parset,
                                       const string& prefix) {
  for (const std::string& direction : settings_.directions) {
    const std::vector<std::string> source_patterns =
        common::ParameterValue(direction).getStringVector();

    steps_.push_back(
        std::make_shared<BdaPredict>(*input, parset, prefix, source_patterns));
    result_steps_.push_back(std::make_shared<BDAResultStep>());
    steps_.back()->setNextStep(result_steps_.back());
  }
}

void BdaDdeCal::updateInfo(const DPInfo& _info) {
  Step::updateInfo(_info);
  for (unsigned int i = 0; i < settings_.directions.size(); i++) {
    steps_[i]->setInfo(_info);
  }
}

bool BdaDdeCal::process(std::unique_ptr<base::BDABuffer> buffer) {
  const size_t data_size = buffer->GetNumberOfElements();

  BDABuffer::Fields fields(false);
  fields.data = true;
  fields.weights = true;
  const BDABuffer::Fields copyfields(false);

  // Feed metadata-only copies of the buffer to the steps.
  for (std::shared_ptr<ModelDataStep>& step : steps_) {
    step->process(boost::make_unique<BDABuffer>(*buffer, fields, copyfields));
  }

  // Combine the results from all steps.
  std::vector<std::unique_ptr<BDABuffer>> model_buffers;
  model_buffers.reserve(result_steps_.size());
  for (std::shared_ptr<BDAResultStep>& result_step : result_steps_) {
    std::vector<std::unique_ptr<BDABuffer>> results = result_step->Extract();

    // Each step should produce a single output buffer with the same shape as
    // the input buffer.
    assert(results.size() == 1);
    assert(results.front()->GetNumberOfElements() == data_size);
    model_buffers.push_back(std::move(results.front()));
  }

  if (settings_.only_predict) {
    // Add all model_buffers and use that as the result.
    std::complex<float> restrict* data = buffer->GetData();
    std::copy_n(model_buffers.front()->GetData(), data_size, data);

    for (size_t i = 1; i < model_buffers.size(); ++i) {
      const std::complex<float> restrict* other_data =
          model_buffers[i]->GetData();
      for (size_t j = 0; j < data_size; ++j) data[j] += other_data[j];
    }
  } else {
    throw std::runtime_error("BdaDdeCal does not support solving yet.");
    // solver_buffer_.AppendAndWeight(*buffer, std::move(model_buffers));

    // if (solver_buffer_.HasCompleteInterval()) {
    // TODO: solver_.Solve(solver_buffer_);
    //  solver_buffer.AdvanceInterval();
    //}
  }

  getNextStep()->process(std::move(buffer));
  return true;
}

void BdaDdeCal::finish() {
  // Solve the remaining data in solver_buffer_.
  // TODO: solver.Solve(solver_buffer_);
  // solver_buffer_.Clear();
}

void BdaDdeCal::show(std::ostream& stream) const {
  stream << "BdaDdeCal " << settings_.name << '\n'
         << "  directions:          " << settings_.directions << '\n';

  for (size_t dir = 0; dir < steps_.size(); ++dir) {
    stream << "Model steps for direction " << settings_.directions[dir] << '\n';
    for (std::shared_ptr<Step> step = steps_[dir]; step;
         step = step->getNextStep()) {
      step->show(stream);
    }
    stream << '\n';
  }
}

}  // namespace steps
}  // namespace dp3
