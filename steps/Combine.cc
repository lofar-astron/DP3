// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Combine.h"

#include <cassert>
#include <iostream>

#include "base/FlagCounter.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3::steps {

Combine::Combine(const common::ParameterSet& parset, const std::string& prefix,
                 MsType ms_type)
    : name_(prefix),
      buffer_name_(parset.getString(prefix + "buffername")),
      ms_type_(ms_type) {
  SetOperation(parset.getString(prefix + "operation", "subtract"));
}

void Combine::updateInfo(const DPInfo& info_in) { Step::updateInfo(info_in); }

void Combine::SetOperation(const std::string& operation) {
  if (operation == "add") {
    operation_ = Operation::kAdd;
  } else if (operation == "subtract") {
    operation_ = Operation::kSubtract;
  } else {
    throw std::invalid_argument("Operation must be 'add' or 'subtract'.");
  }
}

void Combine::show(std::ostream& os) const {
  os << "Combine " << name_ << '\n';
  os << "  buffername: " << buffer_name_ << '\n';
  os << "  operation:  ";
  switch (operation_) {
    case Operation::kAdd:
      os << "add\n";
      break;
    case Operation::kSubtract:
      os << "subtract\n";
      break;
  }
}

void Combine::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " Combine " << name_ << '\n';
}

bool Combine::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  const xt::xtensor<std::complex<float>, 3>& other_data =
      buffer->GetData(buffer_name_);

  switch (operation_) {
    case Operation::kAdd:
      buffer->GetData() += other_data;
      break;
    case Operation::kSubtract:
      buffer->GetData() -= other_data;
      break;
  }

  buffer->RemoveData(buffer_name_);

  timer_.stop();

  getNextStep()->process(std::move(buffer));
  return false;
}

bool Combine::process(std::unique_ptr<base::BdaBuffer> buffer) {
  timer_.start();

  const std::complex<float>* other_data = buffer->GetData(buffer_name_);
  if (other_data == nullptr)
    throw std::runtime_error("Combine step: Could not find buffer with name " +
                             buffer_name_);

  std::complex<float>* destination = buffer->GetData();
  assert(destination != nullptr);

  for (const base::BdaBuffer::Row& row : buffer->GetRows()) {
    const size_t offset = row.offset;
    const size_t data_size = row.GetDataSize();
    switch (operation_) {
      case Operation::kAdd:
        std::transform(&destination[offset], &destination[offset + data_size],
                       &other_data[offset], &destination[offset], std::plus());
        break;
      case Operation::kSubtract:
        std::transform(&destination[offset], &destination[offset + data_size],
                       &other_data[offset], &destination[offset], std::minus());
        break;
    }
  }

  buffer->RemoveData(buffer_name_);
  timer_.stop();

  getNextStep()->process(std::move(buffer));
  return false;
}

void Combine::finish() { getNextStep()->finish(); }

}  // namespace dp3::steps
