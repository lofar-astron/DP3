// DPBuffer.cc: Buffer holding the data of a timeslot/band
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <dp3/base/DPBuffer.h>

#include <algorithm>
#include <cassert>

#include <casacore/casa/version.h>
#include <casacore/casa/BasicSL/Complexfwd.h>

// Casacore < 3.4 does not support move semantics for casacore::Array
// and uses reference semantics in the copy constructor.
#if CASACORE_MAJOR_VERSION > 3 || \
    (CASACORE_MAJOR_VERSION == 3 && CASACORE_MINOR_VERSION >= 4)
#define USE_CASACORE_MOVE_SEMANTICS
#endif

namespace dp3 {
namespace base {

DPBuffer::DPBuffer(double time, double exposure)
    : time_(time),
      exposure_(exposure),
      row_numbers_(),
      data_(),
      extra_data_(),
      flags_(),
      weights_(),
      uvw_(),
      solution_() {}

DPBuffer::DPBuffer(DPBuffer&& that)
    : time_(that.time_),
      exposure_(that.exposure_),
      row_numbers_(std::move(that.row_numbers_)),
      data_(std::move(that.data_)),
      extra_data_(std::move(that.extra_data_)),
      flags_(std::move(that.flags_)),
      weights_(std::move(that.weights_)),
      uvw_(std::move(that.uvw_)),
      solution_(that.solution_) {
#ifndef USE_CASACORE_MOVE_SEMANTICS
  // The copy constructor for casacore::Array creates references. Since
  // moving a buffer does not have reference semantics, clear 'that'.
  that.row_numbers_.assign(decltype(that.row_numbers_)());
#endif
}

DPBuffer::DPBuffer(const DPBuffer& that, const common::Fields& fields)
    : DPBuffer(that.time_, that.exposure_) {
  Copy(that, fields);
}

DPBuffer& DPBuffer::operator=(const DPBuffer& that) {
  if (this != &that) {
    time_ = that.time_;
    exposure_ = that.exposure_;
    solution_ = that.solution_;
    row_numbers_.reference(that.row_numbers_);
    data_ = that.data_;
    extra_data_ = that.extra_data_;
    flags_ = that.flags_;
    weights_ = that.weights_;
    uvw_ = that.uvw_;
  }
  return *this;
}

DPBuffer& DPBuffer::operator=(DPBuffer&& that) {
  if (this != &that) {
    time_ = that.time_;
    exposure_ = that.exposure_;
    data_ = std::move(that.data_);
    extra_data_ = std::move(that.extra_data_);
    flags_ = std::move(that.flags_);
    weights_ = std::move(that.weights_);
    uvw_ = std::move(that.uvw_);
    solution_ = std::move(that.solution_);

    // Casacore < 3.4.0 does not support move semantics for casacore::Array.
    // The copy assignment operator for casacore::Array then creates copies.
#ifdef USE_CASACORE_MOVE_SEMANTICS
    row_numbers_ = std::move(that.row_numbers_);
#else
    // Create references. Since move assignment does not use reference
    // semantics, clear 'that'.
    row_numbers_.reference(that.row_numbers_);
    that.row_numbers_.assign(decltype(that.row_numbers_)());
#endif
  }
  return *this;
}

void DPBuffer::copy(const DPBuffer& that) {
  operator=(that);  // 'Copy' the data in 'that', by referencing it.

  // Ensure that this buffer is independent of the data in 'that'.
  // Do it even if 'this == &that', since this/that may contain references,
  // and the result of 'copy' should always be an independent object.
  row_numbers_.unique();
}

void DPBuffer::Copy(const DPBuffer& that, const common::Fields& fields) {
  if (this != &that) {
    time_ = that.time_;
    exposure_ = that.exposure_;
    row_numbers_.reference(that.row_numbers_);
    if (fields.Data()) data_ = that.data_;
    if (fields.Flags()) flags_ = that.flags_;
    if (fields.Weights()) weights_ = that.weights_;
    if (fields.Uvw()) uvw_ = that.uvw_;
    // TODO(AST-1241): Copy extra data fields, too.
    solution_ = that.solution_;
  }
}

void DPBuffer::AddData(const std::string& name) {
  assert(!name.empty());
  assert(extra_data_.find(name) == extra_data_.end());
  extra_data_[name].resize(data_.shape());
}

void DPBuffer::RemoveData(const std::string& name) {
  if (name.empty()) {
    extra_data_.clear();
  } else {
    extra_data_.erase(name);
  }
}

void DPBuffer::CopyData(const DPBuffer& source, const std::string& source_name,
                        const std::string& target_name) {
  assert(source.HasData(source_name));
  assert(!target_name.empty());

  extra_data_[target_name] = source.GetData(source_name);
}

void DPBuffer::MoveData(DPBuffer& source, const std::string& source_name,
                        const std::string& target_name) {
  assert(source.HasData(source_name));
  assert(!target_name.empty());

  if (source_name.empty()) {
    extra_data_[target_name] = std::move(source.data_);
  } else {
    auto source_data_iterator = source.extra_data_.find(source_name);

    extra_data_.insert_or_assign(target_name,
                                 std::move(source_data_iterator->second));

    source.extra_data_.erase(source_data_iterator);
  }
}

}  // namespace base
}  // namespace dp3
