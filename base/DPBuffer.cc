// DPBuffer.cc: Buffer holding the data of a timeslot/band
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <dp3/base/DPBuffer.h>

#include <cassert>

#include <casacore/casa/version.h>

// Casacore < 3.4 does not support move semantics for casacore::Array
// and uses reference semantics in the copy constructor.
#if CASACORE_MAJOR_VERSION > 3 || \
    (CASACORE_MAJOR_VERSION == 3 && CASACORE_MINOR_VERSION >= 4)
#define USE_CASACORE_MOVE_SEMANTICS
#endif

namespace {

template <typename T>
aocommon::xt::Span<T, 3> CreateSpan(casacore::Cube<T>& cube) {
  const std::array<size_t, 3> shape{static_cast<size_t>(cube.shape()[2]),
                                    static_cast<size_t>(cube.shape()[1]),
                                    static_cast<size_t>(cube.shape()[0])};
  return aocommon::xt::CreateSpan(cube.data(), shape);
}

template <typename T>
aocommon::xt::Span<T, 2> CreateSpan(casacore::Matrix<T>& matrix) {
  const std::array<size_t, 2> shape{static_cast<size_t>(matrix.shape()[1]),
                                    static_cast<size_t>(matrix.shape()[0])};
  return aocommon::xt::CreateSpan(matrix.data(), shape);
}

}  // namespace

namespace dp3 {
namespace base {

DPBuffer::DPBuffer(double time, double exposure)
    : time_(time),
      exposure_(exposure),
      row_numbers_(),
      casa_data_(),
      casa_flags_(),
      casa_uvw_(),
      casa_weights_(),
      casa_full_res_flags_(),
      data_(CreateSpan(casa_data_)),
      extra_data_(),
      extra_data_span_(),
      flags_(CreateSpan(casa_flags_)),
      weights_(CreateSpan(casa_weights_)),
      uvw_(CreateSpan(casa_uvw_)),
      full_res_flags_(CreateSpan(casa_full_res_flags_)),
      solution_() {}

DPBuffer::DPBuffer(DPBuffer&& that)
    : time_(that.time_),
      exposure_(that.exposure_),
      row_numbers_(std::move(that.row_numbers_)),
      casa_data_(std::move(that.casa_data_)),
      casa_flags_(std::move(that.casa_flags_)),
      casa_uvw_(std::move(that.casa_uvw_)),
      casa_weights_(std::move(that.casa_weights_)),
      casa_full_res_flags_(std::move(that.casa_full_res_flags_)),
      data_(std::move(that.data_)),
      extra_data_(std::move(that.extra_data_)),
      extra_data_span_(std::move(that.extra_data_span_)),
      flags_(std::move(that.flags_)),
      weights_(std::move(that.weights_)),
      uvw_(std::move(that.uvw_)),
      full_res_flags_(std::move(that.full_res_flags_)),
      solution_(that.solution_) {
#ifndef USE_CASACORE_MOVE_SEMANTICS
  // The copy constructor for casacore::Array creates references. Since
  // moving a buffer does not have reference semantics, clear 'that'.
  that.row_numbers_.assign(decltype(that.row_numbers_)());
  that.casa_data_.assign(decltype(that.casa_data_)());
  that.casa_flags_.assign(decltype(that.casa_flags_)());
  that.casa_uvw_.assign(decltype(that.casa_uvw_)());
  that.casa_weights_.assign(decltype(that.casa_weights_)());
  that.casa_full_res_flags_.assign(decltype(that.casa_full_res_flags_)());
#endif
  that.CreateSpans();
}

DPBuffer& DPBuffer::operator=(const DPBuffer& that) {
  if (this != &that) {
    time_ = that.time_;
    exposure_ = that.exposure_;
    solution_ = that.solution_;
    row_numbers_.reference(that.row_numbers_);
    casa_data_.reference(that.casa_data_);
    extra_data_ = that.extra_data_;
    casa_flags_.reference(that.casa_flags_);
    casa_weights_.reference(that.casa_weights_);
    casa_uvw_.reference(that.casa_uvw_);
    casa_full_res_flags_.reference(that.casa_full_res_flags_);
    CreateSpans();
  }
  return *this;
}

DPBuffer& DPBuffer::operator=(DPBuffer&& that) {
  if (this != &that) {
    time_ = that.time_;
    exposure_ = that.exposure_;
    extra_data_ = std::move(that.extra_data_);
    solution_ = std::move(that.solution_);

    // Casacore < 3.4.0 does not support move semantics for casacore::Array.
    // The copy assignment operator for casacore::Array then creates copies.
#ifdef USE_CASACORE_MOVE_SEMANTICS
    row_numbers_ = std::move(that.row_numbers_);
    casa_data_ = std::move(that.casa_data_);
    casa_flags_ = std::move(that.casa_flags_);
    casa_weights_ = std::move(that.casa_weights_);
    casa_uvw_ = std::move(that.casa_uvw_);
    casa_full_res_flags_ = std::move(that.casa_full_res_flags_);
#else
    // Create references. Since move assignment does not use reference
    // semantics, clear 'that'.
    row_numbers_.reference(that.row_numbers_);
    casa_data_.reference(that.casa_data_);
    casa_flags_.reference(that.casa_flags_);
    casa_weights_.reference(that.casa_weights_);
    casa_uvw_.reference(that.casa_uvw_);
    casa_full_res_flags_.reference(that.casa_full_res_flags_);
    that.row_numbers_.assign(decltype(that.row_numbers_)());
    that.casa_data_.assign(decltype(that.casa_data_)());
    that.casa_flags_.assign(decltype(that.casa_flags_)());
    that.casa_uvw_.assign(decltype(that.casa_uvw_)());
    that.casa_weights_.assign(decltype(that.casa_weights_)());
    that.casa_full_res_flags_.assign(decltype(that.casa_full_res_flags_)());
#endif
    CreateSpans();
    that.CreateSpans();
  }
  return *this;
}

void DPBuffer::copy(const DPBuffer& that) {
  operator=(that);  // 'Copy' the data in 'that', by referencing it.

  // Ensure that this buffer is independent of the data in 'that'.
  // Do it even if 'this == &that', since this/that may contain references,
  // and the result of 'copy' should always be an independent object.
  row_numbers_.unique();
  casa_data_.unique();
  casa_flags_.unique();
  casa_weights_.unique();
  casa_uvw_.unique();
  casa_full_res_flags_.unique();
  CreateSpans();
}

void DPBuffer::MakeIndependent(const common::Fields& fields) {
  if (fields.Data()) {
    casa_data_.unique();
    data_ = CreateSpan(casa_data_);
  }
  if (fields.Flags()) {
    casa_flags_.unique();
    flags_ = CreateSpan(casa_flags_);
  }
  if (fields.Weights()) {
    casa_weights_.unique();
    weights_ = CreateSpan(casa_weights_);
  }
  if (fields.Uvw()) {
    casa_uvw_.unique();
    uvw_ = CreateSpan(casa_uvw_);
  }
  if (fields.FullResFlags()) {
    casa_full_res_flags_.unique();
    full_res_flags_ = CreateSpan(casa_full_res_flags_);
  }
}

void DPBuffer::CreateSpans() {
  data_ = CreateSpan(casa_data_);
  flags_ = CreateSpan(casa_flags_);
  weights_ = CreateSpan(casa_weights_);
  uvw_ = CreateSpan(casa_uvw_);
  full_res_flags_ = CreateSpan(casa_full_res_flags_);
  extra_data_span_.clear();
  for (auto& extra_element : extra_data_) {
    const std::string& name = extra_element.first;
    extra_data_span_.emplace(
        std::make_pair(name, aocommon::xt::CreateSpan(extra_data_[name])));
  }
}

void DPBuffer::referenceFilled(const DPBuffer& that) {
  if (this != &that) {
    // extra_data does not support the legacy referenceFilled function.
    assert(that.extra_data_.empty());

    time_ = that.time_;
    exposure_ = that.exposure_;
    row_numbers_.reference(that.row_numbers_);
    if (!that.casa_data_.empty()) {
      casa_data_.reference(that.casa_data_);
      data_ = CreateSpan(casa_data_);
    }
    if (!that.casa_flags_.empty()) {
      casa_flags_.reference(that.casa_flags_);
      flags_ = CreateSpan(casa_flags_);
    }
    if (!that.casa_weights_.empty()) {
      casa_weights_.reference(that.casa_weights_);
      weights_ = CreateSpan(casa_weights_);
    }
    if (!that.casa_uvw_.empty()) {
      casa_uvw_.reference(that.casa_uvw_);
      uvw_ = CreateSpan(casa_uvw_);
    }
    if (!that.casa_full_res_flags_.empty()) {
      casa_full_res_flags_.reference(that.casa_full_res_flags_);
      full_res_flags_ = CreateSpan(casa_full_res_flags_);
    }
  }
}

void DPBuffer::MergeFullResFlags() {
  // Flag shape is [nbl, newnchan, ncorr].
  // FullRes shape is [nbl, navgtime, orignchan]
  // where orignchan = navgchan * newnchan.
  const size_t orignchan = full_res_flags_.shape(2);
  const size_t newnchan = flags_.shape(1);
  const size_t navgchan = orignchan / newnchan;
  const size_t navgtime = full_res_flags_.shape(1);
  const size_t nbl = full_res_flags_.shape(0);
  assert(nbl == flags_.shape(0));
  const size_t ncorr = flags_.shape(2);
  bool* fullResPtr = full_res_flags_.data();
  const bool* flagPtr = flags_.data();
  // Loop over all baselines and new channels.
  // Only use the first correlation in the loop.
  for (size_t j = 0; j < nbl; ++j) {
    for (size_t i = 0; i < newnchan; ++i) {
      // If ta data point is flagged, the flags in the corresponding
      // FullRes window have to be set.
      // This is needed in case a data point is further averaged.
      if (*flagPtr) {
        for (size_t i = 0; i < navgtime; ++i) {
          std::fill(fullResPtr, fullResPtr + navgchan, true);
          fullResPtr += orignchan;
        }
        fullResPtr -= navgtime * orignchan;
      }
      flagPtr += ncorr;
      fullResPtr += navgchan;
    }
    // Set pointer to next baseline.
    fullResPtr += (navgtime - 1) * orignchan;
  }
}

void DPBuffer::setData(const casacore::Cube<Complex>& data) {
  casa_data_.reference(data);
  data_ = CreateSpan(casa_data_);
  assert(extra_data_.empty() ||
         extra_data_.begin()->second.shape() == data_.shape());
}

void DPBuffer::AddData(const std::string& name) {
  assert(!name.empty());
  assert(extra_data_.find(name) == extra_data_.end());
  extra_data_[name].resize(data_.shape());
  extra_data_span_.emplace(
      std::make_pair(name, aocommon::xt::CreateSpan(extra_data_[name])));
}

void DPBuffer::RemoveData(const std::string& name) {
  if (name.empty()) {
    extra_data_.clear();
    extra_data_span_.clear();
  } else {
    extra_data_.erase(name);
    extra_data_span_.erase(name);
  }
}

void DPBuffer::ResizeData(size_t n_baselines, size_t n_channels,
                          size_t n_correlations) {
  casa_data_.resize(n_correlations, n_channels, n_baselines);
  data_ = CreateSpan(casa_data_);
  for (auto& extra_element : extra_data_) {
    const std::string& name = extra_element.first;
    extra_element.second.resize({n_baselines, n_channels, n_correlations});
    extra_data_span_.find(name)->second =
        aocommon::xt::CreateSpan(extra_data_[name]);
  }
}

void DPBuffer::setFlags(const casacore::Cube<bool>& flags) {
  casa_flags_.reference(flags);
  flags_ = CreateSpan(casa_flags_);
}

void DPBuffer::ResizeFlags(size_t n_baselines, size_t n_channels,
                           size_t n_correlations) {
  casa_flags_.resize(n_correlations, n_channels, n_baselines);
  flags_ = CreateSpan(casa_flags_);
}

void DPBuffer::setWeights(const casacore::Cube<float>& weights) {
  casa_weights_.reference(weights);
  weights_ = CreateSpan(casa_weights_);
}

void DPBuffer::ResizeWeights(size_t n_baselines, size_t n_channels,
                             size_t n_correlations) {
  casa_weights_.resize(n_correlations, n_channels, n_baselines);
  weights_ = CreateSpan(casa_weights_);
}

void DPBuffer::setUVW(const casacore::Matrix<double>& uvw) {
  casa_uvw_.reference(uvw);
  uvw_ = CreateSpan(casa_uvw_);
}

void DPBuffer::ResizeUvw(size_t n_baselines) {
  casa_uvw_.resize(3, n_baselines);
  uvw_ = CreateSpan(casa_uvw_);
}

void DPBuffer::setFullResFlags(const casacore::Cube<bool>& flags) {
  casa_full_res_flags_.reference(flags);
  full_res_flags_ = CreateSpan(casa_full_res_flags_);
}

void DPBuffer::ResizeFullResFlags(size_t n_baselines, size_t n_averaged_times,
                                  size_t n_full_res_channels) {
  casa_full_res_flags_.resize(n_full_res_channels, n_averaged_times,
                              n_baselines);
  full_res_flags_ = CreateSpan(casa_full_res_flags_);
}

}  // namespace base
}  // namespace dp3
