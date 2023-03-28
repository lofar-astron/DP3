// DPBuffer.h: Buffer holding the data of a timeslot/band
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Buffer holding the data of a timeslot/band
/// @author Ger van Diepen

#ifndef DP3_BASE_DPBUFFER_H_
#define DP3_BASE_DPBUFFER_H_

#include <map>
#include <string>

#include <xtensor/xtensor.hpp>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Cube.h>

#include <aocommon/xt/span.h>
#include <aocommon/xt/utensor.h>

#include <dp3/common/Fields.h>
#include <dp3/common/Types.h>

namespace dp3 {
namespace base {

/// @brief Buffer holding the data of a timeslot/band

/// This class holds the data for one time slot in Array variables.
/// It makes heavy use of reference semantics to avoid data copying
/// when data are pushed from one step to another.
/// This means that a data array can be shared between Step objects.
/// So, if a Step object changes data in a buffer, it has to be sure
/// it can do it. If needed, Array::unique should be called to ensure
/// the array is not shared.
///
/// The following data can be kept in a DPBuffer object.
/// <table>
///  <tr>
///   <td>TIME</td>
///   <td>The time slot center of the current data (in MJD seconds).</td>
///  </tr>
///  <tr>
///   <td>ROWNRS</td>
///   <td>The row numbers of the current time slot. It can be empty
///       when e.g. a time slot is inserted or if data are averaged.</td>
///  </tr>
///  <tr>
///   <td>DATA</td>
///   <td>The visibility data as [ncorr,nchan,nbaseline].</td>
///  </tr>
///  <tr>
///   <td>FLAG</td>
///   <td>The data flags as [ncorr,nchan,nbaseline] (True is bad).
///       Note that the ncorr axis is redundant because NDPPP will always
///       have the same flag for all correlations. The reason all
///       correlations are there is because the MS expects them.</td>
///  </tr>
///  <tr>
///   <td>WEIGHT</td>
///   <td>The data weights as [ncorr,nchan,nbaseline].
///       Similarly to FLAG, the ncorr axis is redundant because the
///       same weight is used for all correlations.</td>
///  </tr>
///  <tr>
///   <td>UVW</td>
///   <td>The UVW coordinates in meters as [3,nbaseline].</td>
///  </tr>
///  <tr>
///   <td>FULLRESFLAG</td>
///   <td>The flags before any averaging was done. In the MS they are kept
///       in column LOFAR_FULL_RES_FLAG. They are used to deal in BBS
///       in a smart way with bandwidth and time smearing.
///       The shape of the array is [nchan, ntimeavg, nbaseline], where
///       ntimeavg is the number of time slots averaged to a single one.
///       The number of channels averaged to a single one can be determined
///       by dividing nchan by the number of channels in the data (or flags).
///   </td>
///  </tr>
/// </table>
/// Each data member (DATA, FLAG, UVW, WEIGHTS, FULLRESFLAGS) is filled in if
/// any Step needs it (the information about the required fields per each step
/// can be read with the getRequiredFields() function). The first Step
/// (MSReader) will read the requested fields from the MS into the DPBuffer. In
/// that way as little memory as needed is used. Note that e.g. the AOFlagger
/// can use a lot of memory if a large time window is used.
///
/// Until early 2015, NDPPP used the strategy of shallow data copies.
/// I.e., a step increased the data reference counter and did not make
/// an actual copy. Only when data were changed, a new data array was made.
/// Thus, MSReader allocated a new array when it read the data.
/// However, it appeared this strategy lead to memory fragmentation and
/// to sudden jumps in memory usage on Linux systems.
/// <br>Therefore the strategy was changed to having each step preallocate
/// its buffers and making deep copies when moving data from one step to
/// the next one. It appeared that it not only improved memory usage,
/// but also improved performance, possible due to far less mallocs.
///
/// The buffer/step guidelines are as follows:
/// 1. If a step keeps a buffer of DPBuffers for later processing (e.g.
///    AOFlagger), it must make a copy of the DPBuffers because the input data
///    arrays might have changed before that step processes the data.
/// 2. A shallow copy of a data member can be used if a step processes
///    the data immediately (e.g. Averager).
class DPBuffer {
 public:
  using Complex = std::complex<float>;

  /// Construct object with empty arrays.
  explicit DPBuffer(double time = 0.0, double exposure = 0.0);

  /// The copy constructor uses reference copies.
  DPBuffer(const DPBuffer&) = default;

  /// The move constructor moves all data without using reference semantics.
  DPBuffer(DPBuffer&&);

  /// Copy assignment uses reference copies.
  DPBuffer& operator=(const DPBuffer&);

  /// Move assignment moves all data without using reference semantics.
  DPBuffer& operator=(DPBuffer&&);

  /// Make a deep copy of all arrays in that to this.
  /// After this call, the buffer is always independent.
  void copy(const DPBuffer& that);

  /// Ensure that this buffer has independent copies of data items / that
  /// the data items do not use reference semantics.
  /// When DPBuffer will store its data in XTensor objects instead of Casacore
  /// objects, a DPBuffer is always independent and this function can go away.
  /// @param fields The fields that should be independent copies.
  void MakeIndependent(const common::Fields& fields);

  /// Reference only the arrays that are filled in that.
  void referenceFilled(const DPBuffer& that);

  /// Set or get the visibility data per corr,chan,baseline.
  /// This legacy function only works if the default data buffer is the
  /// main data buffer.
  void setData(const casacore::Cube<Complex>& data);

  /// Adds an extra visibility buffer and makes it the default data buffer.
  /// The new visibility buffer gets the same shape as the main data buffer.
  /// @param key Key for the new buffer. The key may not be empty. It may
  ///        also not equal the key of an existing extra data buffer.
  void AddData(const std::string& key);

  /// Removes extra visibility buffers.
  /// If the default data buffer is removed, the default data buffer becomes
  /// the main data buffer.
  /// @param key Key of the buffer to remove. If empty, removes all extra data
  /// buffers.
  void RemoveData(const std::string& key = "");

  /// Check if the DPBuffer has a data buffer for the given key.
  /// @param key Data buffer key. When empty, check the main data buffer.
  /// @return If the requested data buffer exists. It can be empty, though.
  bool HasData(const std::string& key = "") const {
    return key.empty() || (extra_data_.find(key) != extra_data_.end());
  }

  /// Set the default data buffer, which GetData() should return.
  /// @param key Buffer key. If empty, the main data buffer becomes the default
  ///        data buffer. If non-empty, a data buffer with that key must exist.
  void SetDefaultData(const std::string& key = "");

  void ResizeData(size_t n_baselines, size_t n_channels, size_t n_correlations);
  const casacore::Cube<Complex>& GetCasacoreData() const { return casa_data_; }
  casacore::Cube<Complex>& GetCasacoreData() { return casa_data_; }

  /// @return A const reference to the default data buffer.
  const aocommon::xt::Span<std::complex<float>, 3>& GetData() const {
    return GetData(default_data_key_);
  }

  /// @return A non-const reference to the default data buffer.
  aocommon::xt::Span<std::complex<float>, 3>& GetData() {
    return GetData(default_data_key_);
  }

  /// @return A const reference to the requested data buffer.
  const aocommon::xt::Span<std::complex<float>, 3>& GetData(
      const std::string& key) const {
    if (key.empty()) {
      return data_;
    } else {
      auto found = extra_data_span_.find(key);
      assert(found != extra_data_span_.end());
      return found->second;
    }
  }

  /// @return A non-const reference to the requested data buffer.
  aocommon::xt::Span<std::complex<float>, 3>& GetData(const std::string& key) {
    if (key.empty()) {
      return data_;
    } else {
      auto found = extra_data_span_.find(key);
      assert(found != extra_data_span_.end());
      return found->second;
    }
  }

  /// Set or get the flags per corr,chan,baseline.
  void setFlags(const casacore::Cube<bool>& flags);
  void ResizeFlags(size_t n_baselines, size_t n_channels,
                   size_t n_correlations);
  const casacore::Cube<bool>& GetCasacoreFlags() const { return casa_flags_; }
  casacore::Cube<bool>& GetCasacoreFlags() { return casa_flags_; }
  const aocommon::xt::Span<bool, 3>& GetFlags() const { return flags_; }
  aocommon::xt::Span<bool, 3>& GetFlags() { return flags_; }

  /// Set or get the weights per corr,chan,baseline.
  void setWeights(const casacore::Cube<float>& weights);
  void ResizeWeights(size_t n_baselines, size_t n_channels,
                     size_t n_correlations);
  const casacore::Cube<float>& GetCasacoreWeights() const {
    return casa_weights_;
  }
  casacore::Cube<float>& GetCasacoreWeights() { return casa_weights_; }
  const aocommon::xt::Span<float, 3>& GetWeights() const { return weights_; }
  aocommon::xt::Span<float, 3>& GetWeights() { return weights_; }

  /// Set or get the flags at the full resolution per chan,timeavg,baseline.
  void setFullResFlags(const casacore::Cube<bool>& flags);
  void ResizeFullResFlags(size_t n_baselines, size_t n_averaged_times,
                          size_t n_full_res_channels);
  const casacore::Cube<bool>& GetCasacoreFullResFlags() const {
    return casa_full_res_flags_;
  }
  casacore::Cube<bool>& GetCasacoreFullResFlags() {
    return casa_full_res_flags_;
  }
  const aocommon::xt::Span<bool, 3>& GetFullResFlags() const {
    return full_res_flags_;
  }
  aocommon::xt::Span<bool, 3>& GetFullResFlags() { return full_res_flags_; }

  /// Get or set the time.
  void setTime(double time) { time_ = time; }
  double getTime() const { return time_; }

  /// Get or set the exposure.
  void setExposure(double exposure) { exposure_ = exposure; }
  double getExposure() const { return exposure_; }

  /// Get or set the row numbers used by the InputStep class.
  /// It can be empty (e.g. when MSReader inserted a dummy time slot).
  void setRowNrs(const casacore::Vector<common::rownr_t>& rownrs) {
    row_numbers_.reference(rownrs);
  }
  const casacore::Vector<common::rownr_t>& getRowNrs() const {
    return row_numbers_;
  }
  casacore::Vector<common::rownr_t>& getRowNrs() { return row_numbers_; }

  /// Get or set the UVW coordinates per baseline.
  void setUVW(const casacore::Matrix<double>& uvw);
  void ResizeUvw(size_t n_baselines);
  const casacore::Matrix<double>& GetCasacoreUvw() const { return casa_uvw_; }
  casacore::Matrix<double>& GetCasacoreUvw() { return casa_uvw_; }
  const aocommon::xt::Span<double, 2>& GetUvw() const { return uvw_; }
  aocommon::xt::Span<double, 2>& GetUvw() { return uvw_; }

  /// Merge the flags into the pre-average flags.
  /// For each flagged point, the corresponding pre-average flags are set.
  void MergeFullResFlags();

  void SetSolution(
      const std::vector<std::vector<std::complex<double>>>& solution) {
    solution_ = solution;
  }
  const std::vector<std::vector<std::complex<double>>>& GetSolution() const {
    return solution_;
  }

 private:
  /// (Re)create the XTensor views on the Casacore data.
  void CreateSpans();

  double time_;
  double exposure_;
  casacore::Vector<common::rownr_t> row_numbers_;

  // The casa_ objects hold the actual data. The aocommon::xt::Span objects
  // below provide XTensor views to the data in the casa_ objects.
  // In the future, the XTensor views will be replaced by xt::xtensor objects
  // that hold the data. The casa_ objects then provide a Casacore view to the
  // data. When the casa_ objects are no longer used, they will be removed.
  casacore::Cube<Complex> casa_data_;         ///< ncorr,nchan,nbasel
  casacore::Cube<bool> casa_flags_;           ///< ncorr,nchan,nbasel
  casacore::Matrix<double> casa_uvw_;         ///< 3,nbasel
  casacore::Cube<float> casa_weights_;        ///< ncorr,nchan,nbasel
  casacore::Cube<bool> casa_full_res_flags_;  ///< fullres_nchan,ntimeavg,nbl

  /// Visibilities (n_baselines x n_channels x n_correlations)
  aocommon::xt::Span<std::complex<float>, 3> data_;
  /// Extra visibilities, e.g., containing predictions for different directions.
  std::map<std::string, aocommon::xt::UTensor<std::complex<float>, 3>>
      extra_data_;
  /// For now, create spans for the extra data, so GetData() can always return
  /// a reference to a span.
  std::map<std::string, aocommon::xt::Span<std::complex<float>, 3>>
      extra_data_span_;
  /// Key (name) of the default data buffer returned by GetData(). An empty
  /// string indicates the main data_ buffer.
  std::string default_data_key_;

  /// Flags (n_baselines x n_channels x n_correlations)
  aocommon::xt::Span<bool, 3> flags_;
  /// Weights (n_baselines x n_channels x n_correlations)
  aocommon::xt::Span<float, 3> weights_;
  /// UVW coordinates (n_baselines x 3)
  aocommon::xt::Span<double, 2> uvw_;
  /// LOFAR full resolution flags
  /// (n_baselines x n_averaged_times x n_full_resolution_channels)
  aocommon::xt::Span<bool, 3> full_res_flags_;

  std::vector<std::vector<std::complex<double>>>
      solution_;  ///< nchan,nant*npol
};

}  // namespace base
}  // namespace dp3

#endif  // DP3_BASE_DPBUFFER_H_
