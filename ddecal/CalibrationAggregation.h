// See the LICENSE file at the top-level directory of this distribution.

#ifndef DDECAL_CALIBRATIONAGGREGATION_H_
#define DDECAL_CALIBRATIONAGGREGATION_H_

#include <algorithm>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <schaapcommon/threading/staticfor.h>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/views/xview.hpp>

#include <aocommon/logger.h>

#include "ska-sdp-func/utility/sdp_errors.h"
#include "ska-sdp-func/utility/sdp_mem.h"
#include "ska-sdp-func/visibility/sdp_cal_aggregation.h"

namespace dp3::ddecal {

/**
 * Memory-reducing calibration aggregation for the regular (non-BDA) DDECal
 * path.
 * see ska-sdp-func for more explanation
 * The CalibrationAggregation class holds the partials that are incremented at
 * each time sample.
 */
class CalibrationAggregation {
 public:
  using Complex = std::complex<float>;

  /**
   * @param n_batches Number of independent systems to aggregate. In DDECal
   * this is n_baselines * n_channels: one batch per (baseline, frequency)
   * pair, accumulated over the timesteps of the solution interval.
   * @param n_directions Number of directions, ndir.
   * @param n_pols Number of polarizations, npol. DDECal uses 4, in the order
   *               [XX, XY, YX, YY].
   */
  CalibrationAggregation(size_t n_batches, size_t n_directions, size_t n_pols)
      : n_batches_(n_batches),
        n_pols_(n_pols),
        lin(n_directions * n_pols),
        mm(xt::xtensor<Complex, 3>::from_shape({n_batches, lin, lin})),
        vm(xt::xtensor<Complex, 3>::from_shape({n_batches, lin, n_pols})),
        q(xt::xtensor<Complex, 4>::from_shape(
            {n_batches, lin, n_directions, n_pols})),
        q_inv(xt::xtensor<Complex, 3>::from_shape({n_batches, lin, lin})),
        vis(xt::xtensor<Complex, 3>::from_shape({n_batches, lin, n_pols})) {}

  CalibrationAggregation(const CalibrationAggregation&) = delete;
  CalibrationAggregation& operator=(const CalibrationAggregation&) = delete;
  CalibrationAggregation(CalibrationAggregation&&) = delete;
  CalibrationAggregation& operator=(CalibrationAggregation&&) = delete;

  /**
   * Zero the MM and VM accumulators. Must be called before the first
   * `Accumulate()` of a new solution interval.
   */
  void Reset() {
    mm.fill(Complex(0.0f, 0.0f));
    vm.fill(Complex(0.0f, 0.0f));
  }

  void ReplaceWithIdentity(size_t batch) {
    if (batch >= n_batches_) {
      throw std::out_of_range("CalibrationAggregation: batch out of range");
    }
    for (size_t row = 0; row < lin; ++row) {
      for (size_t column = 0; column < lin; ++column) {
        mm(batch, row, column) = row == column ? Complex(1.0f) : Complex(0.0f);
      }
      for (size_t pol = 0; pol < n_pols_; ++pol) {
        vm(batch, row, pol) = Complex(0.0f);
      }
    }
  }

  /**
   * Fold one batch of samples into MM and VM: MM += M^H M, VM += M^H V.
   *
   * Accumulates in place and is safe to call repeatedly with different numbers
   * of samples, which is what lets DDECal feed timesteps one at a time. Samples
   * must already be weighted and their flagged entries zeroed, exactly as
   * ddecal::AssignAndWeight produces them, and all samples of a batch must
   * belong to the same solution index.
   *
   * Taking strides rather than tensors is what makes the aggregation
   * worthwhile: DDECal feeds one timestep at a time, so the sample axis is a
   * single timestep, and the (baseline, frequency) batch axis is a strided view
   * over the existing DPBuffer layout. Peak memory then stays independent of
   * the solution interval length instead of needing a contiguous staging copy.
   *
   * The batches are independent, so they are accumulated in parallel over the
   * (baseline, frequency) axis.
   *
   * @param model Model matrices [n_batches, n_samples, n_directions, n_pols],
   *        given as a base pointer plus shape plus strides.
   * @param model_shape Element count of each model dimension.
   * @param model_strides Element stride of each model dimension.
   * @param vis Observed visibilities [n_batches, n_samples, n_pols].
   * @param vis_shape Element count of each visibility dimension.
   * @param vis_strides Element stride of each visibility dimension.
   *
   * The batch, direction and polarization extents are not re-checked against
   * the values this object was configured with: the caller derives them from
   * the same numbers it passed to the constructor, so a disagreement could only
   * be a caller bug, and sdp reports the resulting bad system as a singular
   * square root rather than silently. What is checked is the rank and the
   * innermost stride, because those are read out of the vectors and out of the
   * caller's buffer respectively.
   */
  void Accumulate(const Complex* model, const std::vector<size_t>& model_shape,
                  const std::vector<size_t>& model_strides, const Complex* vis,
                  const std::vector<size_t>& vis_shape,
                  const std::vector<size_t>& vis_strides) {
    // The kernels index the shape array directly and the stride of the last
    // axis is read below, so both vectors must have exactly the rank expected:
    // a short stride vector would be read out of bounds by the check that
    // follows, and a short shape by sdp itself.
    if (model_shape.size() != 4 || model_strides.size() != 4 ||
        vis_shape.size() != 3 || vis_strides.size() != 3) {
      throw std::runtime_error(
          "CalibrationAggregation::Accumulate: expected a rank-4 model "
          "and a rank-3 visibility, each with as many strides.");
    }
    // The last dimension must be contiguous for the BLAS kernels to see the
    // polarization axis correctly; a strided inner axis is not supported.
    if (model_strides[3] != 1 || vis_strides[2] != 1) {
      throw std::runtime_error(
          "CalibrationAggregation::Accumulate: the polarization axis "
          "must be contiguous (stride 1).");
    }

    // sdp_Mem expects the stride of each dimension in BYTES (not elements), so
    // convert the element strides the caller gave us. Passing the element
    // strides straight through understates every stride by a factor of
    // sizeof(Complex), which reads outside the intended rows and silently
    // yields a non-finite reduced system.
    const std::vector<int64_t> model_shape_i = Extent(model_shape);
    const std::vector<int64_t> vis_shape_i = Extent(vis_shape);
    const std::vector<int64_t> model_stride_i = Bytes(model_strides);
    const std::vector<int64_t> vis_stride_i = Bytes(vis_strides);

    const bool outer_product_path = model_strides[2] == n_pols_ &&
                                    model_strides[3] == 1 &&
                                    vis_strides[2] == 1;

    // Every batch is independent, so the (baseline, frequency) batch axis is
    // split into chunks that are accumulated in parallel. This is where the
    // bulk of the calibration aggregation time is spent.
    RunBatches(
        model_shape[0], "CalibrationAggregation::Accumulate",
        [&](size_t begin, size_t end, sdp_Error* error) {
          if (outer_product_path) {
            AccumulateOuterProducts(model, model_strides, vis, vis_strides,
                                    model_shape[1], begin, end);
            return;
          }
          const size_t count = end - begin;
          // The caller's batch stride is kept; only the base pointer
          // moves to the first batch of this chunk and the batch extent
          // shrinks.
          std::vector<int64_t> model_view_shape = model_shape_i;
          model_view_shape[0] = static_cast<int64_t>(count);
          std::vector<int64_t> vis_view_shape = vis_shape_i;
          vis_view_shape[0] = static_cast<int64_t>(count);

          const MemHandle model_mem(model + begin * model_strides[0],
                                    model_view_shape, &model_stride_i);
          const MemHandle vis_mem(vis + begin * vis_strides[0], vis_view_shape,
                                  &vis_stride_i);
          // Accumulate into the batch range [begin, end) of MM and VM.
          // The view keeps the underlying batch stride, and MemHandle
          // takes the view's data offset into account, so the kernel
          // sees one batch block starting at the right place.
          auto mm_view = xt::view(mm, xt::range(begin, end));
          auto vm_view = xt::view(vm, xt::range(begin, end));
          const MemHandle mm_mem(mm_view.data() + mm_view.data_offset(),
                                 Extent(mm_view.shape()), nullptr);
          const MemHandle vm_mem(vm_view.data() + vm_view.data_offset(),
                                 Extent(vm_view.shape()), nullptr);
          // Allocate a thread_local variable which is alive per thread
          // this avoids reallocation but it means that the workspace
          // stays allocated if part of a thread_pool
          thread_local xt::xtensor<Complex, 1> workspace;
          if (workspace.size() != lin + n_pols_) {
            workspace = xt::xtensor<Complex, 1>::from_shape({lin + n_pols_});
          }
          const MemHandle workspace_mem(workspace);

          sdp_cal_aggregate_accumulate(model_mem, vis_mem, mm_mem, vm_mem,
                                       workspace_mem, error);
        });
  }

  /**
   * Compute Q = (M^H M)^(1/2) and Q_inv = (M^H M)^(-1/2) from the accumulated
   * MM. Destroys MM, so this must follow the last `Accumulate()` and the
   * interval cannot be extended afterwards.
   */
  void SquareRoot() {
    // Each batch has its own MM, so its eigendecomposition is independent and
    // the batches can be split over threads exactly like the accumulation.
    RunBatches(n_batches_, "CalibrationAggregation::SquareRoot",
               [&](size_t begin, size_t end, sdp_Error* error) {
                 auto mm_view = xt::view(mm, xt::range(begin, end));
                 auto q_view = xt::view(q, xt::range(begin, end));
                 auto q_inv_view = xt::view(q_inv, xt::range(begin, end));
                 const MemHandle mm_mem(mm_view.data() + mm_view.data_offset(),
                                        Extent(mm_view.shape()), nullptr);
                 const MemHandle q_mem(q_view.data() + q_view.data_offset(),
                                       Extent(q_view.shape()), nullptr);
                 const MemHandle q_inv_mem(
                     q_inv_view.data() + q_inv_view.data_offset(),
                     Extent(q_inv_view.shape()), nullptr);

                 // Allocate a thread_local variable which is alive per thread
                 // this avoids reallocation but it means that the workspace
                 // stays allocated if part of a thread_pool
                 thread_local xt::xtensor<Complex, 1> workspace;
                 if (workspace.size() != lin * lin + lin) {
                   workspace =
                       xt::xtensor<Complex, 1>::from_shape({lin * lin + lin});
                 }
                 const MemHandle workspace_mem(workspace);

                 sdp_cal_aggregate_square_root(mm_mem, q_mem, q_inv_mem,
                                               workspace_mem, error);
               });
  }

  /** Compute aggregated_vis = Q_inv * VM. */
  void Apply() {
    // As in SquareRoot(), each batch's Q_inv * VM is independent of the others.
    RunBatches(n_batches_, "CalibrationAggregation::Apply",
               [&](size_t begin, size_t end, sdp_Error* error) {
                 auto vm_view = xt::view(vm, xt::range(begin, end));
                 auto q_inv_view = xt::view(q_inv, xt::range(begin, end));
                 auto vis_view = xt::view(vis, xt::range(begin, end));
                 const MemHandle vm_mem(vm_view.data() + vm_view.data_offset(),
                                        Extent(vm_view.shape()), nullptr);
                 const MemHandle q_inv_mem(
                     q_inv_view.data() + q_inv_view.data_offset(),
                     Extent(q_inv_view.shape()), nullptr);
                 const MemHandle vis_mem(
                     vis_view.data() + vis_view.data_offset(),
                     Extent(vis_view.shape()), nullptr);
                 sdp_cal_aggregate_apply(vm_mem, q_inv_mem, vis_mem, error);
               });
  }

  /**
   * The reduced model Q, [n_batches, ndir * npol, n_directions, n_pols]. Valid
   * after SquareRoot().
   */
  const xt::xtensor<Complex, 4>& ReducedModel() const { return q; }

  /**
   * The reduced visibilities, [n_batches, ndir * npol, n_pols]. Valid after
   * Apply().
   */
  const xt::xtensor<Complex, 3>& ReducedVisibility() const { return vis; }

 private:
  /**
   * RAII holder for an sdp_Mem* that does not own the underlying memory. It
   * converts to sdp_Mem* so a view can be handed straight to the kernels.
   *
   * There is deliberately no default constructor: a handle always describes
   * some buffer, so there is no null state to accidentally pass to a kernel.
   */
  class MemHandle {
   public:
    ~MemHandle() { reset(); }
    MemHandle(const MemHandle&) = delete;
    MemHandle& operator=(const MemHandle&) = delete;

    /**
     * Describe the storage of an xtensor or a view of one. Shape, strides and
     * the data offset are taken from the tensor, so the wrapper describes the
     * buffer that is actually there rather than an assumed dense one. That
     * matters because the tensor's allocator may align the base pointer
     * (xtensor uses an xsimd aligned allocator when XTENSOR_USE_XSIMD is set),
     * because the layout is a build-time choice, and because a view's first
     * element sits at data() + data_offset() rather than at data().
     */
    template <class TensorType>
    explicit MemHandle(TensorType& tensor) {
      const std::vector<int64_t> stride = Bytes(tensor.strides());
      // The kernels address the last axis with a unit stride, so a tensor
      // that is not contiguous there would give wrong answers rather than an
      // error.
      // An axis of extent one is exempt: xtensor gives those a stride of zero,
      // which is harmless because the index along them is always zero and so
      // the stride never takes part in the address.
      if (tensor.shape().back() > 1 && tensor.strides().back() != 1) {
        throw std::runtime_error(
            "CalibrationAggregation: the tensor handed to the aggregation is "
            "not contiguous in its last dimension.");
      }
      Create(tensor.data() + tensor.data_offset(), Extent(tensor.shape()),
             &stride);
    }

    /**
     * Describe a caller-supplied buffer. When `stride` is null the view is
     * dense, which is what the kernels assume for the fast path. Otherwise
     * `stride` gives the byte stride of each dimension. The buffer may be
     * read-only, in which case it must only be handed to a kernel as a const
     * sdp_Mem.
     */
    MemHandle(const void* data, const std::vector<int64_t>& shape,
              const std::vector<int64_t>* stride) {
      // sdp has no const-aware sdp_Mem: it always stores a plain void*. This
      // constructor describes the caller's model and visibilities, which the
      // kernels only read, so casting away const at this one boundary is safe.
      Create(const_cast<void*>(data), shape, stride);
    }

    operator sdp_Mem*() const { return handle_; }

   private:
    void Create(void* data, const std::vector<int64_t>& shape,
                const std::vector<int64_t>* stride) {
      sdp_Error error = SDP_SUCCESS;
      handle_ = sdp_mem_create_wrapper(
          data, SDP_MEM_COMPLEX_FLOAT, SDP_MEM_CPU,
          static_cast<int32_t>(shape.size()), shape.data(),
          stride ? stride->data() : nullptr, &error);
      Check(&error, "sdp_mem_create_wrapper");
    }

    void reset() {
      if (handle_) {
        sdp_mem_free(handle_);
        handle_ = nullptr;
      }
    }

    sdp_Mem* handle_ = nullptr;
  };

  static const char* ErrorName(sdp_Error error) {
    switch (error) {
      case SDP_SUCCESS:
        return "success";
      case SDP_ERR_RUNTIME:
        return "runtime error";
      case SDP_ERR_INVALID_ARGUMENT:
        return "invalid argument";
      case SDP_ERR_DATA_TYPE:
        return "unsupported data type";
      case SDP_ERR_MEM_ALLOC_FAILURE:
        return "memory allocation failure";
      case SDP_ERR_MEM_COPY_FAILURE:
        return "memory copy failure";
      case SDP_ERR_MEM_LOCATION:
        return "unsupported memory location";
    }
    return "unknown error";
  }

  static void Check(sdp_Error* error, const char* function) {
    if (*error != SDP_SUCCESS) {
      throw std::runtime_error(std::string("sdp_cal_aggregate: ") + function +
                               " failed: " + ErrorName(*error));
    }
  }

  /**
   * Outer-product form of Accumulate() for the contiguous DDECal layout:
   *   MM[i,j] += sum_a conj(M[a,i]) * M[a,j]
   *   VM[i,p] += sum_a conj(M[a,i]) * V[a,p]
   * which is exactly what the generic path's two gemmHN(Hermitian, NoTrans)
   * calls compute. `model` and `vis` point at batch 0; each batch's
   * direction/polarization block and polarization block are contiguous, which
   * is what the caller of Accumulate() has already been checked for.
   */
  void AccumulateOuterProducts(const Complex* model,
                               const std::vector<size_t>& model_strides,
                               const Complex* vis,
                               const std::vector<size_t>& vis_strides,
                               size_t n_samples, size_t begin, size_t end) {
    Complex* const mm_data = mm.data();
    Complex* const vm_data = vm.data();
    for (size_t batch = begin; batch != end; ++batch) {
      const float* const M_batch =
          reinterpret_cast<const float*>(model + batch * model_strides[0]);
      const float* const V_batch =
          reinterpret_cast<const float*>(vis + batch * vis_strides[0]);
      float* const mm_batch =
          reinterpret_cast<float*>(mm_data + batch * lin * lin);
      float* const vm_batch =
          reinterpret_cast<float*>(vm_data + batch * lin * n_pols_);
      for (size_t sample = 0; sample != n_samples; ++sample) {
        // The per-sample block is lin (model) / n_pols (vis) complex values,
        // and the sample stride counts complex elements, so double it for the
        // float view.
        const float* const M = M_batch + 2 * sample * model_strides[1];
        const float* const V = V_batch + 2 * sample * vis_strides[1];
        for (size_t i = 0; i != lin; ++i) {
          const float real_i = M[2 * i];
          const float imag_i = M[2 * i + 1];
          float* const mm_row = mm_batch + 2 * i * lin;
          for (size_t j = 0; j != lin; ++j) {
            const float real_j = M[2 * j];
            const float imag_j = M[2 * j + 1];
            mm_row[2 * j] += real_i * real_j + imag_i * imag_j;
            mm_row[2 * j + 1] += real_i * imag_j - imag_i * real_j;
          }
          float* const vm_row = vm_batch + 2 * i * n_pols_;
          for (size_t p = 0; p != n_pols_; ++p) {
            const float real_p = V[2 * p];
            const float imag_p = V[2 * p + 1];
            vm_row[2 * p] += real_i * real_p + imag_i * imag_p;
            vm_row[2 * p + 1] += real_i * imag_p - imag_i * real_p;
          }
        }
      }
    }
  }

  template <class Body>
  static void RunBatches(size_t n_batches, const char* what, Body&& body) {
    if (n_batches == 0) return;
    const size_t max_threads =
        schaapcommon::ThreadPool::GetInstance().NThreads();
    std::vector<sdp_Error> errors(max_threads, SDP_SUCCESS);
    schaapcommon::RunConstrainedStaticFor<size_t>(
        0, n_batches, max_threads,
        [&](size_t begin, size_t end, size_t thread) {
          body(begin, end, &errors[thread]);
        });
    if (std::any_of(errors.begin(), errors.end(), [](const sdp_Error error) {
          return error != SDP_SUCCESS;
        })) {
      for (const sdp_Error error : errors) {
        if (error != SDP_SUCCESS) {
          aocommon::Logger::Error << what << ": " << ErrorName(error) << '\n';
        }
      }
      throw std::runtime_error(std::string(what) +
                               ": one or more chunks failed");
    }
  }

  /** Convert a shape to the int64_t extents sdp_Mem expects. */
  template <class Shape>
  static std::vector<int64_t> Extent(const Shape& shape) {
    std::vector<int64_t> result;
    result.reserve(shape.size());
    for (const size_t value : shape) {
      result.push_back(static_cast<int64_t>(value));
    }
    return result;
  }

  /** Convert element strides to the byte strides sdp_Mem expects. */
  template <class Strides>
  static std::vector<int64_t> Bytes(const Strides& strides) {
    std::vector<int64_t> result;
    result.reserve(strides.size());
    for (const size_t stride : strides) {
      result.push_back(static_cast<int64_t>(stride) *
                       static_cast<int64_t>(sizeof(Complex)));
    }
    return result;
  }

  const size_t n_batches_;
  const size_t n_pols_;
  const size_t lin;

  // MM and Q_inv are [batch, lin, lin]; VM and the aggregated visibilities are
  // [batch, lin, npol]; Q is [batch, lin, ndir, npol]. Declared in constructor
  // initialisation order, so `lin` is set before any of the shapes that use it.
  xt::xtensor<Complex, 3, xt::layout_type::row_major> mm;
  xt::xtensor<Complex, 3, xt::layout_type::row_major> vm;
  xt::xtensor<Complex, 4> q;
  xt::xtensor<Complex, 3> q_inv;
  xt::xtensor<Complex, 3> vis;
};

}  // namespace dp3::ddecal

#endif  // DDECAL_CALIBRATIONAGGREGATION_H_