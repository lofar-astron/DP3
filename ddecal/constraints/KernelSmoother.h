// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_KERNEL_SMOOTHER_H_
#define DP3_DDECAL_KERNEL_SMOOTHER_H_

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>
#include <limits>

namespace dp3 {
namespace ddecal {

/**
 * \brief Smooths a series of possibly irregularly gridded values by a
 * given kernel.
 *
 * The class is optimized to smooth many series which are all placed on the same
 * grid. This is the case when smoothing the solutions on a (a possibly
 * irregular) channel grid.
 *
 * This class uses internally stored scratch space, and is therefore not thread
 * safe. To smooth with multiple threads, instantiate a KernelSmoother for each
 * thread.
 */
template <typename DataType, typename NumType>
class KernelSmoother {
 public:
  enum KernelType {
    RectangularKernel,
    TriangularKernel,
    /** Gaussian, trimmed off at 3 sigma */
    GaussianKernel,
    /** The Epanechnikov kernel is a quadratic kernel, given by 3/4 (1 - x^2) */
    EpanechnikovKernel
  };

  /**
   * Construct and initialize kernel smoother.
   * With se the @c spectralExponent, s the @c kernelBandwidth,
   * nu the frequency at which to evaluate the kernel
   * the kernel size is calculated as:
   * s(nu) = (nu_0 / nu) ^ se
   * @param frequencies Vector defining the channel frequencies: frequencies[i]
   * specifies the frequency of channel i in Hz. The size of this vector (number
   * of channels) defines the size of the grid @c n.
   * @param kernelType Type of kernel to use for smoothing
   * @param kernelBandwidth Size of the kernel (smoothing strength) in frequency
   * units (Hz). May be 0.0 to disable frequency correction of the kernel size.
   * @param bandwidthRefFrequency If non-zero, sets the frequency value at which
   * the kernel size equals the specified size with @c kernelBandwidth.
   * @param spectralExponent A value that specifies the exponent in the relative
   * spectral factor of the kernel size.
   */
  KernelSmoother(const std::vector<NumType>& frequencies, KernelType kernelType,
                 NumType kernelBandwidth, NumType bandwidthRefFrequency,
                 NumType spectralExponent, bool kernel_truncation)
      : _frequencies(frequencies),
        _scratch(frequencies.size()),
        _kernelType(kernelType),
        _bandwidth(kernelBandwidth),
        _bandwidthRefFrequency(bandwidthRefFrequency),
        _spectralExponent(spectralExponent),
        _truncate(kernel_truncation) {}

  /**
   * Evaluate the kernel for a given position.
   * @param distance Distance (positive or negative) from centre of the kernel
   * to evaluate the kernel for, in units of frequency (Hz).
   * @returns Unnormalized kernel value (i.e., integral of kernel is not
   * necessarily unity).
   */
  NumType Kernel(NumType distance) const {
    NumType x = distance / _bandwidth;
    switch (_kernelType) {
      case RectangularKernel:
      default:
        if (x < NumType(-1.0) || x > NumType(1.0))
          return NumType(0.0);
        else
          return NumType(0.5);
      case TriangularKernel:
        if (x < NumType(-1.0) || x > NumType(1.0))
          return NumType(0.0);
        else
          return x >= NumType(0.0) ? (NumType(1.0) - x) : (NumType(1.0) + x);
      case GaussianKernel:
        /// e^(-x^2 / sigma^2), sigma = bandwidth / 3.
        if (_truncate && (x < NumType(-1.0) || x > NumType(1.0)))
          return 0.0;
        else
          return std::exp(-x * x * NumType(9.0));
      case EpanechnikovKernel:
        /// 3/4 * (1-x)^2;
        if (_truncate && (x < NumType(-1.0) || x > NumType(1.0))) {
          return 0.0;
        } else {
          x = NumType(1.0) - x;
          return (NumType(3.0) / NumType(4.0)) * x * x;
        }
    }
  }

  /**
   * Replaces the data with a smoothed version of the data.
   * @param data Data array of size @c n (as specified in constructor) that is
   * smoothed on output.
   * @param weight Associated weights, array of size @c n.
   * @param kernelSizeFactor A factor that is applied to the frequency axis,
   * such that a larger factors causes less smoothing.
   */
  void Smooth(DataType* data, const NumType* weight, NumType kernelSizeFactor) {
    size_t n = _frequencies.size();

    size_t bandLeft = 0;
    /// find right-most kernel value corresponding to the first element of data
    const size_t lbound =
        std::lower_bound(_frequencies.begin(), _frequencies.end(),
                         _frequencies[0] + _bandwidth * 0.5) -
        _frequencies.begin();
    size_t bandRight = std::min(lbound + 1u, n);
    for (size_t i = 0; i != n; ++i) {
      size_t start = 0;
      size_t end = n;
      const NumType localBandwidth =
          _bandwidthRefFrequency == 0.0
              ? _bandwidth
              : _bandwidth * _frequencies[i] / _bandwidthRefFrequency;
      if (_truncate) {
        /// If a boundary is further than half the bandwidth away, move boundary
        while (_frequencies[bandLeft] < _frequencies[i] - localBandwidth * 0.5)
          ++bandLeft;
        while (bandRight != n &&
               _frequencies[bandRight] < _frequencies[i] + localBandwidth * 0.5)
          ++bandRight;

        /// A value of 1 is added to make sure we are not skipping a value
        /// because of rounding errors (kernel will correctly evaluate past
        /// boundaries, so including an unnecessary value has no effect)
        start = bandLeft > 0 ? bandLeft - 1 : 0;
        end = bandRight < n ? bandRight + 1 : n;
      }

      DataType sum(0.0);
      NumType weightSum(0.0);
      const NumType frequencyCorrection =
          std::pow(localBandwidth / _bandwidth, _spectralExponent);
      const NumType kernelCorrection = frequencyCorrection * kernelSizeFactor;
      for (size_t j = start; j != end; ++j) {
        NumType distance = _frequencies[i] - _frequencies[j];
        double w = Kernel(distance * kernelCorrection) * weight[j];
        sum += data[j] * w;
        weightSum += w;
      }
      if (weightSum == 0.0)
        _scratch[i] = quiet_NaN(_scratch[i]);
      else
        _scratch[i] = sum / weightSum;
    }
    std::copy(_scratch.begin(), _scratch.end(), data);
  }

 private:
  double quiet_NaN(double) { return std::numeric_limits<double>::quiet_NaN(); }

  std::complex<double> quiet_NaN(std::complex<double>) {
    return std::complex<double>(std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN());
  }

  std::vector<NumType> _frequencies;
  std::vector<DataType> _scratch;
  enum KernelType _kernelType;
  NumType _bandwidth;
  NumType _bandwidthRefFrequency;
  NumType _spectralExponent;
  bool _truncate = true;
};

}  // namespace ddecal
}  // namespace dp3

#endif
