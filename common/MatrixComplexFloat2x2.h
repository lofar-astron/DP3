// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef AOCOMMON_MATRIX_COMPLEX_FLOAT_2X2_H
#define AOCOMMON_MATRIX_COMPLEX_FLOAT_2X2_H

#include "common/scalar/MatrixComplexFloat2x2.h"
#include "common/avx256/MatrixComplexFloat2x2.h"

#include <functional>

namespace aocommon {

#ifndef __AVX2__
// AVX isn't present at build time.
//
// It's not possible to generate AVX code since the compiler issue an ABI issue.
//   warning: AVX vector return without AVX enabled changes the ABI [-Wpsabi]
// TODO See whether it's still possible to generate function multi versioning.
using MatrixComplexFloat2x2 = Scalar::MatrixComplexFloat2x2;

#elif !defined(PORTABLE_BUILD)  // __AVX2__

// AVX is present at build time and non-portable builds.
//
// Always use AVX.
using MatrixComplexFloat2x2 = Avx256::MatrixComplexFloat2x2;

#else  // __AVX2__

// AVX is present at build time and non-portable builds.
//
// Use function multiversioning to conditionally support AVX. This isn't as
// fast as always using AVX. Still when present it's still faster than the
// scalar code. Platforms that don't support AVX at runtime will run slower
// than before.

/**
 * MatrixComplexFloat2x2 using Function Multiversioning (FMV).
 *
 * This class can be used in environments with and without AVX. It uses the
 * compiler's FMV. The class contains comments explaining its usage in the
 * class itself.
 *
 * Based on benchmarking the FMV has little overhead, depending on how it's
 * used in the class. A different approach would be to duplicate all functions
 * and call the enum's @ref scalar_ or @ref avx_ element. This has two
 * disadvantages
 * - The number of functions grows hard and it has a lot of boiler-plate.
 * - The overhead of all the FMV dispatching is measurable. Tested with the
 *   multiply function:
 *   - Scalar 24 ns
 *   - FMV     7 ns
 *   - AVX     6 ns
 * The approach is still used for the constructors and @ref GetDispatch.
 */
class MatrixComplexFloat2x2 {
  /**
   * The dispatch method used.
   *
   * This can be extended in the future with other SIMD variants.
   */
  enum class Dispatch {
    /// Uses scalar code.
    Scalar,
    /// Uses AVX-256 and FMA3 instructions.
    Avx256
  };

  /**
   * This contains the "magic" to determine the active value of the union.
   *
   * Based on whether or not the CPU has AVX-256 and FMA3 support it will
   * return the appropriate dispatcher.
   *
   * @note When caching this value in a static class variable the code became
   * slower.
   */
  [[nodiscard]] [[gnu::target("default")]] static Dispatch GetDispatch() {
    return Dispatch::Scalar;
  }

  [[nodiscard]] [[gnu::target("avx2,fma")]] static Dispatch GetDispatch() {
    return Dispatch::Avx256;
  }

  /// Converting constructor using a scalar matrix.
  MatrixComplexFloat2x2(Scalar::MatrixComplexFloat2x2 scalar)
      : scalar_{scalar} {}

  /// Converting constructor using an AVX  matrix.
  MatrixComplexFloat2x2(Avx256::MatrixComplexFloat2x2 avx) : avx_{avx} {}

  /**
   * Executes a function call on the proper underlying matrix.
   *
   * This is used to convert a member functions to the proper member function of
   * the underlying matrix. The converted function can have additional
   * arguments that are forwarded. For example,
   * @code std::complex<double> operator[](size_t index); @endcode uses an
   * argument.
   *
   * The naming in the scalar matrix doesn't match Google style, the code in
   * the AVX matrix does. When the names differ some adjustments need to be
   * made. For example
@code
[[nodiscard]] MatrixComplexFloat2x2 HermitianTranspose() const noexcept {
  // Use generic lambda's to accept all underlying types.
  return ExecuteCall([](auto object)
         // The functions return their underlying matrix. These are different
         // types, specify the return type to get the proper type.
     -> MatrixComplexFloat2x2
  {
    // Is the underlying type the scalar matrix?
    if constexpr (std::is_same_v<decltype(object),
                                 Scalar::MatrixComplexFloat2x2>)
          // Call the scalar matrix' member function.
      return object.HermTranspose();
    else
          // Call the AVX matrix' member function.
      return object.HermitianTranspose();
  });
}
@endcode
   */
  template <class Function, class... Args>
  [[nodiscard]] auto ExecuteCall(Function function, Args&&... args) const {
    switch (GetDispatch()) {
      case Dispatch::Scalar:
        return function(scalar_, std::forward<Args>(args)...);
      case Dispatch::Avx256:
        return function(avx_, std::forward<Args>(args)...);
    }
    __builtin_unreachable();
  }

  /**
   * Executes an operation on two matrices.
   *
   * For example @c operator*, @c operator+, etc.
   */
  template <class BinaryOperator>
  [[nodiscard]] friend MatrixComplexFloat2x2 ExecuteBinaryOperator(
      MatrixComplexFloat2x2& lhs, MatrixComplexFloat2x2 rhs,
      BinaryOperator binary_operator) {
    switch (GetDispatch()) {
      case Dispatch::Scalar:
        return binary_operator(lhs.scalar_, rhs.scalar_);
      case Dispatch::Avx256:
        return binary_operator(lhs.avx_, rhs.avx_);
    }
    __builtin_unreachable();
  }

  /**
   * Handles a (compound) assignment.
   *
   * The function should take the first argument as a reference to it can be
   * modified. For example @c operator=, @c operator+=, etc.
   */
  template <class AssignmentOperator>
  [[nodiscard]] MatrixComplexFloat2x2& ExecuteAssign(
      MatrixComplexFloat2x2 rhs, AssignmentOperator assignment_operator) {
    switch (GetDispatch()) {
      case Dispatch::Scalar:
        assignment_operator(scalar_, rhs.scalar_);
        break;
      case Dispatch::Avx256:
        assignment_operator(avx_, rhs.avx_);
        break;
    }
    return *this;
  }

 public:
  //
  // Constructors
  //
  // These need to use FMV to store the data in the proper matrix.
  //
  [[gnu::target("default")]] explicit MatrixComplexFloat2x2(
      std::complex<double> a, std::complex<double> b, std::complex<double> c,
      std::complex<double> d) noexcept
      : scalar_{a, b, c, d} {}

  [[gnu::target("avx2,fma")]] explicit MatrixComplexFloat2x2(
      std::complex<double> a, std::complex<double> b, std::complex<double> c,
      std::complex<double> d) noexcept
      : avx_{a, b, c, d} {}

  [[gnu::target("default")]] explicit MatrixComplexFloat2x2(
      const std::complex<float> matrix[4]) noexcept
      : scalar_{matrix} {}

  [[gnu::target("avx2,fma")]] explicit MatrixComplexFloat2x2(
      const std::complex<float> matrix[4]) noexcept
      : avx_{matrix} {}

  [[gnu::target("default")]] explicit MatrixComplexFloat2x2(
      const std::complex<double> matrix[4]) noexcept
      : scalar_{matrix} {}

  [[gnu::target("avx2,fma")]] explicit MatrixComplexFloat2x2(
      const std::complex<double> matrix[4]) noexcept
      : avx_{matrix} {}

  [[nodiscard]] MatrixComplexFloat2x2& operator=(
      const MatrixComplexFloat2x2& rhs) {
    return ExecuteAssign(rhs, [](auto& lhs, auto rhs) { lhs = rhs; });
  }

  //
  // Operations
  //
  // The use the private helper functions to dispatch the call to the proper
  // underlying type.
  //

  [[nodiscard]] std::complex<double> operator[](size_t index) const noexcept {
    return ExecuteCall(
        [](const auto& object, size_t index) { return object[index]; }, index);
  }

  MatrixComplexFloat2x2& operator+=(MatrixComplexFloat2x2 rhs) noexcept {
    return ExecuteAssign(rhs, [](auto& lhs, auto rhs) { lhs += rhs; });
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator*(
      MatrixComplexFloat2x2 lhs, MatrixComplexFloat2x2 rhs) noexcept {
    return ExecuteBinaryOperator(lhs, rhs, std::multiplies{});
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator*(
      MatrixComplexFloat2x2 lhs, std::complex<double> rhs) const noexcept {
    return ExecuteBinaryOperator(lhs, rhs, std::multiplies{});
  }

  [[nodiscard]] friend MatrixComplexFloat2x2 operator*(
      std::complex<double> lhs, MatrixComplexFloat2x2 rhs) const noexcept {
    return ExecuteBinaryOperator(lhs, rhs, std::multiplies{});
  }

  [[nodiscard]] MatrixComplexFloat2x2 Transpose() const noexcept {
    return ExecuteCall([](auto object) -> MatrixComplexFloat2x2 {
      return object.Transpose();
    });
  }

  [[nodiscard]] MatrixComplexFloat2x2 HermitianTranspose() const noexcept {
    return ExecuteCall([](auto object) -> MatrixComplexFloat2x2 {
      if constexpr (std::is_same_v<decltype(object),
                                   Scalar::MatrixComplexFloat2x2>)
        return object.HermTranspose();
      else
        return object.HermitianTranspose();
    });
  }

  [[nodiscard]] MatrixComplexFloat2x2 Conjugate() const noexcept {
    return ExecuteCall([](auto object) -> MatrixComplexFloat2x2 {
      return object.Conjugate();
    });
  }

  [[nodiscard]] double Norm() const noexcept {
    return ExecuteCall([](const auto& object) {
      if constexpr (std::is_same_v<std::remove_cv_t<std::remove_reference_t<
                                       decltype(object)>>,
                                   Scalar::MatrixComplexFloat2x2>)
        return Norm(object);
      else
        return object.Norm();
    });
  }

  [[nodiscard]] std::complex<float> Trace() const noexcept {
    switch (GetDispatch()) {
      case Dispatch::Scalar:
        return Trace(*this);
      case Dispatch::Avx256:
        return Avx256::MatrixComplexFloat2x2::Trace();
    }
  }

  /** Return 2x2 complex identity matrix */
  [[nodiscard]] static MatrixComplexFloat2x2 Unity() noexcept {
    switch (GetDispatch()) {
      case Dispatch::Scalar:
        return Scalar::MatrixComplexFloat2x2::Unity();
      case Dispatch::Avx256:
        return Avx256::MatrixComplexFloat2x2::Unity();
    }
  }

  void AssignTo(std::complex<float>* destination) const noexcept {
    ExecuteCall([](auto object) -> MatrixComplexFloat2x2 {
      object.AssignTo(destination)();
    });
  }

  void AssignTo(std::complex<double>* destination) const noexcept {
    ExecuteCall([](auto object) -> MatrixComplexFloat2x2 {
      object.AssignTo(destination)();
    });
  }

  [[nodiscard]] DiagnoalMatrixComplexFloat2x2 Diagonal() const noexcept {
    switch (GetDispatch()) {
      case Dispatch::Scalar:
        return Diagonal(*this);
      case Dispatch::Avx256:
        return Avx256::MatrixComplexFloat2x2::Diagonal();
    }
  }

  //
  // Deprecated operations
  //
  // The are resembling operations but use names not conforming to Google
  // Style or use named operations instead of operator overloading.
  //
  // RAP-133 enabled diagnostic [[deprecated("Use HermitianTranspose")]]
  [[nodiscard]] MatrixComplexFloat2x2 HermTranspose() const noexcept {
    return HermitianTranspose();
  }

 private:
  union {
    Scalar::MatrixComplexFloat2x2 scalar_;
    Avx256::MatrixComplexFloat2x2 avx_;
  };
};

// RAP-133 enabled diagnostic
// [[deprecated("Use MatrixComplexFloat2x2::HermitianTranspose")]]
[[nodiscard]] inline MatrixComplexFloat2x2 HermTranspose(
    MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.HermitianTranspose();
}

[[nodiscard]] inline std::complex<float> Trace(
    MatrixComplexFloat2x2 matrix) noexcept {
  return matrix.Trace();
}

[[nodiscard]] inline DiagnoalMatrixComplexFloat2x2 Diagonal(
    MatrixComplexFloat2x2 matrix) const noexcept {
  return matrix.Diagnoal();
}

#endif  // __AVX2__

}  // namespace aocommon

#endif  // AOCOMMON_MATRIX_COMPLEX_FLOAT_2X2_H
