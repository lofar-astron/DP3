#include <algorithm>
#include <cassert>
#include <complex>
#include <vector>

namespace detail {

template <typename T>
void FlagInvalid(const std::vector<T>& data, std::vector<bool>& flags) {
  assert(data.size() == flags.size());

  for (size_t i = 0; i < data.size(); ++i) {
    if (!std::isfinite(data[i])) {
      flags[i] = true;
    }
  }
}

}  // namespace detail

template <typename T>
T ComputeMean(const std::vector<T>& data, const std::vector<bool>& flags) {
  assert(data.size() == flags.size());

  T sum = 0;
  size_t unflagged_count = 0;

  for (size_t i = 0; i < data.size(); ++i) {
    if (!flags[i]) {
      sum += data[i];
      ++unflagged_count;
    }
  }

  return unflagged_count > 0 ? sum / unflagged_count : 0;
}

/*
 * Computes the median of the data with length n, ignoring all values
 * where the corresponding flag is set to True.
 * This code is based on:
 *  https://medium.com/@nxtchg/calculating-median-without-sorting-eaa639cedb9f
 * with a number of notable changes:
 *  - The input is an array with values of type T (instead of objects with value
 *    and weight)
 *   - Skips values for which its flag is true, or which are NaN or Inf
 *   - Handles even-length input as a special case
 */
template <typename T>
T ComputeMedianInplace(const std::vector<T>& data,
                       const std::vector<bool>& flags) {
  assert(data.size() == flags.size());

  T below = 0;
  T above = 0;
  T equal = ComputeMean(data, flags);

  while (true) {
    size_t n_below = 0;
    size_t n_equal = 0;
    size_t n_above = 0;

    for (size_t i = 0; i < data.size(); ++i) {
      T value = data[i];
      if (flags[i] || !IsValid(value)) {
        continue;
      }
      if (value > equal) {
        if (n_above == 0 || value < above) {
          above = value;
        }
        ++n_above;
      } else if (value < equal) {
        if (n_below == 0 || value > below) {
          below = value;
        }
        ++n_below;
      } else {
        ++n_equal;
      }
    }

    T mid = n_below + n_equal + n_above;

    if (mid > 0) {
      mid /= static_cast<T>(2);
    } else {
      return 0;
    }

    if (mid >= n_below) {
      if (mid <= n_below + n_equal) {
        if (data.size() % 2 == 1) {
          return equal;
        } else {
          if (mid == n_below) {
            for (size_t j = 0; j < data.size(); ++j) {
              if (flags[j]) {
                continue;
              }
              T value = data[j];
              if ((value > below) && (value <= above)) {
                above = value;
              }
            }
            return (below + above) / static_cast<T>(2);
          } else {
            return (equal + above) / static_cast<T>(2);
          }
        }
      }

      equal = above;
    } else {
      equal = below;
    }
  }
}

template <typename T>
T ComputeMedian(std::vector<T> v) {
  auto target = v.begin() + v.size() / static_cast<T>(2);
  std::nth_element(v.begin(), target, v.end());
  if (v.size() % 2 != 0) {
    return *target;
  } else {
    auto before = target - 1;
    std::nth_element(v.begin(), before, v.end());
    return (*target + *before) / static_cast<T>(2);
  }
}

template <typename T>
T ComputeSTD(const std::vector<T>& data, const std::vector<bool>& flags) {
  assert(data.size() == flags.size());

  T sum = 0;
  size_t m = 0;

  T mean = ComputeMean(data, flags);

  for (size_t i = 0; i < data.size(); ++i) {
    if (!flags[i]) {
      T value = mean - data[i];
      sum += value * value;
      ++m;
    }
  }

  return m > 0 ? std::sqrt(sum / m) : static_cast<T>(0);
}

template <typename T>
T ComputeSUMP2(const std::vector<T>& data, const std::vector<bool>& flags) {
  assert(data.size() == flags.size());

  T sum = 0;

  for (size_t i = 0; i < data.size(); ++i) {
    if (!flags[i]) {
      const T value = data[i];
      sum += value * value;
    }
  }

  return sum;
}

template <typename T>
std::vector<bool> FindOutliers(float sigma, int maxiters,
                               const std::vector<T>& data) {
  std::vector<bool> flags(data.size(), {false});
  detail::FlagInvalid(data, flags);

  for (int i = 0; i < maxiters; ++i) {
    const T median = ComputeMedian(data);
    const T std = ComputeSTD(data, flags);

    // Find outliers
    size_t outlier_count = 0;
    const T lower_bound = median - (sigma * std);
    const T upper_bound = median + (sigma * std);

    for (size_t j = 0; j < data.size(); ++j) {
      const T value = data[j];
      if (!flags[j] && (value < lower_bound || value > upper_bound)) {
        flags[j] = true;
        ++outlier_count;
      }
    }

    // Exit early when no new outliers were found
    if (outlier_count == 0) {
      break;
    }
  }

  return flags;
}