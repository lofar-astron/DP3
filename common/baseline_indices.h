/*
 * File taken from Aartfaac2ms repository:
 * https://git.astron.nl/RD/aartfaac-tools/-/blob/master/common/baseline_indices.h
 */

#ifndef DP3_COMMON_BASELINE_INDICES_H_
#define DP3_COMMON_BASELINE_INDICES_H_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace dp3 {
namespace common {

/*
 * Baselines are typically ordered in either column-major (e.g. for LOFAR), or
 * in row-major order (e.g. for AARTFAAC). In both cases, for any baseline
 * (antenna1,antenna2): antenna1 < antenna2 (i.e. there is a baseline (0,1), but
 * not its mirror (1,0)).
 *
 * This is an example of column-major ordering with
 * 5 antennas in the format "index:(antenna1,antenna2)":
 *   0:(0,0)
 *   1:(0,1)  2:(1,1)
 *   3:(0,2)  4:(1,2)  5:(2,2)
 *   6:(0,3)  7:(1,3)  8:(2,3)  9:(3,3)
 *  10:(0,4) 11:(1,4) 12:(2,4) 13:(3,4) 14:(4,4)
 *  The first column has all the baselines with antenna1 == 0
 *
 * And similarly for the row-major ordering:
 *   0:(0,0)  1:(0,1)  2:(0,2)  3:(0,3)  4:(0,4)
 *   5:(1,1)  6:(1,2)  7:(1,3)  8:(1,4)
 *   9:(2,2) 10:(2,3) 11:(2,4)
 *  12:(3,3) 13:(3,4)
 *  14:(4,4)
 * The first row has all the baselines with antenna1 == 0
 *
 * Next, we list all baselines that contain a given antenna.
 * The format is as follows: "antenna: [indices] [baselines]".
 *
 * First, column-major ordering:
 * 0: [ 0,  1,  3,  6, 10] [(0,0), (0,1), (0,2), (0,3), (0,4)]
 * 1: [ 1,  2,  4,  7, 11] [(0,1), (1,1), (1,2), (1,3), (1,4)]
 * 2: [ 3,  4,  5,  8, 12] [(0,2), (1,2), (2,2), (2,3), (2,4)]
 * 3: [ 6,  7,  8,  9, 13] [(0,3), (1,3), (2,3), (3,3), (3,4)]
 * 4: [10, 11, 12, 13, 14] [(0,4), (1,4), (2,4), (3,4), (4,4)]
 *
 * Next, row-major ordering:
 * 0: [0, 1,  2,  3,  4] [(0,0], (0,1), (0,2), (0,3), (0,4)]
 * 1: [1, 5,  6,  7,  8] [(0,1], (1,1), (1,2), (1,3), (1,4)]
 * 2: [2, 6,  9, 10, 11] [(0,2], (1,2), (2,2), (2,3), (2,4)]
 * 3: [3, 7, 10, 12, 13] [(0,3], (1,3), (2,3), (3,3), (3,4)]
 * 4: [4, 8, 11, 13, 14] [(0,4], (1,4), (2,4), (3,4), (4,4)]
 *
 * This file contains several helper functions to go from index
 * to baseline, from baseline to index or to get all the baselines
 * for a given antenna.
 */

enum class BaselineOrder { kColumnMajor, kRowMajor };

/*
 * Computes the number of baselines, including auto-correlations.
 */
inline size_t ComputeNBaselines(size_t n_antennas) {
  const size_t n_correlations = (n_antennas * (n_antennas - 1)) / 2;
  const size_t n_autocorrelations = n_antennas;
  return n_correlations + n_autocorrelations;
}

/*
 * Computes the baseline index corresponding to antenna a and b, given the
 * number of antennas. It assumes that auto-correlations of antennas are
 * included.
 */
inline size_t ComputeBaselineIndex(size_t antenna_a, size_t antenna_b,
                                   size_t nAntennas, BaselineOrder order) {
  size_t row = std::min(antenna_a, antenna_b);
  size_t col = std::max(antenna_a, antenna_b);
  if (order == BaselineOrder::kRowMajor) {
    return size_t((row * nAntennas) + col - row * (row + 1) / 2);
  } else {
    return size_t((col * (col + 1)) / 2) + row;
  }
}

/*
 * Computes the antenna pair (baseline) given the baseline index. It assumes
 * that auto-correlations of antennas are included.
 */
inline std::pair<size_t, size_t> ComputeBaseline(size_t baseline_index,
                                                 size_t n_antennas,
                                                 BaselineOrder order) {
  if (order == BaselineOrder::kRowMajor) {
    const size_t n_baselines = ComputeNBaselines(n_antennas);
    const size_t antenna =
        (1 + std::sqrt(1 + 8 * (n_baselines - baseline_index - 1))) / 2;
    const size_t antenna1 = n_antennas - antenna;
    // n is the number of baselines (a,b) with a < antenna
    const size_t n = n_baselines - n_antennas - (antenna * (antenna - 1)) / 2;
    const size_t antenna2 = baseline_index - n;
    return {antenna1, antenna2};
  } else {
    const size_t antenna2 = (1 + std::sqrt(1 + 8 * baseline_index)) / 2;
    // n is the number of baselines (a,b) with b < antenna2
    const size_t n = int((antenna2 * (antenna2 - 1)) / 2);
    const size_t antenna1 = baseline_index - n;
    return {antenna1, antenna2 - 1};
  }
}

/*
 * Computes the list of baselines for a given antenna. It assumes that
 * auto-correlations are included.
 */
inline std::vector<size_t> ComputeBaselineList(size_t antenna,
                                               size_t n_antennas,
                                               BaselineOrder order) {
  std::vector<size_t> indices;
  indices.reserve(n_antennas);

  for (size_t i = 0; i < n_antennas; i++) {
    indices.push_back(ComputeBaselineIndex(antenna, i, n_antennas, order));
  }

  return indices;
}

}  // namespace common
}  // namespace dp3

#endif  // DP3_COMMON_BASELINE_INDICES_H_
