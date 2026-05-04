#include "Constraint.h"

#include <xtensor/views/xview.hpp>

namespace dp3::ddecal {

void Constraint::ApplyReferenceAntenna(SolutionSpan& solutions) {
  const size_t n_polarizations = solutions.shape()[3];
  for (size_t pol_index = 0; pol_index != n_polarizations; ++pol_index) {
    // Choose reference antenna that has at least 20% channels unflagged
    size_t ref_antenna = 0;
    for (; ref_antenna != NAntennas(); ++ref_antenna) {
      // Only check flagged state for first direction
      const size_t n_unflagged_channels = xt::sum(xt::isfinite(
          xt::view(solutions, xt::all(), ref_antenna, 0, pol_index)))();
      if (n_unflagged_channels > 0.2 * NChannelBlocks())
        // Choose this refAntenna;
        break;
    }
    // All antennas are flagged, use first one (will lead to NaNs for this
    // solint)
    if (ref_antenna == NAntennas()) ref_antenna = 0;

    for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
      auto reference_view =
          xt::view(solutions, ch, ref_antenna, xt::all(), pol_index);
      for (size_t antenna_index = 0; antenna_index != NAntennas();
           ++antenna_index) {
        if (antenna_index != ref_antenna) {
          xt::view(solutions, ch, antenna_index, xt::all(), pol_index) /=
              reference_view;
        }
      }
      reference_view.fill(1.0);
    }
  }
}

}  // namespace dp3::ddecal
