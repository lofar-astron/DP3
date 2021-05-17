#include "RegularSolverBase.h"

namespace dp3 {
namespace ddecal {

void RegularSolverBase::Initialize(size_t n_antennas, size_t n_directions,
                                   size_t n_channels, size_t n_channel_blocks,
                                   const std::vector<int>& ant1,
                                   const std::vector<int>& ant2) {
  SolverBase::Initialize(n_antennas, n_directions, n_channel_blocks);
  n_channels_ = n_channels;
  ant1_ = ant1;
  ant2_ = ant2;
}

}  // namespace ddecal
}  // namespace dp3
