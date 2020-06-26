#include "BDABuffer.h"

namespace DP3 {
  namespace DPPP {

    BDABuffer::BDABuffer(std::size_t poolSize)
    {
    }

    bool BDABuffer::addLine(double time,
                            double exposure,
                            std::size_t nChannels,
                            std::size_t nCorrelations,
                            const std::complex<double> *measurements)
    {
      return false;
    }
  }
}
