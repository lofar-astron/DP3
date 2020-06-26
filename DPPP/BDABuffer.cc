#include "BDABuffer.h"

namespace DP3 {
  namespace DPPP {

    BDABuffer::Row::Row(const double time,
                        const double exposure,
                        const rownr_t rowNr,
                        const std::size_t nChannels,
                        const std::size_t nCorrelations,
                        std::complex<float>* const data,
                        const std::vector<bool>::iterator flags,
                        float* const weights,
                        const std::vector<bool>::iterator fullResFlags)
    : itsTime(time)
    , itsExposure(exposure)
    , itsRowNr(rowNr)
    , itsNChannels(nChannels)
    , itsNCorrelations(nCorrelations)
    , itsData(data)
    , itsFlags(flags)
    , itsWeights(weights)
    , itsFullResFlags(fullResFlags)
    , itsUVW{ 0.0, 0.0, 0.0 }
    {}

    BDABuffer::BDABuffer(const std::size_t poolSize)
        : itsData(), itsFlags(), itsWeights(), itsFullResFlags(), itsRows()
    {
      itsData.reserve(poolSize);
      itsFlags.reserve(poolSize);
      itsWeights.reserve(poolSize);
      itsFullResFlags.reserve(poolSize);
    }

    bool BDABuffer::addRow(const double time,
                           const double exposure,
                           const rownr_t rowNr,
                           const std::size_t nChannels,
                           const std::size_t nCorrelations,
                           const std::complex<double> *const data,
                           const bool *const flags,
                           const float *const weights,
                           const bool *const fullResFlags,
                           const double *const UVW)
    {
      const std::size_t dataSize = nChannels * nCorrelations;

      // Check if there is enough capacity left.
      if ((itsData.capacity() - itsData.size()) < dataSize) {
        return false;
      }

      itsRows.emplace_back(time,
                           exposure,
                           rowNr,
                           nChannels,
                           nCorrelations,
                           itsData.data() + itsData.size(),
                           itsFlags.begin() + itsFlags.size(),
                           itsWeights.data() + itsWeights.size(),
                           itsFullResFlags.begin() + itsFullResFlags.size());
      if (data) {
        itsData.insert(itsData.end(),
                       data,
                       data + dataSize);
      } else {
        itsData.insert(itsData.end(), dataSize, 0.0);
      }

      return true;
    }
  }
}
