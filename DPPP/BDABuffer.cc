#include "BDABuffer.h"

namespace DP3 {
  namespace DPPP {

    BDABuffer::Row::Row(const double time,
                        const double exposure,
                        const rownr_t rowNr,
                        const std::size_t nChannels,
                        const std::size_t nCorrelations,
                        const std::vector<std::complex<float>>::iterator data,
                        const std::vector<bool>::iterator flags,
                        const std::vector<float>::iterator weights,
                        const std::vector<bool>::iterator fullResFlags,
                        const double *const UVW)
    : itsTime(time)
    , itsExposure(exposure)
    , itsRowNr(rowNr)
    , itsNChannels(nChannels)
    , itsNCorrelations(nCorrelations)
    , itsData(data)
    , itsFlags(flags)
    , itsWeights(weights)
    , itsFullResFlags(fullResFlags)
    , itsUVW{ UVW ? UVW[0] : std::nan(""),
              UVW ? UVW[1] : std::nan(""),
              UVW ? UVW[2] : std::nan("") }
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
                           const std::complex<float>* const data,
                           const bool* const flags,
                           const float* const weights,
                           const bool* const fullResFlags,
                           const double* const UVW)
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
                           itsData.end(),
                           itsFlags.end(),
                           itsWeights.end(),
                           itsFullResFlags.end(),
                           UVW);
      if (data) {
        itsData.insert(itsData.end(), data, data + dataSize);
      } else {
        const std::complex<float> nan(std::nanf(""), std::nanf(""));
        itsData.insert(itsData.end(), dataSize, nan);
      }

      if (flags) {
        itsFlags.insert(itsFlags.end(), flags, flags + dataSize);
      } else {
        itsFlags.insert(itsFlags.end(), dataSize, false);
      }

      if (weights) {
        itsWeights.insert(itsWeights.end(), weights, weights + dataSize);
      } else {
        const float nan = std::nanf("");
        itsWeights.insert(itsWeights.end(), dataSize, nan);
      }

      if (fullResFlags) {
        itsFullResFlags.insert(itsFullResFlags.end(),
                               fullResFlags,
                               fullResFlags + dataSize);
      } else {
        itsFullResFlags.insert(itsFullResFlags.end(), dataSize, false);
      }

      return true;
    }
  }
}
