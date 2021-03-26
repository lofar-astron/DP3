// BlobArray.cc: Blob handling for arrays
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BlobArray.h"
#include "BlobHeader.h"

#include <vector>

namespace dp3 {
namespace blob {

uint64_t putBlobArrayHeader(BlobOStream& bs, bool useBlobHeader,
                            const std::string& headerName,
                            const uint64_t* shape, uint16_t ndim,
                            bool fortranOrder, unsigned int alignment) {
  if (useBlobHeader) {
    bs.putStart(headerName, 1);  // version 1
  }
  unsigned char nalign = 0;
  if (alignment > 1) {
    int64_t pos = bs.tellPos();
    if (pos > 0) {
      nalign = (pos + 4 + ndim * sizeof(uint64_t)) % alignment;
      if (nalign > 0) {
        nalign = alignment - nalign;
      }
    }
  }
  bs << fortranOrder << nalign << ndim;
  bs.put(shape, ndim);
  uint64_t n = (ndim == 0 ? 0 : 1);
  for (int i = 0; i < ndim; i++) {
    n *= shape[i];
  }
  if (nalign > 0) {
    bs.put("        ", nalign);
  }
  return n;
}

// Get the shape of an array from the blob.
// This is a helper function for the functions reading an array.
uint64_t getBlobArrayShape(BlobIStream& bs, uint64_t* shape, unsigned int ndim,
                           bool swapAxes, unsigned int nalign) {
  bs.get(shape, ndim);
  if (swapAxes) {
    std::vector<uint64_t> shp(shape, shape + ndim);
    for (unsigned int i = 0; i < ndim; i++) {
      shape[i] = shp[ndim - i - 1];
    }
  }
  unsigned int n = 1;
  for (unsigned int i = 0; i < ndim; i++) {
    n *= shape[i];
  }
  char buf[32];
  while (nalign > 0) {
    int nb = std::min(32u, nalign);
    bs.get(buf, nb);
    nalign -= nb;
  }
  return n;
}

void convertArrayHeader(common::DataFormat fmt, char* header,
                        bool useBlobHeader) {
  char* buf = header;
  if (useBlobHeader) {
    BlobHeader* hdr = (BlobHeader*)header;
    hdr->setLocalDataFormat();
    buf += hdr->getHeaderLength();
  }
  // Skip the first 2 characters that do not need to be converted.
  buf += +2;
  // Get ndim and convert it in the buffer.
  uint16_t ndim;
  dataConvert16(fmt, &ndim, buf);
  dataConvert16(fmt, buf);
  buf += 2;
  // Convert all dimensions.
  dataConvert32(fmt, buf, ndim);
}

BlobOStream& operator<<(BlobOStream& bs, const std::vector<bool>& vec) {
  uint64_t size = vec.size();
  putBlobArrayHeader(bs, true, common::typeName((const bool**)0), &size, 1,
                     true, 1);
  bs.putBoolVec(vec);
  bs.putEnd();
  return bs;
}

BlobIStream& operator>>(BlobIStream& bs, std::vector<bool>& vec) {
  bs.getStart(common::typeName((const bool**)0));
  bool fortranOrder;
  uint16_t ndim;
  unsigned int nalign = getBlobArrayStart(bs, fortranOrder, ndim);
  assert(ndim == 1);
  uint64_t size;
  getBlobArrayShape(bs, &size, 1, false, nalign);
  bs.getBoolVec(vec, size);
  return bs;
}

#if defined(HAVE_AIPSPP)
BlobOStream& operator<<(BlobOStream& bs, const casacore::IPosition& ipos) {
  uint64_t size = ipos.size();
  putBlobArrayHeader(bs, true, DP3::typeName((const int64_t**)0), &size, 1,
                     true, 1);
  std::vector<int64_t> ivec(ipos.begin(), ipos.end());
  bs.put(&(ivec[0]), size);
  bs.putEnd();
  return bs;
}

BlobIStream& operator>>(BlobIStream& bs, casacore::IPosition& ipos) {
  bs.getStart(DP3::typeName((const int64_t**)0));
  bool fortranOrder;
  uint16_t ndim;
  unsigned int nalign = getBlobArrayStart(bs, fortranOrder, ndim);
  ASSERT(ndim == 1);
  uint64_t size;
  getBlobArrayShape(bs, &size, 1, false, nalign);
  std::vector<int64_t> ivec(size);
  ipos.resize(size, false);
  bs.get(&(ivec[0]), size);
  for (unsigned int i = 0; i < size; ++i) {
    ipos[i] = ivec[i];
  }
  bs.getEnd();
  return bs;
}
#endif

}  // namespace blob
}  // namespace dp3
