// BlobOBufStream.cc: Output buffer for a blob using an ostream
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BlobOBufStream.h"

#include <iostream>

namespace dp3 {
namespace blob {

BlobOBufStream::BlobOBufStream(std::ostream& os) : itsStream(os.rdbuf()) {}

BlobOBufStream::~BlobOBufStream() {}

uint64_t BlobOBufStream::put(const void* buffer, uint64_t nbytes) {
  return itsStream->sputn((const char*)buffer, nbytes);
}

int64_t BlobOBufStream::tellPos() const {
  return itsStream->pubseekoff(0, std::ios::cur);
}

int64_t BlobOBufStream::setPos(int64_t pos) {
  return itsStream->pubseekoff(pos, std::ios::beg);
}

}  // namespace blob
}  // namespace dp3
