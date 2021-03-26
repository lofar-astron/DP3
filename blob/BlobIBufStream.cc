// BlobIBufStream.cc: Input buffer for a blob using an istream
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BlobIBufStream.h"

#include <iostream>

namespace dp3 {
namespace blob {

BlobIBufStream::BlobIBufStream(std::istream& is) : itsStream(is.rdbuf()) {}

BlobIBufStream::~BlobIBufStream() {}

uint64_t BlobIBufStream::get(void* buffer, uint64_t nbytes) {
  return itsStream->sgetn((char*)buffer, nbytes);
}

int64_t BlobIBufStream::tellPos() const {
  return itsStream->pubseekoff(0, std::ios::cur);
}

int64_t BlobIBufStream::setPos(int64_t pos) {
  return itsStream->pubseekoff(pos, std::ios::beg);
}

}  // namespace blob
}  // namespace dp3
