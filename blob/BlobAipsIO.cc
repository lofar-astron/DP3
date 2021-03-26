// BlobAipsIO.cc: A Blob buffer for Aips++ ByteIO
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BlobAipsIO.h"

using namespace casacore;

namespace dp3 {
namespace blob {

BlobAipsIO::BlobAipsIO(BlobOStream& os) : itsOBuf(&os), itsIBuf(0) {
  itsOBuf->putStart("BlobAipsIO", 1);
}

BlobAipsIO::BlobAipsIO(BlobIStream& is) : itsOBuf(0), itsIBuf(&is) {
  itsIBuf->getStart("BlobAipsIO");
}

BlobAipsIO::~BlobAipsIO() {
  if (itsOBuf) {
    itsOBuf->putEnd();
  } else {
    itsIBuf->getEnd();
  }
}

void BlobAipsIO::write(Int64 size, const void* buf) {
  itsOBuf->put(static_cast<const unsigned char*>(buf), size);
}

void BlobAipsIO::write(uInt size, const void* buf) {
  itsOBuf->put(static_cast<const unsigned char*>(buf), size);
}

Int64 BlobAipsIO::read(Int64 size, void* buf, Bool) {
  itsIBuf->get(static_cast<unsigned char*>(buf), size);
  return size;
}

Int BlobAipsIO::read(uInt size, void* buf, Bool) {
  itsIBuf->get(static_cast<unsigned char*>(buf), size);
  return size;
}

Int64 BlobAipsIO::length() { return -1; }

Bool BlobAipsIO::isReadable() const { return itsIBuf != 0; }

Bool BlobAipsIO::isWritable() const { return itsOBuf != 0; }

Bool BlobAipsIO::isSeekable() const { return false; }

Int64 BlobAipsIO::doSeek(Int64, ByteIO::SeekOption) { return 0; }

}  // namespace blob
}  // namespace dp3
