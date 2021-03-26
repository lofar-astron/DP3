// BlobOBuffer.h: Abstract base class for output buffer for a blob
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBOBUFFER_H
#define LOFAR_BLOB_BLOBOBUFFER_H

#include <cinttypes>

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief Abstract base class for output buffer for a blob

/// @{

/// BlobOBuffer is the abstract base class for the sink of a
/// BlobOStream object. In this way the destination of a BlobOStream
/// can be chosen as needed. Currently three sink types are possible:
/// <ul>
/// <li> BlobOBufStream uses an ostream as the sink. It is mainly
///      meant to write blobs into a file, socket, etc.
/// <li> BlobOBufVector<T> uses a vector<T> (where T is char or uchar)
///      as the sink. It is mainly meant for creating blobs to be
///      stored using the classes in the LCS Persistent Object Layer.
/// <li> BlobOBufChar uses a plain C-array as the sink. It serves as
///      a helper class for BlobOBufVector, but can also be used in itself.
/// </ul>

class BlobOBuffer {
 public:
  BlobOBuffer() {}

  virtual ~BlobOBuffer() {}

  /// Put the requested nr of bytes.
  virtual uint64_t put(const void* buffer, uint64_t nbytes) = 0;

  /// Get the position in the stream.
  /// -1 is returned if the stream is not seekable.
  virtual int64_t tellPos() const = 0;

  /// Set the position in the stream.
  /// It returns the new position which is -1 if the stream is not seekable.
  virtual int64_t setPos(int64_t pos) = 0;
};

/// @}

}  // namespace blob
}  // namespace dp3

#endif
