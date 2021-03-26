// BlobIBuffer.h: Abstract base class for input buffer for a blob
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cstdint>

#ifndef LOFAR_BLOB_BLOBIBUFFER_H
#define LOFAR_BLOB_BLOBIBUFFER_H

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief Abstract base class for input buffer for a blob

/// @{

/// BlobIBuffer is the abstract base class for the source of a
/// BlobIStream object. In this way the source of a BlobIStream
/// can be chosen as needed. Currently three source types are possible:
/// <ul>
/// <li> BlobIBufStream uses an istream as the source. It is mainly
///      meant to read blobs from a file, socket, etc.
/// <li> BlobIBufChar uses a plain C-array as the source. Its mainly used
///      as a helper class, but can also be used in itself.
/// <li> BlobIBufVector<T> uses a vector<T> (where T is char or uchar)
///      as the source.
///      It uses BlobIBufChar to do the actual work.
/// <li> BlobIBufString uses a BlobString as the source.
///      It uses BlobIBufChar to do the actual work.
///      It is mainly meant for reading blobs stored using the LCS
///      Persistency Layer.
/// </ul>

class BlobIBuffer {
 public:
  /// Constructor.
  BlobIBuffer() {}

  /// Destructor.
  virtual ~BlobIBuffer() {}

  /// Get the requested nr of bytes.
  virtual uint64_t get(void* buffer, uint64_t nbytes) = 0;

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
