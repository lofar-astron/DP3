// BlobIBufStream.h: Input buffer for a blob using an istream
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBIBUFSTREAM_H
#define LOFAR_BLOB_BLOBIBUFSTREAM_H

#include "BlobIBuffer.h"

#include <iosfwd>

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief Input buffer for a blob using an istream

/// @{

/// This class is the BlobIBuffer that makes use of an istream object.
/// The istream can be any type (ifstream, istringstream, ...).
/// It can, for instance, be used to read from a file or a socket.

class BlobIBufStream : public BlobIBuffer {
 public:
  /// Construct it with the underlying istream object.
  explicit BlobIBufStream(std::istream&);

  /// Destructor.
  virtual ~BlobIBufStream();

  /// Get the requested nr of bytes.
  virtual uint64_t get(void* buffer, uint64_t nbytes);

  /// Get the position in the stream.
  /// -1 is returned if the stream is not seekable.
  virtual int64_t tellPos() const;

  /// Set the position in the stream.
  /// It returns the new position which is -1 if the stream is not seekable.
  virtual int64_t setPos(int64_t pos);

 private:
  std::streambuf* itsStream;
};

/// @}

}  // namespace blob
}  // namespace dp3

#endif
