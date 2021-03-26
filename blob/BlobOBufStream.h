// BlobOBufStream.h: Output buffer for a blob using an ostream
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBOBUFSTREAM_H
#define LOFAR_BLOB_BLOBOBUFSTREAM_H

#include "BlobOBuffer.h"

#include <iosfwd>

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief Output buffer for a blob using an ostream

/// @{

/// This class is the BlobOBuffer that makes use of an ostream object.
/// The ostream can be any type (ofstream, ostringstream, ...)

class BlobOBufStream : public BlobOBuffer {
 public:
  /// Construct it with the underlying ostream object.
  explicit BlobOBufStream(std::ostream&);

  /// Destructor.
  virtual ~BlobOBufStream();

  /// Put the requested nr of bytes.
  virtual uint64_t put(const void* buffer, uint64_t nbytes);

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
