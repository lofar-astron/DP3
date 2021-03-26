// BlobAipsIO.h: A Blob buffer for Aips++ ByteIO
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBAIPSIO_H
#define LOFAR_BLOB_BLOBAIPSIO_H

#include <casacore/casa/IO/ByteIO.h>

#include "BlobOStream.h"
#include "BlobIStream.h"

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief A Blob buffer for Aips++ ByteIO

/// @{

/**
 * This class makes it possible to use the AIPS++ AipsIO class on top of
 * a Blob buffer. In this way it is possible to write an AIPS++ LSQFit
 * object directly into a Blob buffer.
 * Note that it is even possible to mix AipsIO and Blob puts;
 * the same is true for gets.
 * Of course, the order of reading must be the same as the order of writing.
 */
class BlobAipsIO : public casacore::ByteIO {
 public:
  /// Construct from a Blob buffer.
  /// @{
  explicit BlobAipsIO(BlobOStream&);
  explicit BlobAipsIO(BlobIStream&);
  /// @}

  virtual ~BlobAipsIO();

  /**
   * Write the number of bytes to the Blob stream.
   * \throw An exception if the stream is not writable.
   * The 2nd write and read function are defined for the read/write signatures
   * in older casacore versions (pre v2.0).
   */
  virtual void write(casacore::Int64 size, const void* buf);
  virtual void write(casacore::uInt size, const void* buf);

  /**
   * Read \a size bytes from the Blob stream. Returns the number of
   * bytes actually read.
   * \throw AipsError if the requested number of bytes could not be read,
   *        unless \a throwException is set to False.
   */
  virtual casacore::Int64 read(casacore::Int64 size, void* buf,
                               bool throwException);
  virtual casacore::Int read(casacore::uInt size, void* buf,
                             bool throwException);

  /// Get the length of the Blob stream. Returns 0 if the stream is not
  /// seekable.
  virtual casacore::Int64 length();

  /// Is the Blob stream readable?
  virtual bool isReadable() const;

  /// Is the Blob stream writable?
  virtual bool isWritable() const;

  /// Is the Blob stream seekable?
  virtual bool isSeekable() const;

 private:
  /// Make copy constructor and assignment private, so a user cannot
  /// use them.
  /// @{
  BlobAipsIO(const BlobAipsIO&);
  BlobAipsIO& operator=(const BlobAipsIO&);
  /// @}

  /// Reset the position pointer to the given value. It returns the
  /// new position.
  virtual casacore::Int64 doSeek(casacore::Int64 offset,
                                 casacore::ByteIO::SeekOption);

  BlobOStream* itsOBuf;
  BlobIStream* itsIBuf;
};

/// @}

}  // namespace blob
}  // namespace dp3

#endif
