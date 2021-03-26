// BlobStreamable.h: Interface for classes that can be streamed using blobs.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBSTREAMABLE_H
#define LOFAR_BLOB_BLOBSTREAMABLE_H

#include <string>

namespace dp3 {
namespace blob {

class BlobIStream;
class BlobOStream;

/// \addtogroup Blob

/// @{

/// \brief Interface for classes that can be streamed using blobs.

/// Classes that want to be blob-streamable must implement the methods
/// declared in this interface.
class BlobStreamable {
 public:
  /// Destructor.
  virtual ~BlobStreamable() {}

  /// Create a new BlobStreamable object by deserializing the contents of
  /// the blob input stream \a bis. The first element in the input stream
  /// \a bis shall be a string containing the class type of the object to
  /// be constructed.
  static BlobStreamable* deserialize(BlobIStream& bis);

  /// Serialize \c *this by writing its contents into the blob output
  /// stream \a bos.
  void serialize(BlobOStream& bos) const;

 protected:
  /// Return the class type of \c *this as a string.
  virtual const std::string& classType() const = 0;

 private:
  /// Read the contents from the blob input stream \a bis into \c *this.
  virtual void read(BlobIStream& bis) = 0;

  /// Write the contents of \c *this into the blob output stream \a bos.
  virtual void write(BlobOStream& bos) const = 0;

  /// The private methods must be accessible for DH_BlobStreamable.
  friend class DH_BlobStreamable;
};

/// @}

}  // namespace blob
}  // namespace dp3

#endif
