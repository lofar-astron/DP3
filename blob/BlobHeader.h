// BlobHeader.h: Standard header for a blob
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBHEADER_H
#define LOFAR_BLOB_BLOBHEADER_H

#include "../common/DataConvert.h"

#include <cstring>

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief  Standard header for a blob.

/// @{

/// Blob stands for binary large object.
/// The %LOFAR Common software provides classes to serialize one or more
/// objects into a blob and to de-serialize the blob to objects.
/// To be sure that a blob is interpreted in the correct way, each object
/// in it will be preceeded by a header. The header contains the
/// following information:
/// <ul>
/// <li> A magic value is used to indicate the start of an object.
/// <li> The length defines the total size of the object in the blob
///   (including the header). It may not always be possible to store
///   the length. In that case the length is 0.
/// <li> The version defines the object version. Because blobs can be
///   persistent (e.g. in a database), it makes it possible to handle
///   older instances of a class.
/// <li> The name defines the object type. In principle it is the name
///   of the class. Functions in TypeNames.h can be used to generate
///   the name of a templated class.
/// <li> The name length gives the actual length of the name.
/// <li> The reserved name length can be somewhat more in order to
///   achieve that the data thereafter is aligned on 8-byte boundary.
/// </ul>
/// This class is meant for handling blobs in a static and in a dynamic way.
/// Handling blobs in a dynamic way is done by means of the BlobOStream and
/// BlobIStream classes.
/// A blob can be created statically by putting a BlobHeader object ahead
/// of the data. This mode will be used in the CEPFrame environment.
/// For this mode the class is templated with the length of the name as the
/// template parameter.
/// For example:
/// \code
/// class SomeClass {
///   BlobHeader<9> itsHeader;
///   dcomplex      itsData[10][20];
/// };
///
/// SomeClass::SomeClass() : itsHeader("SomeClass", 1)
///    { itsHeader.setLength (sizeof(SomeClass)); }
/// \endcode
///
/// Because blobs can be created on one machine and retrieved on another,
/// care has to be taken that data type sizes and alignment are the same
/// everywhere. For this reason standard data types are defined in
/// LofarTypedefs.h. Data in a CEPFrame DataPacket should be declared such
/// that longer data types are declared first. In this way alignment should
/// never be a problem.
/// Care has been taken that a BlobHeader object does not disturb alignment.
/// So it is always aligned on a double boundary and its length is always
/// a multiple of 8.

class BlobHeader {
  friend class BlobOStream;
  friend class BlobIStream;

 public:
  /// Construct for the given name and version.
  BlobHeader(int version = 0, unsigned int level = 0);

  /// Get the data format.
  dp3::common::DataFormat getDataFormat() const {
    return dp3::common::DataFormat(itsDataFormat);
  }

  /// Set the data format to local. Useful after data is converted in place.
  void setLocalDataFormat();

  /// Get the version. Data will be converted if needed.
  int getVersion() const {
    return (mustConvert()
                ? dp3::common::dataConvert(getDataFormat(), itsVersion)
                : itsVersion);
  }

  /// Get the length of the blob. Data will be converted if needed.
  uint64_t getLength() const {
    return (mustConvert() ? dp3::common::dataConvert(getDataFormat(), itsLength)
                          : itsLength);
  }

  /// Set the length of the blob.
  void setLength(uint64_t length) { itsLength = length; }

  /// Test if the data format in the header mismatches the data format of
  /// this machine, thus if data have to be converted.
  bool mustConvert() const {
    return itsDataFormat != dp3::common::dataFormat();
  }

  /// Get the offset of the length.
  unsigned int lengthOffset() const { return 0; }

  /// Get the name length.
  unsigned int getNameLength() const { return itsNameLength; }

  /// Get the begin-of-blob magic value.
  static uint32_t bobMagicValue() { return 0xbebebebe; }

  /// Get the end-of-blob magic value.
  static uint32_t eobMagicValue() { return 0xbfbfbfbf; }

  bool checkMagicValue() const { return itsMagicValue == bobMagicValue(); }

  /// Get the length of the total header (thus including objecttype).
  unsigned int getHeaderLength() const {
    return sizeof(BlobHeader) + itsNameLength;
  }

 private:
  uint64_t itsLength;
  uint32_t itsMagicValue;
  char itsVersion;
  char itsDataFormat;
  unsigned char itsLevel;
  unsigned char itsNameLength;
};

/// @}

}  // namespace blob
}  // namespace dp3

#endif
