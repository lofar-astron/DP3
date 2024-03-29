// BlobOStream.h: Output stream for a blob
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBOSTREAM_H
#define LOFAR_BLOB_BLOBOSTREAM_H

#include "BlobOBuffer.h"

#include <stack>
#include <vector>
#include <string>
#include <complex>
#include <cstring>

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief Output stream for a blob

/// @{

/// This class makes it possible to create a blob.
/// It creates a header (in the putStart function) using the
/// BlobHeader definition.
/// The user can define overloaded operator<< functions to be able
/// to put objects of a given class into the blob stream.
///
/// The blob is written into a BlobOBuffer object that can be a memory
/// buffer or an ostream object. The BlobIStream class can be used to
/// retrieve objects from a blob.
///
/// See LOFAR document
/// <a
/// href="http://www.lofar.org/forum/document.php?action=match&docname=LOFAR-ASTRON-MAN-006">
/// LOFAR-ASTRON-MAN-006</a> for more information.

class BlobOStream {
 public:
  /// Construct it with the underlying buffer object.
  /// It keeps a pointer to the buffer, so be sure that the BlobOBuffer
  /// is not deleted before this object.
  explicit BlobOStream(BlobOBuffer&);

  /// Destructor.
  ~BlobOStream();

  /// Clear the object. I.e., reset the current level and length.
  void clear();

  /// Get the total size.
  uint64_t size() const { return itsCurLength; }

  /// Start putting a blob. It writes the header containing data that are
  /// checked when reading the blob back in BlobIStream::getStart.
  /// Data in nested objects can be put without an intermediate putStart.
  /// However, for complex objects it is recommended to do a putStart
  /// to have a better checking.
  /// <br>
  /// After all values (including nested objects) of the object have
  /// been put, a call to putEnd has to be done.
  //
  /// If no or an empty objecttype is given, the header is
  /// written without the objecttype and the associated length fields.
  /// In that case the function getStart in BlobIStream should also be
  /// called that way.
  //
  /// The function returns the nesting level.
  /// @{
  unsigned int putStart(int objectVersion);
  unsigned int putStart(const std::string& objectType, int objectVersion);
  unsigned int putStart(const char* objectType, int objectVersion);
  /// @}

  /// End putting an object. It returns the object length (including
  /// possible nested objects).
  uint64_t putEnd();

  /// Put a single value.
  /// A string will be stored as a length followed by the characters.
  /// @{
  BlobOStream& operator<<(const bool& value);
  BlobOStream& operator<<(const char& value);
  BlobOStream& operator<<(const int8_t& value);
  BlobOStream& operator<<(const uint8_t& value);
  BlobOStream& operator<<(const int16_t& value);
  BlobOStream& operator<<(const uint16_t& value);
  BlobOStream& operator<<(const int32_t& value);
  BlobOStream& operator<<(const uint32_t& value);
  BlobOStream& operator<<(const int64_t& value);
  BlobOStream& operator<<(const uint64_t& value);
  BlobOStream& operator<<(const float& value);
  BlobOStream& operator<<(const double& value);
  BlobOStream& operator<<(const std::complex<float>& value);
  BlobOStream& operator<<(const std::complex<double>& value);
  BlobOStream& operator<<(const std::string& value);
  BlobOStream& operator<<(const char* value);
  /// @}

  /// Put an array of values with the given number of values.
  /// Bool values are stored as bits.
  /// @{
  void put(const bool* values, uint64_t nrval);
  void put(const char* values, uint64_t nrval);
  void put(const int8_t* values, uint64_t nrval);
  void put(const uint8_t* values, uint64_t nrval);
  void put(const int16_t* values, uint64_t nrval);
  void put(const uint16_t* values, uint64_t nrval);
  void put(const int32_t* values, uint64_t nrval);
  void put(const uint32_t* values, uint64_t nrval);
  void put(const int64_t* values, uint64_t nrval);
  void put(const uint64_t* values, uint64_t nrval);
  void put(const float* values, uint64_t nrval);
  void put(const double* values, uint64_t nrval);
  void put(const std::complex<float>* values, uint64_t nrval);
  void put(const std::complex<double>* values, uint64_t nrval);
  void put(const std::string* values, uint64_t nrval);
  /// @}

  /// Put a vector of values. First it puts the size of the vector.
  /// Specialise for bool because a vector of bools uses bits.
  /// @{
  void put(const std::vector<bool>& values);
  template <class T>
  void put(const std::vector<T>& values);
  /// @}

  /// Put a vector of bools (without putting the size).
  void putBoolVec(const std::vector<bool>& values);

  /// Reserve the given amount of space (the opposite of
  /// BlobIStream::getSpace).
  /// This is useful when creating a static blob in a dynamic way.
  /// It returns the position of the skipped space in the stream.
  /// It is meant for use with the BlobOBufString buffer. The function
  /// getPointer in that class (in fact, in its base class BlobOBufChar)
  /// can be used to turn the position into a pointer.
  int64_t setSpace(uint64_t nbytes);

  /// Add filler bytes as needed to make the total length a multiple of n.
  /// In this way the next data are aligned properly.
  /// It returns the number of filler bytes used.
  /// It is only useful for seekable buffers.
  unsigned int align(unsigned int n);

  /// Get the current stream position.
  /// It returns -1 if the stream is not seekable.
  int64_t tellPos() const;

 private:
  /// Function to do the actual putStart.
  unsigned int doPutStart(const char* objectType, unsigned int nrc,
                          int objectVersion);

  /// Write the buffer, increment itsCurLength, and check
  /// if everything is written.
  void putBuf(const void* buf, uint64_t sz);

  /// Throw an exception if a put cannot be done.
  /// @{
  void checkPut() const;
  void throwPut() const;
  /// @}

  bool itsSeekable;
  uint64_t itsCurLength;
  unsigned int itsLevel;
  /// Object length at each level
  std::stack<uint64_t> itsObjLen;
  /// Offset of length at each level
  std::stack<int64_t> itsObjPtr;
  /// The underlying stream object.
  BlobOBuffer* itsStream;
};

/// @}

inline void BlobOStream::clear() {
  itsCurLength = 0;
  itsLevel = 0;
}

inline unsigned int BlobOStream::putStart(int objectVersion) {
  return doPutStart("", 0, objectVersion);
}

inline unsigned int BlobOStream::putStart(const std::string& objectType,
                                          int objectVersion) {
  return doPutStart(objectType.data(), objectType.size(), objectVersion);
}

inline unsigned int BlobOStream::putStart(const char* objectType,
                                          int objectVersion) {
  return doPutStart(objectType, strlen(objectType), objectVersion);
}

inline int64_t BlobOStream::tellPos() const { return itsStream->tellPos(); }

template <class T>
inline void BlobOStream::put(const std::vector<T>& vec) {
  operator<<(uint64_t(vec.size()));
  put(&vec[0], vec.size());
}
inline void BlobOStream::put(const std::vector<bool>& vec) {
  operator<<(uint64_t(vec.size()));
  putBoolVec(vec);
}

inline void BlobOStream::checkPut() const {
  if (itsLevel == 0) throwPut();
}

}  // namespace blob
}  // namespace dp3

#endif
