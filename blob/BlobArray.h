// BlobArray.h: Blob handling for arrays
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBARRAY_H
#define LOFAR_BLOB_BLOBARRAY_H

#include "BlobOStream.h"
#include "BlobIStream.h"
#include <vector>
#if defined(HAVE_BLITZ)
#include <blitz/array.h>
#endif
#if defined(HAVE_AIPSPP)
#include <casacore/casa/Arrays/Array.h>
#endif

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief Blob handling for arrays

/// \addtogroup BlobArray Global BlobArray functions
/// Define functions to write N-dimensional arrays into a blob and to
/// read them back from a blob.
/// The arrays can be:
/// <ul>
/// <li> A plain N-dimensional C-array.
/// <li> A blitz array in Fortran order (minor axis first) or in C order.
/// <li> An AIPS++ array (which is in Fortran order).
/// </ul>
/// Special functions exist to read or write a vector (1-dimensional array).
/// Because all array types are written in a standard way, it is possible
/// to write, for example, a blitz array and read it back as an AIPS++ array.
/// If the axes ordering is different, the axes are reversed.
///
/// The write functions follow the same standard as the static array header
/// defined in BlobArrayHeader.h, so it is possible to read a static array
/// back in a dynamic way.
/// @{

/// \name General function to write a data array
/// Usually it is used by the other functions, but it can be used on
/// its own to write, say, a C-style array.
/// A 1-dim C-style array can be written with putBlobVector.
/// @{
template <typename T>
BlobOStream& putBlobArray(BlobOStream& bs, const T* data, const uint64_t* shape,
                          uint16_t ndim, bool fortranOrder);
template <typename T>
BlobOStream& putBlobVector(BlobOStream& bs, const T* data, uint64_t size);
/// @}

/// \name Reserve space for an array with the given shape
/// The axes ordering (Fortran or C-style) has to be given.
/// It returns the offset of the array in the blob.
/// It is useful for allocating a static blob in a dynamic way.
/// It is only possible if the underlying buffer is seekable.
/// It is meant for use with the BlobOBufString buffer. The function
/// getPointer in that class can be used to turn the position into a pointer.
/// The data will be aligned on the given alignment which defaults to
/// sizeof(T) bytes with a maximum of 8.
/// @{
template <typename T>
uint64_t setSpaceBlobArray1(BlobOStream& bs, bool useBlobHeader, uint64_t size0,
                            unsigned int alignment = 0);
template <typename T>
uint64_t setSpaceBlobArray2(BlobOStream& bs, bool useBlobHeader, uint64_t size0,
                            uint64_t size1, bool fortranOrder,
                            unsigned int alignment = 0);
template <typename T>
uint64_t setSpaceBlobArray3(BlobOStream& bs, bool useBlobHeader, uint64_t size0,
                            uint64_t size1, uint64_t size2, bool fortranOrder,
                            unsigned int alignment = 0);
template <typename T>
uint64_t setSpaceBlobArray4(BlobOStream& bs, bool useBlobHeader, uint64_t size0,
                            uint64_t size1, uint64_t size2, uint64_t size3,
                            bool fortranOrder, unsigned int alignment = 0);
template <typename T>
uint64_t setSpaceBlobArray(BlobOStream& bs, bool useBlobHeader,
                           const std::vector<uint64_t>& shape,
                           bool fortranOrder, unsigned int alignment = 0);
template <typename T>
uint64_t setSpaceBlobArray(BlobOStream& bs, bool useBlobHeader,
                           const uint64_t* shape, uint16_t ndim,
                           bool fortranOrder, unsigned int alignment = 0);
/// @}

#if defined(HAVE_BLITZ)
/// Write a blitz array (which can be non-contiguous).
template <typename T, unsigned int N>
BlobOStream& operator<<(BlobOStream&, const blitz::Array<T, N>&);

/// Read back a blitz array.
/// The dimensionality found in the stream has to match N.
/// If the shape mismatches, the array is resized.
/// If the shape matches, the array can be non-contiguous.
template <typename T, unsigned int N>
BlobIStream& operator>>(BlobIStream&, blitz::Array<T, N>&);
#endif

#if defined(HAVE_AIPSPP)
/// Write an AIPS++ array (which can be non-contiguous).
template <typename T>
BlobOStream& operator<<(BlobOStream&, const casacore::Array<T>&);

/// Read back an AIPS++ array.
/// If the shape mismatches, the array is resized.
/// If the shape matches, the array can be non-contiguous.
template <typename T>
BlobIStream& operator>>(BlobIStream&, casacore::Array<T>&);

/// Write/read the shape of an AIPS++ array.
/// @{
BlobOStream& operator<<(BlobOStream&, const casacore::IPosition&);
BlobIStream& operator>>(BlobIStream&, casacore::IPosition&);
/// @}
#endif

/// \name Write a vector of objects
/// @{
BlobOStream& operator<<(BlobOStream&, const std::vector<bool>&);
template <typename T>
BlobOStream& operator<<(BlobOStream&, const std::vector<T>&);
/// @}

/// \name Read back a vector of objects
/// The dimensionality found in the stream has to be 1.
/// The vector is resized as needed.
/// @{
BlobIStream& operator>>(BlobIStream&, std::vector<bool>&);
template <typename T>
BlobIStream& operator>>(BlobIStream&, std::vector<T>&);
/// @}

/// Read back as a C-array with the axes in Fortran or C-style order.
/// It allocates the required storage and puts the pointer to it in arr.
/// The shape is stored in the vector.
/// The user is responsible for deleting the data.
template <typename T>
BlobIStream& getBlobArray(BlobIStream& bs, T*& arr,
                          std::vector<uint64_t>& shape, bool fortranOrder);

/// Find an array in the blob with the axes in Fortran or C-style order.
/// The shape is stored in the vector.
/// It returns the offset of the array in the buffer and treats the array
/// as being read (thus skips over the data).
/// It is only possible if the underlying buffer is seekable.
/// It is meant for use with the BlobIBufString buffer. The function
/// getPointer in that class can be used to turn the position into a pointer.
/// The data are assumed to be aligned on the given alignment which
/// defaults to sizeof(T) bytes with a maximum of 8.
template <typename T>
uint64_t getSpaceBlobArray(BlobIStream& bs, bool useBlobHeader,
                           std::vector<uint64_t>& shape, bool fortranOrder);

/// Reserve space for a 1-dim array of the given size.
template <typename T>
inline uint64_t setSpaceBlobArray1(BlobOStream& bs, bool useBlobHeader,
                                   uint64_t size0, unsigned int alignment) {
  return setSpaceBlobArray<T>(bs, useBlobHeader, &size0, 1, true, alignment);
}

/// Put a vector object as an array.
template <typename T>
inline BlobOStream& operator<<(BlobOStream& bs, const std::vector<T>& vec) {
  return putBlobVector(bs, vec.empty() ? 0 : &(vec[0]), vec.size());
}

/// Put a C-style vector of values as an array.
template <typename T>
inline BlobOStream& putBlobVector(BlobOStream& bs, const T* vec,
                                  uint64_t size) {
  return putBlobArray(bs, vec, &size, 1, true);
}

/// Put a blob array header. It returns the number of elements in the array.
/// This is a helper function for the functions writing an array.
/// After writing the shape it aligns the stream on the given alignment.
uint64_t putBlobArrayHeader(BlobOStream& bs, bool useBlobHeader,
                            const std::string& headerName,
                            const uint64_t* shape, uint16_t ndim,
                            bool fortranOrder, unsigned int alignment);

/// Get the ordering and dimensionality.
/// This is a helper function for the functions reading an array.
/// It returns the number of alignment bytes used.
inline unsigned int getBlobArrayStart(BlobIStream& bs, bool& fortranOrder,
                                      uint16_t& ndim) {
  unsigned char nalign;
  bs >> fortranOrder >> nalign >> ndim;
  return nalign;
}

/// Convert the array header data.
void convertArrayHeader(common::DataFormat, char* header, bool useBlobHeader);

/// Get the shape of an array from the blob.
/// This is a helper function for the functions reading an array.
/// It returns the number of elements in the array.
uint64_t getBlobArrayShape(BlobIStream& bs, uint64_t* shape, unsigned int ndim,
                           bool swapAxes, unsigned int nalign);

/// Helper function to put an array of data.
/// It is specialized for the standard types (including complex and string).
template <typename T>
void putBlobArrayData(BlobOStream& bs, const T* data, uint64_t nr);

/// Helper function to get an array of data.
template <typename T>
void getBlobArrayData(BlobIStream& bs, T* data, uint64_t nr);

/// Specializations for the standard types (including complex and string).
#define BLOBARRAY_PUTGET_SPEC(TP)                                              \
  template <>                                                                  \
  inline void putBlobArrayData(BlobOStream& bs, const TP* data, uint64_t nr) { \
    bs.put(data, nr);                                                          \
  }                                                                            \
  template <>                                                                  \
  inline void getBlobArrayData(BlobIStream& bs, TP* data, uint64_t nr) {       \
    bs.get(data, nr);                                                          \
  }
BLOBARRAY_PUTGET_SPEC(bool)
BLOBARRAY_PUTGET_SPEC(int8_t)
BLOBARRAY_PUTGET_SPEC(uint8_t)
BLOBARRAY_PUTGET_SPEC(int16_t)
BLOBARRAY_PUTGET_SPEC(uint16_t)
BLOBARRAY_PUTGET_SPEC(int32_t)
BLOBARRAY_PUTGET_SPEC(uint32_t)
BLOBARRAY_PUTGET_SPEC(int64_t)
BLOBARRAY_PUTGET_SPEC(uint64_t)
BLOBARRAY_PUTGET_SPEC(float)
BLOBARRAY_PUTGET_SPEC(double)
BLOBARRAY_PUTGET_SPEC(std::complex<float>)
BLOBARRAY_PUTGET_SPEC(std::complex<double>)
BLOBARRAY_PUTGET_SPEC(std::string)

/// @}

}  // namespace blob
}  // namespace dp3

#include "BlobArray.tcc"

using dp3::blob::operator<<;
using dp3::blob::operator>>;
using dp3::blob::getBlobArray;
using dp3::blob::getSpaceBlobArray;
using dp3::blob::putBlobArray;
using dp3::blob::putBlobVector;
using dp3::blob::setSpaceBlobArray;
using dp3::blob::setSpaceBlobArray1;
using dp3::blob::setSpaceBlobArray2;
using dp3::blob::setSpaceBlobArray3;
using dp3::blob::setSpaceBlobArray4;

#endif
