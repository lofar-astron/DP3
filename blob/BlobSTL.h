// BlobSTL.h: Blob handling for STL sequences
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBSTL_H
#define LOFAR_BLOB_BLOBSTL_H

#include "BlobOStream.h"
#include "BlobIStream.h"

#include <map>
#include <list>
#include <set>
#include <queue>
#include <deque>

namespace dp3 {
namespace blob {

/// \ingroup Blob
/// \brief Blob handling for STL sequences

/// @{

/// Define functions to write a map into a blob and to read it back.
/// The map is preceeded by the header 'map<T,U>', where T and U are
/// the type names as defined in TypeNames.h.
/// Type names are only defined for the basic types. Other types
/// are set to 'unknown'.
///   @{
/// \name Write a map.
template <typename T, typename U>
BlobOStream& operator<<(BlobOStream&, const std::map<T, U>&);

/// \name Read back a map.
template <typename T, typename U>
BlobIStream& operator>>(BlobIStream&, std::map<T, U>&);
///   @}

/// Define helper functions to write any STL sequence into a blob and to
/// read it back.
/// The sequence is preceeded by the header 'array<T>', where T is
/// the type name of Seq::value_type as defined in TypeNames.h.
/// Type names are only defined for the basic types. Other types
/// are set to 'unknown'.
/// All sequences are written in the same way (as 1-dim arrays). It means
/// that they can be read back using any other sequence type.
///   @{
/// \name Write a sequence.
template <typename Seq>
void sequenceToBlob(BlobOStream&, const Seq&);

/// \name Read back a sequence.
template <typename Seq>
void sequenceFromBlob(BlobOStream&, Seq&);
/// Specialize for a set.
template <typename T>
void sequenceFromBlob(BlobOStream&, std::set<T>&);
///   @}

/// Define helper functions to write an STL sequence into a blob and to
/// read it back.
/// The sequence is preceeded by the header 'array<T>', where T is
/// the type name of T as defined in TypeNames.h.
/// Type names are only defined for the basic types. Other types
/// are set to 'unknown'.
/// All sequences are written in the same way (as 1-dim arrays). It means
/// that they can be read back in any other sequence type (including
/// AIPS++ and Blitz arrays).
/// \note STL vectors are handled by BlobArray.h.
///   @{
/// \name Write a list.
template <typename T>
BlobOStream& operator<<(BlobOStream& bs, const std::list<T>& seq) {
  sequenceToBlob(bs, seq);
  return bs;
}

/// \name Read back a list.
template <typename T>
BlobIStream& operator>>(BlobIStream& bs, std::list<T>& seq) {
  sequenceFromBlob(bs, seq);
  return bs;
}

/// \name Write a set.
template <typename T>
BlobOStream& operator<<(BlobOStream& bs, const std::set<T>& seq) {
  sequenceToBlob(bs, seq);
  return bs;
}

/// \name Read back a set.
template <typename T>
BlobIStream& operator>>(BlobIStream& bs, std::set<T>& seq) {
  sequenceFromBlob(bs, seq);
  return bs;
}

/// \name Write a queue.
template <typename T>
BlobOStream& operator<<(BlobOStream& bs, const std::queue<T>& seq) {
  sequenceToBlob(bs, seq);
  return bs;
}

/// \name Read back a queue.
template <typename T>
BlobIStream& operator>>(BlobIStream& bs, std::queue<T>& seq) {
  sequenceFromBlob(bs, seq);
  return bs;
}

/// \name Write a deque.
template <typename T>
BlobOStream& operator<<(BlobOStream& bs, const std::deque<T>& seq) {
  sequenceToBlob(bs, seq);
  return bs;
}

/// \name Read back a deque.
template <typename T>
BlobIStream& operator>>(BlobIStream& bs, std::deque<T>& seq) {
  sequenceFromBlob(bs, seq);
  return bs;
}
///   @}

/// @}

}  // namespace blob
}  // namespace dp3

#include "BlobSTL.tcc"

using dp3::blob::operator<<;
using dp3::blob::operator>>;

#endif
