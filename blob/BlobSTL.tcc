// BlobSTL.tcc: Blob handling for STLs
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later


#ifndef LOFAR_BLOB_BLOBSTL_TCC
#define LOFAR_BLOB_BLOBSTL_TCC

#include "BlobSTL.h"
#include "BlobArray.h"

#include "../common/TypeNames.h"

namespace dp3 {
namespace blob {

template<typename T, typename U>
BlobOStream& operator<< (BlobOStream& bs, const std::map<T,U>& m)
{
  bs.putStart ("map<" + typeName((T*)0) + ',' + typeName((U*)0) + '>', 1);
  bs << static_cast<uint64_t>(m.size());
  for (typename std::map<T,U>::const_iterator it=m.begin();
        it!=m.end();
        ++it) {
    bs << it->first << it->second;
  }
  bs.putEnd();
  return bs;
}

template<typename T, typename U>
BlobIStream& operator>> (BlobIStream& bs, std::map<T,U>& m)
{
  bs.getStart ("map<" + typeName((T*)0) + ',' + typeName((U*)0) + '>');
  m.clear();
  uint64_t size;
  bs >> size;
  T t;
  U u;
  for (uint64_t i=0; i<size; ++i) {
    bs >> t >> u;
    m[t] = u;
  }
  bs.getEnd();
  return bs;
}

template<typename Seq>
void sequenceToBlob (BlobOStream& bs, const Seq& s)
{
  uint64_t n = s.size();
  putBlobArrayHeader (bs, true,
                      common::typeName((const typename Seq::value_type**)0),
                      &n, 1, true, 1);
  for (typename Seq::const_iterator it=s.begin();
        it!=s.end();
        ++it) {
    bs << *it;
  }
  bs.putEnd();
}

template<typename Seq>
void sequenceFromBlob (BlobIStream& bs, Seq& s)
{
  bs.getStart (common::typeName((const typename Seq::value_type**)0));
  bool fortranOrder;
  uint16_t ndim;
  unsigned int nalign = getBlobArrayStart (bs, fortranOrder, ndim);
  assert(ndim == 1);
  uint64_t size;
  getBlobArrayShape (bs, &size, 1, false, nalign);
  typename Seq::value_type t;
  s.clear();
  for (uint64_t i=0; i<size; ++i) {
    bs >> t;
    s.push_back (t);
  }
  bs.getEnd();
}

template<typename T>
void sequenceFromBlob (BlobIStream& bs, std::set<T>& s)
{
  bs.getStart (common::typeName((const T**)0));
  bool fortranOrder;
  uint16_t ndim;
  unsigned int nalign = getBlobArrayStart (bs, fortranOrder, ndim);
  assert(ndim == 1);
  uint64_t size;
  getBlobArrayShape (bs, &size, 1, false, nalign);
  T t;
  s.clear();
  for (uint64_t i=0; i<size; ++i) {
    bs >> t;
    s.insert (t);
  }
  bs.getEnd();
}

}
} // end namespace LOFAR

#endif
