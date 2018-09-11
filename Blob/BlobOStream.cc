//# BlobOStream.cc: Output stream for a blob
//#
//# Copyright (C) 2003
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: BlobOStream.cc 14057 2009-09-18 12:26:29Z diepen $

//# Always #include <lofar_config.h> first!
#include "BlobOStream.h"
#include "BlobHeader.h"

#include "../Common/DataConvert.h"

#include "BlobException.h"

#include <cassert>
#include <iostream>

namespace DP3 {

BlobOStream::BlobOStream (BlobOBuffer& bb)
: itsCurLength (0),
  itsLevel     (0),
  itsStream    (&bb)
{
  itsSeekable = (bb.tellPos() != -1);
}

BlobOStream::~BlobOStream()
{
  // We could check that isLevel==0, but in case this destructor is called
  // by an exception, it gives another exception causing an abort.
  // So it should only be done if it is known that the destructor is called
  // in the normal way, but I don't know how to test that.
  // Maybe uncaught_exception is the way to go.
  /////  ASSERT (itsLevel == 0);
}

// putStart starts writing an object.
// It puts the object type and version and reserves space for the length.
// It increases the level for each object to hold the length.
uint BlobOStream::doPutStart (const char* type, uint nrc, int version)
{
  BlobHeader hdr(version, itsLevel);
  assert (nrc < 256);
  hdr.itsNameLength = nrc;
  itsObjLen.push (itsCurLength);         // length of outer blob
  itsObjPtr.push (tellPos()+hdr.lengthOffset()); // remember where to put len
  itsCurLength = 0;                      // initialize object length
  // Need to increment here, otherwise putBuf gives an exception.
  itsLevel++;
  // Put header plus objecttype.
  putBuf (&hdr, sizeof(hdr));
  if (nrc > 0) {
    putBuf (type, nrc);
  }
  return itsLevel;
}

// putend ends putting an object. It decreases the level and writes
// the object length if the file is seekable.
uint64_t BlobOStream::putEnd()
{
  assert (itsLevel > 0);
  uint32_t eob = BlobHeader::eobMagicValue();
  *this << eob;                        // write end-of-blob
  uint64_t len = itsCurLength;           // length of this object
  itsCurLength = itsObjLen.top();      // length of parent object
  int64_t pos    = itsObjPtr.top();
  itsObjLen.pop();
  itsObjPtr.pop();
  if (itsSeekable) {
    int64_t curpos = tellPos();
    itsStream->setPos (pos);
    itsStream->put (&len, sizeof(len));
    itsStream->setPos (curpos);
  }
  itsLevel--;
  if (itsLevel > 0) {
    itsCurLength += len;               // add length to parent object
  }
  return len;
}

void BlobOStream::putBuf (const void* buf, uint64_t sz)
{
  checkPut();
  uint64_t sz1 = itsStream->put (static_cast<const char*>(buf), sz);
  if (sz1 != sz) {
    throw BlobException(
	   "BlobOStream::putBuf - " + std::to_string(sz) + " bytes asked, but only "
	   + std::to_string(sz1) + " could be written (pos=" + std::to_string(tellPos()) + ")");
  }
  itsCurLength += sz1;
}

BlobOStream& BlobOStream::operator<< (const bool& var)
{
  char v = (var ? 1 : 0);
  putBuf (&v, 1);
  return *this;
}
BlobOStream& BlobOStream::operator<< (const char& var)
{
  putBuf (&var, 1);
  return *this;
}
BlobOStream& BlobOStream::operator<< (const int8_t& var)
{
  putBuf (&var, 1);
  return *this;
}
BlobOStream& BlobOStream::operator<< (const uint8_t& var)
{
  putBuf (&var, 1);
  return *this;
}
BlobOStream& BlobOStream::operator<< (const int16_t& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const uint16_t& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const int32_t& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const uint32_t& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const int64_t& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const uint64_t& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const float& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const double& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const fcomplex& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const dcomplex& var)
{
  putBuf (&var, sizeof(var));
  return *this;
}
BlobOStream& BlobOStream::operator<< (const std::string& var)
{
  operator<< (int64_t(var.size()));
  putBuf (var.data(), var.size());
  return *this;
}
BlobOStream& BlobOStream::operator<< (const char* var)
{
  int64_t sz = strlen(var);
  operator<< (sz);
  putBuf (var, sz);
  return *this;
}

void BlobOStream::put (const bool* values, uint64_t nrval)
{
  unsigned char buf[256];
  while (nrval > 0) {
    uint nr = std::min(nrval, uint64_t(8*256));
    // Convert to bits and put.
    uint nrb = DP3::boolToBit (buf, values, nr);
    putBuf (buf, nrb);
    nrval -= nr;
    values += nr;
  }
}
void BlobOStream::put (const char* values, uint64_t nrval)
{
  putBuf (values, nrval);
}
void BlobOStream::put (const int8_t* values, uint64_t nrval)
{
  putBuf (values, nrval);
}
void BlobOStream::put (const uint8_t* values, uint64_t nrval)
{
  putBuf (values, nrval);
}
void BlobOStream::put (const int16_t* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(int16_t));
}
void BlobOStream::put (const uint16_t* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(uint16_t));
}
void BlobOStream::put (const int32_t* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(int32_t));
}
void BlobOStream::put (const uint32_t* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(uint32_t));
}
void BlobOStream::put (const int64_t* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(int64_t));
}
void BlobOStream::put (const uint64_t* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(uint64_t));
}
void BlobOStream::put (const float* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(float));
}
void BlobOStream::put (const double* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(double));
}
void BlobOStream::put (const fcomplex* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(fcomplex));
}
void BlobOStream::put (const dcomplex* values, uint64_t nrval)
{
  putBuf (values, nrval*sizeof(dcomplex));
}
void BlobOStream::put (const std::string* values, uint64_t nrval)
{
  for (uint64_t i=0; i<nrval; i++) {
    *this << values[i];
  }
}

void BlobOStream::putBoolVec (const std::vector<bool>& values)
{
  uint64_t sz = values.size();
  // Convert to bools and put as such.
  bool buf[256];
  uint64_t inx=0;
  while (sz > 0) {
    uint nr = std::min(sz, uint64_t(256));
    for (uint i=0; i<nr; i++) {
      buf[i] = values[inx++];
    }
    put (buf, nr);
    sz -= nr;
  }
}

int64_t BlobOStream::setSpace (uint64_t nbytes)
{
  checkPut();
  int64_t pos = tellPos();
  if (pos == -1) {
    throw BlobException(
	   "BlobOStream::setSpace cannot be done; "
	   "its BlobOBuffer is not seekable");
  }
  itsStream->setPos (pos+nbytes);
  itsCurLength += nbytes;
  return pos;
}

uint BlobOStream::align (uint n)
{
  uint nfill = 0;
  if (n > 1) {
    int64_t pos = tellPos();
    if (pos > 0) {
      nfill = pos % n;
    }
  }
  if (nfill > 0) {
    char fill=0;
    nfill = n-nfill;
    for (uint i=0; i<nfill; i++) {
      uint sz1 = itsStream->put (&fill, 1);
      if (sz1 != 1) {
	throw BlobException(
	       "BlobOStream::align - could not write fill (pos="
	       + std::to_string(tellPos()) + ")");
      }
      itsCurLength++;
    }
  }
  return nfill;
}

void BlobOStream::throwPut() const
{
  throw BlobException("BlobOStream: putStart should be done first");
}


} // end namespace
