//# BlobIBuffer.h: Abstract base class for input buffer for a blob
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
//# $Id: BlobIBuffer.h 14057 2009-09-18 12:26:29Z diepen $

#include <cstdint>

#ifndef LOFAR_BLOB_BLOBIBUFFER_H
#define LOFAR_BLOB_BLOBIBUFFER_H

// \file
// Abstract base class for input buffer for a blob

namespace DP3 {

// \ingroup %pkgname%
  // @{
  
  // BlobIBuffer is the abstract base class for the source of a
  // BlobIStream object. In this way the source of a BlobIStream
  // can be chosen as needed. Currently three source types are possible:
  // <ul>
  // <li> BlobIBufStream uses an istream as the source. It is mainly
  //      meant to read blobs from a file, socket, etc.
  // <li> BlobIBufChar uses a plain C-array as the source. Its mainly used
  //      as a helper class, but can also be used in itself.
  // <li> BlobIBufVector<T> uses a vector<T> (where T is char or uchar)
  //      as the source.
  //      It uses BlobIBufChar to do the actual work.
  // <li> BlobIBufString uses a BlobString as the source.
  //      It uses BlobIBufChar to do the actual work.
  //      It is mainly meant for reading blobs stored using the LCS
  //      Persistency Layer.
  // </ul>
  
  class BlobIBuffer
    {
    public:
      // Constructor.
      BlobIBuffer()
	{};
      
      // Destructor.
      virtual ~BlobIBuffer()
	{};
      
      // Get the requested nr of bytes.
      virtual uint64_t get (void* buffer, uint64_t nbytes) = 0;
      
      // Get the position in the stream.
      // -1 is returned if the stream is not seekable.
      virtual int64_t tellPos() const = 0;
      
      // Set the position in the stream.
      // It returns the new position which is -1 if the stream is not seekable.
      virtual int64_t setPos (int64_t pos) = 0;
    };

  // @}

} // end namespace

#endif
