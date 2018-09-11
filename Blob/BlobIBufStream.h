//# BlobIBufStream.h: Input buffer for a blob using an istream
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
//# $Id: BlobIBufStream.h 14057 2009-09-18 12:26:29Z diepen $

#ifndef LOFAR_BLOB_BLOBIBUFSTREAM_H
#define LOFAR_BLOB_BLOBIBUFSTREAM_H

// \file
// Input buffer for a blob using an istream

#include "BlobIBuffer.h"

#include <iosfwd>

namespace DP3 {

// \ingroup %pkgname%
  // @{
  
  // This class is the BlobIBuffer that makes use of an istream object.
  // The istream can be any type (ifstream, istringstream, ...).
  // It can, for instance, be used to read from a file or a socket.
  
  class BlobIBufStream : public BlobIBuffer
    {
    public:
      // Construct it with the underlying istream object.
      explicit BlobIBufStream (std::istream&);
      
      // Destructor.
      virtual ~BlobIBufStream();
      
      // Get the requested nr of bytes.
      virtual uint64_t get (void* buffer, uint64_t nbytes);
      
      // Get the position in the stream.
      // -1 is returned if the stream is not seekable.
      virtual int64_t tellPos() const;
      
      // Set the position in the stream.
      // It returns the new position which is -1 if the stream is not seekable.
      virtual int64_t setPos (int64_t pos);
      
    private:
      std::streambuf* itsStream;
    };

  // @}

} // end namespace

#endif
