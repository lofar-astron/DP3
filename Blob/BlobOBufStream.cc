//# BlobOBufStream.cc: Output buffer for a blob using an ostream
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
//# $Id: BlobOBufStream.cc 14057 2009-09-18 12:26:29Z diepen $

#include "BlobOBufStream.h"

#include <iostream>

namespace DP3 {

BlobOBufStream::BlobOBufStream (std::ostream& os)
: itsStream (os.rdbuf())
{}

BlobOBufStream::~BlobOBufStream()
{}

uint64_t BlobOBufStream::put (const void* buffer, uint64_t nbytes)
{
  return itsStream->sputn ((const char*)buffer, nbytes);
}

int64_t BlobOBufStream::tellPos() const
{
  return itsStream->pubseekoff (0, std::ios::cur);
}

int64_t BlobOBufStream::setPos (int64_t pos)
{
  return itsStream->pubseekoff (pos, std::ios::beg);
}

} // end namespace
