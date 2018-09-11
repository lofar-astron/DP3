//# BlobIBufStream.cc: Input buffer for a blob using an istream
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
//# $Id: BlobIBufStream.cc 14057 2009-09-18 12:26:29Z diepen $

#include "BlobIBufStream.h"

#include <iostream>

namespace DP3 {

BlobIBufStream::BlobIBufStream (std::istream& is)
: itsStream (is.rdbuf())
{}

BlobIBufStream::~BlobIBufStream()
{}

uint64_t BlobIBufStream::get (void* buffer, uint64_t nbytes)
{
  return itsStream->sgetn ((char*)buffer, nbytes);
}

int64_t BlobIBufStream::tellPos() const
{
  return itsStream->pubseekoff (0, std::ios::cur);
}

int64_t BlobIBufStream::setPos (int64_t pos)
{
  return itsStream->pubseekoff (pos, std::ios::beg);
}

} // end namespace
