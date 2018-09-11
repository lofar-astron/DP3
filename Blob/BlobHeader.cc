//# BlobHeader.tcc: Standard header for a blob
//#
//# Copyright (C) 2003
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the DP3 software suite.
//# The DP3 software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The DP3 software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the DP3 software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: BlobHeader.cc 14057 2009-09-18 12:26:29Z diepen $

//# Always #include <lofar_config.h> first!

#include "BlobHeader.h"

#include "../Common/DataFormat.h"

#include <cassert>

DP3::BlobHeader::BlobHeader (int version, uint level)
: itsLength         (0),
  itsMagicValue     (bobMagicValue()),
  itsVersion        (version),
  itsDataFormat     (DP3::dataFormat()),
  itsLevel          (level),
  itsNameLength     (0)
{
  assert (version > -128  &&  version < 128);
  assert (level < 256);
}
    
void DP3::BlobHeader::setLocalDataFormat()
{
  itsLength     = DP3::dataConvert (getDataFormat(), itsLength);
  itsDataFormat = DP3::dataFormat();
}
