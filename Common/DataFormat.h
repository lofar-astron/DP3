//# DataFormat.h: Get the data format (endian type)
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
//# $Id: DataFormat.h 14057 2009-09-18 12:26:29Z diepen $

#ifndef LOFAR_COMMON_DATAFORMAT_H
#define LOFAR_COMMON_DATAFORMAT_H

// \file
// Get the data format (endian type).
  // This file defines an enum for the possible machine data formats.
  // Currently only little and big endian is possible with floating point
  // numbers as IEEE and characters in the ASCII representation.
  // It is used in the Blob classes and the DataConvert functions.
  //
  // Furthermore it contains a function giving the data format in use on
  // the machine in use.

//# Never #include <config.h> or #include <lofar_config.h> in a header file!

namespace LOFAR
{
  
  enum DataFormat {LittleEndian=0, BigEndian=1};
  
  // Get the endian type on this machine.
  inline DataFormat dataFormat()
#if defined(WORDS_BIGENDIAN)
    {return BigEndian; }
#else
  {return LittleEndian; }
#endif
  
}


#endif
