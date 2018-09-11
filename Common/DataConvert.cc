//# DataConvert.cc: Global functions to convert data values
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
//# $Id: DataConvert.cc 14057 2009-09-18 12:26:29Z diepen $

//# Always #include <lofar_config.h> first!
#include "DataConvert.h"

void DP3::byteSwap16 (void* val, uint nrval)
{
  char* v = (char*)val;
  for (uint i=0; i<nrval; i++) {
    DP3::byteSwap16 (v);
    v += 2;
  }
}

void DP3::byteSwap16 (void* out, const void* in, uint nrval)
{
  char* vout = (char*)out;
  const char* vin = (const char*)in;
  for (uint i=0; i<nrval; i++) {
    DP3::byteSwap16 (vout, vin);
    vout += 2;
    vin += 2;
  }
}

void DP3::byteSwap32 (void* val, uint nrval)
{
  char* v = (char*)val;
  for (uint i=0; i<nrval; i++) {
    DP3::byteSwap32 (v);
    v += 4;
  }
}

void DP3::byteSwap32 (void* out, const void* in, uint nrval)
{
  char* vout = (char*)out;
  const char* vin = (const char*)in;
  for (uint i=0; i<nrval; i++) {
    DP3::byteSwap32 (vout, vin);
    vout += 4;
    vin += 4;
  }
}

void DP3::byteSwap64 (void* val, uint nrval)
{
  char* v = (char*)val;
  for (uint i=0; i<nrval; i++) {
    DP3::byteSwap64 (v);
    v += 8;
  }
}

void DP3::byteSwap64 (void* out, const void* in, uint nrval)
{
  char* vout = (char*)out;
  const char* vin = (const char*)in;
  for (uint i=0; i<nrval; i++) {
    DP3::byteSwap64 (vout, vin);
    vout += 8;
    vin += 8;
  }
}

uint DP3::boolToBit (void* to, const void* from, uint nvalues, uint startbit)
{
  if (nvalues == 0) {
    return 0;
  }
  const bool* data = (const bool*)from;
  unsigned char* bits = (unsigned char*)to + startbit/8;
  startbit %= 8;
  //# Fill as many bytes as needed.
  uint nbytes = (nvalues + startbit + 7) / 8;
  uint i,j;
  uint index = 0;
  {
    unsigned char mask = 1;
    mask <<= startbit;
    unsigned char& ch = bits[0];
    //# Take care of correct number of bits in first byte.
    uint nbits = (nvalues-index < 8-startbit  ?  nvalues-index : 8-startbit);
    for (j=0; j<nbits; j++) {
      if (data[index++]) {
	ch |= mask;
      } else {
	ch &= ~mask;
      }
      mask <<= 1;
    }
  }
  for (i=1; i<nbytes; ++i) {
    unsigned char mask = 1;
    unsigned char& ch = bits[i];
    ch = 0;
    //# Take care of correct number of bits in last byte.
    uint nbits = (nvalues-index < 8  ?  nvalues-index : 8);
    for (j=0; j<nbits; j++) {
      if (data[index++]) {
	ch |= mask;
      }
      mask <<= 1;
    }
  }
  return nbytes;
}

uint DP3::bitToBool (void* to, const void* from, uint nvalues, uint startbit)
{
  bool* data = (bool*)to;
  const unsigned char* bits = (const unsigned char*)from + startbit/8;
  startbit %= 8;
  //# Fill as many bytes as needed.
  uint nbytes = (nvalues + startbit + 7) / 8;
  uint i,j;
  uint index = 0;
  {
    unsigned char mask = 1;
    mask <<= startbit;
    const unsigned char ch = bits[0];
    //# Take care of correct number of bits in first byte.
    uint nbits = (nvalues-index < 8-startbit  ?  nvalues-index : 8-startbit);
    for (j=0; j<nbits; j++) {
      data[index++] = ((ch & mask) != 0);
      mask <<= 1;
    }
  }
  for (i=1; i<nbytes; ++i) {
    unsigned char mask = 1;
    const unsigned char ch = bits[i];
    //# Take care of correct number of bits in last byte.
    uint nbits = (nvalues-index < 8  ?  nvalues-index : 8);
    for (j=0; j<nbits; j++) {
      data[index++] = ((ch & mask) != 0);
      mask <<= 1;
    }
  }
  return nbytes;
}
