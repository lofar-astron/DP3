// DataConvert.cc: Global functions to convert data values
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Always #include <lofar_config.h> first!
#include "DataConvert.h"

void DP3::byteSwap16(void* val, unsigned int nrval) {
  char* v = (char*)val;
  for (unsigned int i = 0; i < nrval; i++) {
    DP3::byteSwap16(v);
    v += 2;
  }
}

void DP3::byteSwap16(void* out, const void* in, unsigned int nrval) {
  char* vout = (char*)out;
  const char* vin = (const char*)in;
  for (unsigned int i = 0; i < nrval; i++) {
    DP3::byteSwap16(vout, vin);
    vout += 2;
    vin += 2;
  }
}

void DP3::byteSwap32(void* val, unsigned int nrval) {
  char* v = (char*)val;
  for (unsigned int i = 0; i < nrval; i++) {
    DP3::byteSwap32(v);
    v += 4;
  }
}

void DP3::byteSwap32(void* out, const void* in, unsigned int nrval) {
  char* vout = (char*)out;
  const char* vin = (const char*)in;
  for (unsigned int i = 0; i < nrval; i++) {
    DP3::byteSwap32(vout, vin);
    vout += 4;
    vin += 4;
  }
}

void DP3::byteSwap64(void* val, unsigned int nrval) {
  char* v = (char*)val;
  for (unsigned int i = 0; i < nrval; i++) {
    DP3::byteSwap64(v);
    v += 8;
  }
}

void DP3::byteSwap64(void* out, const void* in, unsigned int nrval) {
  char* vout = (char*)out;
  const char* vin = (const char*)in;
  for (unsigned int i = 0; i < nrval; i++) {
    DP3::byteSwap64(vout, vin);
    vout += 8;
    vin += 8;
  }
}

unsigned int DP3::boolToBit(void* to, const void* from, unsigned int nvalues,
                            unsigned int startbit) {
  if (nvalues == 0) {
    return 0;
  }
  const bool* data = (const bool*)from;
  unsigned char* bits = (unsigned char*)to + startbit / 8;
  startbit %= 8;
  // Fill as many bytes as needed.
  unsigned int nbytes = (nvalues + startbit + 7) / 8;
  unsigned int i, j;
  unsigned int index = 0;
  {
    unsigned char mask = 1;
    mask <<= startbit;
    unsigned char& ch = bits[0];
    // Take care of correct number of bits in first byte.
    unsigned int nbits =
        (nvalues - index < 8 - startbit ? nvalues - index : 8 - startbit);
    for (j = 0; j < nbits; j++) {
      if (data[index++]) {
        ch |= mask;
      } else {
        ch &= ~mask;
      }
      mask <<= 1;
    }
  }
  for (i = 1; i < nbytes; ++i) {
    unsigned char mask = 1;
    unsigned char& ch = bits[i];
    ch = 0;
    // Take care of correct number of bits in last byte.
    unsigned int nbits = (nvalues - index < 8 ? nvalues - index : 8);
    for (j = 0; j < nbits; j++) {
      if (data[index++]) {
        ch |= mask;
      }
      mask <<= 1;
    }
  }
  return nbytes;
}

unsigned int DP3::bitToBool(void* to, const void* from, unsigned int nvalues,
                            unsigned int startbit) {
  bool* data = (bool*)to;
  const unsigned char* bits = (const unsigned char*)from + startbit / 8;
  startbit %= 8;
  // Fill as many bytes as needed.
  unsigned int nbytes = (nvalues + startbit + 7) / 8;
  unsigned int i, j;
  unsigned int index = 0;
  {
    unsigned char mask = 1;
    mask <<= startbit;
    const unsigned char ch = bits[0];
    // Take care of correct number of bits in first byte.
    unsigned int nbits =
        (nvalues - index < 8 - startbit ? nvalues - index : 8 - startbit);
    for (j = 0; j < nbits; j++) {
      data[index++] = ((ch & mask) != 0);
      mask <<= 1;
    }
  }
  for (i = 1; i < nbytes; ++i) {
    unsigned char mask = 1;
    const unsigned char ch = bits[i];
    // Take care of correct number of bits in last byte.
    unsigned int nbits = (nvalues - index < 8 ? nvalues - index : 8);
    for (j = 0; j < nbits; j++) {
      data[index++] = ((ch & mask) != 0);
      mask <<= 1;
    }
  }
  return nbytes;
}
