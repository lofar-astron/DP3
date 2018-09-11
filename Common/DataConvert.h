//# DataConvert.h: Global functions to convert data values
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
//# $Id: DataConvert.h 14057 2009-09-18 12:26:29Z diepen $

#ifndef LOFAR_COMMON_DATACONVERT_H
#define LOFAR_COMMON_DATACONVERT_H

// \file
// Global functions to convert data values

// Global functions to convert data values

#include "DataFormat.h"

#include <complex>

//# If std::complex is used for the complex types, their functions are
//# template specialisations, so they need template<>.
#ifndef LOFAR_BUILTIN_COMPLEXINT
# define LFDC_TMPL_INT template<>
# else
# define LFDC_TMPL_INT
#endif
#ifndef LOFAR_BUILTIN_COMPLEXFP
# define LFDC_TMPL_FP template<>
# else
# define LFDC_TMPL_FP
#endif

namespace DP3
{
// \ingroup Common
// \addtogroup DataConvert Data conversion functions
//
  // This file declares functions to convert data from one representation to
  // another, in particular from little endian to big endian (or vice-versa).
  //
  // The functions are defined in a general way for each standard data type,
  // so in principle every conceivable conversion could be done (for example,
  // from the old VAX format to IEEE format). However, currently byte swap
  // is the only conversion needed, so only that one is implemented.
  //
  // Furthermore it contains a function to convert bool values to bits
  // and vice-versa.
  // <group>
  
  // \name Convert the possible native types.
  // These functions can be used in templates.
  // <group>
  void dataConvert (DataFormat, char* inout, uint nrval);
  void dataConvert (DataFormat, int8_t* inout, uint nrval);
  void dataConvert (DataFormat, uint8_t* inout, uint nrval);
  void dataConvert (DataFormat, int16_t* inout, uint nrval);
  void dataConvert (DataFormat, uint16_t* inout, uint nrval);
  void dataConvert (DataFormat, int32_t* inout, uint nrval);
  void dataConvert (DataFormat, uint32_t* inout, uint nrval);
  void dataConvert (DataFormat, int64_t* inout, uint nrval);
  void dataConvert (DataFormat, uint64_t* inout, uint nrval);
  void dataConvert (DataFormat, float* inout, uint nrval);
  void dataConvert (DataFormat, double* inout, uint nrval);
  template<class T> void dataConvert (DataFormat, std::complex<T>* inout,
				      uint nrval);
  LFDC_TMPL_FP  void dataConvert (DataFormat, std::complex<float>* inout, uint nrval);
  LFDC_TMPL_FP  void dataConvert (DataFormat, std::complex<double>* inout, uint nrval);
  // </group>

  // \name Convert char, int8, or uint8.
  // Currently it simply returns the input.
  // <group>
  char dataConvert (DataFormat, char in);
  int8_t dataConvert (DataFormat, int8_t in);
  uint8_t dataConvert (DataFormat, uint8_t in);
  // </group>

  // \name Convert 16 bit integers.
  // <group>
  int16_t dataConvert (DataFormat, int16_t in);
  uint16_t dataConvert (DataFormat, uint16_t in);
  void dataConvert16 (DataFormat, void* out, const void* in);
  void dataConvert16 (DataFormat, void* inout);
  void dataConvert16 (DataFormat, void* out, const void* in, uint nrval);
  void dataConvert16 (DataFormat, void* inout, uint nrval);
  // </group>

  // \name Convert 32 bit integers.
  // <group>
  int32_t dataConvert (DataFormat, int32_t in);
  uint32_t dataConvert (DataFormat, uint32_t in);
  void dataConvert32 (DataFormat, void* out, const void* in);
  void dataConvert32 (DataFormat, void* inout);
  void dataConvert32 (DataFormat, void* out, const void* in, uint nrval);
  void dataConvert32 (DataFormat, void* inout, uint nrval);
  // </group>

  // \name Convert 64 bit integers.
  // <group>
  int64_t dataConvert (DataFormat, int64_t in);
  uint64_t dataConvert (DataFormat, uint64_t in);
  void dataConvert64 (DataFormat, void* out, const void* in);
  void dataConvert64 (DataFormat, void* inout);
  void dataConvert64 (DataFormat, void* out, const void* in, uint nrval);
  void dataConvert64 (DataFormat, void* inout, uint nrval);
  // </group>

  // \name Convert 32 bit floats.
  // <group>
  void dataConvertFloat (DataFormat, void* out, const void* in);
  void dataConvertFloat (DataFormat, void* inout);
  void dataConvertFloat (DataFormat, void* out, const void* in, uint nrval);
  void dataConvertFloat (DataFormat, void* inout, uint nrval);
  // </group>

  // \name Convert 64 bit floats.
  // <group>
  void dataConvertDouble (DataFormat, void* out, const void* in);
  void dataConvertDouble (DataFormat, void* inout);
  void dataConvertDouble (DataFormat, void* out, const void* in, uint nrval);
  void dataConvertDouble (DataFormat, void* inout, uint nrval);
  // </group>

  // \name Swap bytes in 16 bit values.
  // <group>
  int16_t byteSwap (int16_t in);
  uint16_t byteSwap (uint16_t in);
  void byteSwap16 (void* out, const void* in);
  void byteSwap16 (void* inout);
  void byteSwap16 (void* out, const void* in, uint nrval);
  void byteSwap16 (void* inout, uint nrval);
  // </group>

  // \name Swap bytes in 32 bit values.
  // <group>
  int32_t byteSwap (int32_t in);
  uint32_t byteSwap (uint32_t in);
  void byteSwap32 (void* out, const void* in);
  void byteSwap32 (void* inout);
  void byteSwap32 (void* out, const void* in, uint nrval);
  void byteSwap32 (void* inout, uint nrval);
  // </group>

  // \name Swap bytes in 64 bit values.
  // <group>
  int64_t byteSwap (int64_t in);
  uint64_t byteSwap (uint64_t in);
  void byteSwap64 (void* out, const void* in);
  void byteSwap64 (void* inout);
  void byteSwap64 (void* out, const void* in, uint nrval);
  void byteSwap64 (void* inout, uint nrval);
  // </group>

  // Convert bools to bits.
  // startbit gives to first bit to use in the to buffer.
  // It returns the number of bytes used.
  uint boolToBit (void* to, const void* from, uint nvalues, uint startbit=0);

  // Convert bits to bools.
  // startbit gives to first bit to use in the from buffer.
  // It returns the number of bytes used.
  uint bitToBool (void* to, const void* from, uint nvalues, uint startbit=0);

  // </group>

} // end namespace LOFAR


namespace DP3
{
  template<class T>
  inline void dataConvert (DataFormat fmt, std::complex<T>* inout, uint nrval)
    { dataConvert (fmt, (T*)inout, 2*nrval); }

  inline void dataConvert (DataFormat, char*, uint)
    {}
  inline void dataConvert (DataFormat, int8_t*, uint)
    {}
  inline void dataConvert (DataFormat, uint8_t*, uint)
    {}
  inline void dataConvert (DataFormat fmt, int16_t* inout, uint nrval)
    { dataConvert16 (fmt, inout, nrval); }
  inline void dataConvert (DataFormat fmt, uint16_t* inout, uint nrval)
    { dataConvert16 (fmt, inout, nrval); }
  inline void dataConvert (DataFormat fmt, int32_t* inout, uint nrval)
    { dataConvert32 (fmt, inout, nrval); }
  inline void dataConvert (DataFormat fmt, uint32_t* inout, uint nrval)
    { dataConvert32 (fmt, inout, nrval); }
  inline void dataConvert (DataFormat fmt, int64_t* inout, uint nrval)
    { dataConvert64 (fmt, inout, nrval); }
  inline void dataConvert (DataFormat fmt, uint64_t* inout, uint nrval)
    { dataConvert64 (fmt, inout, nrval); }
  inline void dataConvert (DataFormat fmt, float* inout, uint nrval)
    { dataConvert32 (fmt, inout, nrval); }
  inline void dataConvert (DataFormat fmt, double* inout, uint nrval)
    { dataConvert64 (fmt, inout, nrval); }
  inline void dataConvert (DataFormat fmt, std::complex<int16_t>* inout, uint nrval)
    { dataConvert16 (fmt, inout, 2*nrval); }
  LFDC_TMPL_FP inline void dataConvert (DataFormat fmt, std::complex<float>* inout, uint nrval)
    { dataConvertFloat (fmt, inout, 2*nrval); }
  LFDC_TMPL_FP inline void dataConvert (DataFormat fmt, std::complex<double>* inout, uint nrval)
    { dataConvertDouble (fmt, inout, 2*nrval); }

  inline char dataConvert (DataFormat, char in)
    { return in; }
  inline int8_t dataConvert (DataFormat, int8_t in)
    { return in; }
  inline uint8_t dataConvert (DataFormat, uint8_t in)
    { return in; }

  inline int16_t dataConvert (DataFormat, int16_t in)
    { return byteSwap (in); }
  inline uint16_t dataConvert (DataFormat, uint16_t in)
    { return byteSwap (in); }
  inline void dataConvert16 (DataFormat, void* out, const void* in)
    { byteSwap16 (out, in); }
  inline void dataConvert16 (DataFormat, void* inout)
    { byteSwap16 (inout); }
  inline void dataConvert16 (DataFormat, void* out, const void* in, uint nrval)
    { byteSwap16 (out, in, nrval); }
  inline void dataConvert16 (DataFormat, void* inout, uint nrval)
    { byteSwap16 (inout, nrval); }

  inline int32_t dataConvert (DataFormat, int32_t in)
    { return byteSwap (in); }
  inline uint32_t dataConvert (DataFormat, uint32_t in)
    { return byteSwap (in); }
  inline void dataConvert32 (DataFormat, void* out, const void* in)
    { byteSwap32 (out, in); }
  inline void dataConvert32 (DataFormat, void* inout)
    { byteSwap32 (inout); }
  inline void dataConvert32 (DataFormat, void* out, const void* in, uint nrval)
    { byteSwap32 (out, in, nrval); }
  inline void dataConvert32 (DataFormat, void* inout, uint nrval)
    { byteSwap32 (inout, nrval); }

  inline int64_t dataConvert (DataFormat, int64_t in)
    { return byteSwap (in); }
  inline uint64_t dataConvert (DataFormat, uint64_t in)
    { return byteSwap (in); }
  inline void dataConvert64 (DataFormat, void* out, const void* in)
    { byteSwap64 (out, in); }
  inline void dataConvert64 (DataFormat, void* inout)
    { byteSwap64 (inout); }
  inline void dataConvert64 (DataFormat, void* out, const void* in, uint nrval)
    { byteSwap64 (out, in, nrval); }
  inline void dataConvert64 (DataFormat, void* inout, uint nrval)
    { byteSwap64 (inout, nrval); }

  inline void dataConvertFloat (DataFormat, void* out, const void* in)
    { byteSwap32 (out, in); }
  inline void dataConvertFloat (DataFormat, void* inout)
    { byteSwap32 (inout); }
  inline void dataConvertFloat (DataFormat, void* out, const void* in, uint nrval)
    { byteSwap32 (out, in, nrval); }
  inline void dataConvertFloat (DataFormat, void* inout, uint nrval)
    { byteSwap32 (inout, nrval); }

  inline void dataConvertDouble (DataFormat, void* out, const void* in)
    { byteSwap64 (out, in); }
  inline void dataConvertDouble (DataFormat, void* inout)
    { byteSwap64 (inout); }
  inline void dataConvertDouble (DataFormat, void* out, const void* in, uint nrval)
    { byteSwap64 (out, in, nrval); }
  inline void dataConvertDouble (DataFormat, void* inout, uint nrval)
    { byteSwap64 (inout, nrval); }


  inline int16_t byteSwap (int16_t in)
  {
    int16_t v;
    byteSwap16 (&v, &in);
    return v;
  }

  inline uint16_t byteSwap (uint16_t in)
  {
    uint16_t v;
    byteSwap16 (&v, &in);
    return v;
  }

  inline void byteSwap16 (void* out, const void* in)
  {
    ((char*)(out))[0] = ((const char*)(in))[1];
    ((char*)(out))[1] = ((const char*)(in))[0];
  }

  inline void byteSwap16 (void* inout)
  {
    char v0 = ((const char*)(inout))[0];
    ((char*)(inout))[0] = ((const char*)(inout))[1];
    ((char*)(inout))[1] = v0;
  }


  inline int32_t byteSwap (int32_t in)
  {
    int32_t v;
    byteSwap32 (&v, &in);
    return v;
  }

  inline uint32_t byteSwap (uint32_t in)
  {
    uint32_t v;
    byteSwap32 (&v, &in);
    return v;
  }

  inline void byteSwap32 (void* out, const void* in)
  {
    ((char*)(out))[0] = ((const char*)(in))[3];
    ((char*)(out))[1] = ((const char*)(in))[2];
    ((char*)(out))[2] = ((const char*)(in))[1];
    ((char*)(out))[3] = ((const char*)(in))[0];
  }

  inline void byteSwap32 (void* inout)
  {
    char v0 = ((const char*)(inout))[0];
    char v1 = ((const char*)(inout))[1];
    ((char*)(inout))[0] = ((const char*)(inout))[3];
    ((char*)(inout))[1] = ((const char*)(inout))[2];
    ((char*)(inout))[2] = v1;
    ((char*)(inout))[3] = v0;
  }


  inline int64_t byteSwap (int64_t in)
  {
    int64_t v;
    byteSwap64 (&v, &in);
    return v;
  }

  inline uint64_t byteSwap (uint64_t in)
  {
    uint64_t v;
    byteSwap64 (&v, &in);
    return v;
  }

  inline void byteSwap64 (void* out, const void* in)
  {
    ((char*)(out))[0] = ((const char*)(in))[7];
    ((char*)(out))[1] = ((const char*)(in))[6];
    ((char*)(out))[2] = ((const char*)(in))[5];
    ((char*)(out))[3] = ((const char*)(in))[4];
    ((char*)(out))[4] = ((const char*)(in))[3];
    ((char*)(out))[5] = ((const char*)(in))[2];
    ((char*)(out))[6] = ((const char*)(in))[1];
    ((char*)(out))[7] = ((const char*)(in))[0];
  }

  inline void byteSwap64 (void* inout)
  {
    char v0 = ((const char*)(inout))[0];
    char v1 = ((const char*)(inout))[1];
    char v2 = ((const char*)(inout))[2];
    char v3 = ((const char*)(inout))[3];
    ((char*)(inout))[0] = ((const char*)(inout))[7];
    ((char*)(inout))[1] = ((const char*)(inout))[6];
    ((char*)(inout))[2] = ((const char*)(inout))[5];
    ((char*)(inout))[3] = ((const char*)(inout))[4];
    ((char*)(inout))[4] = v3;
    ((char*)(inout))[5] = v2;
    ((char*)(inout))[6] = v1;
    ((char*)(inout))[7] = v0;
  }

} // end namespace LOFAR


#undef LFDC_TMPL_FP
#undef LFDC_TMPL_INT

#endif
