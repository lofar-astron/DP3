// BlobIStream.h: Input stream for a blob
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Always #include <lofar_config.h> first!
#include "BlobIStream.h"
#include "BlobHeader.h"

#include "../common/DataConvert.h"

#include "BlobException.h"

#include <cassert>

namespace dp3 {
namespace blob {

BlobIStream::BlobIStream(BlobIBuffer& bb)
    : itsMustConvert(false),
      itsHasCachedType(false),
      itsCurLength(0),
      itsLevel(0),
      itsVersion(0),
      itsStream(&bb) {
  itsSeekable = (bb.tellPos() != -1);
}

BlobIStream::~BlobIStream() {
  // We could check that isLevel==0, but in case this destructor is called
  // by an exception, it gives another exception causing an abort.
  // So it should only be done if it is known that the destructor is called
  // in the normal way, but I don't know how to test that.
  // Maybe uncaught_exception is the way to go.
  /////  assert (itsLevel == 0);
}

// getNextType gets the object type of the next piece of
// information to read. It can only be used if a file has been
// opened and if no put is in operation.
// It checks if it finds the correct magic value preceeding
// the object type.
const std::string& BlobIStream::getNextType() {
  uint64_t size;
  return getNextType(size);
}

const std::string& BlobIStream::getNextType(uint64_t& size) {
  // Return current type if already cached.
  if (itsHasCachedType) {
    return itsObjectType;
  }
  // Read header and check the magic value and the level.
  BlobHeader hdr;
  itsStream->get((char*)(&hdr), sizeof(hdr));
  assert(hdr.checkMagicValue());
  // The level does not need to be equal in case of BlobFieldSet::putExtraBlob.
  assert(itsLevel >= hdr.itsLevel);
  // Determine if data has to be converted (in case data format mismatches).
  if (itsLevel == 0) {
    itsMustConvert = hdr.mustConvert();
    itsDataFormat = common::DataFormat(hdr.itsDataFormat);
  } else {
    assert(hdr.itsDataFormat == itsDataFormat);
  }
  // Keep the current length read.
  itsLevel++;
  itsObjLen.push(itsCurLength);
  uint64_t sz = hdr.getLength();
  itsObjTLN.push(sz);
  size = sz;
  itsVersion = hdr.getVersion();
  itsObjectType.resize(hdr.itsNameLength);  // resize string adding trailing 0
  char* ptr = &(itsObjectType[0]);
  itsCurLength = sizeof(hdr);  // length read
  if (hdr.itsNameLength > 0) {
    getBuf(ptr, hdr.itsNameLength);  // read objecttype
  }
  itsHasCachedType = true;
  return itsObjectType;
}

int BlobIStream::getStart(const std::string& type) {
  // Read the header and check if the type matches.
  if (type != getNextType()) {
    throw BlobException("BlobIStream::getStart: found object type " +
                        getNextType() + ", expected " + type);
  }
  itsHasCachedType = false;  // type is not cached anymore
  return itsVersion;
}

uint64_t BlobIStream::getEnd() {
  assert(itsLevel > 0);
  uint32_t eob;
  *this >> eob;
  if (eob != BlobHeader::eobMagicValue()) {
    throw BlobException("BlobIStream::getEnd - no end-of-blob value found");
  }
  uint64_t toRead = itsObjTLN.top();
  uint64_t len = itsCurLength;
  itsCurLength = itsObjLen.top();
  itsObjTLN.pop();
  itsObjLen.pop();
  if (itsLevel > 0) {
    if (!(len == toRead || toRead == 0)) {
      throw BlobException("BlobIStream::getEnd: part of object not read");
    }
  }
  if (--itsLevel > 0) {
    itsCurLength += len;
  }
  return len;
}

void BlobIStream::getBuf(void* buf, uint64_t sz) {
  checkGet();
  uint64_t sz1 = itsStream->get(static_cast<char*>(buf), sz);
  if (sz1 != sz) {
    throw BlobException("BlobIStream::getBuf - " + std::to_string(sz) +
                        " bytes asked, but only " + std::to_string(sz1) +
                        " could be read (pos=" + std::to_string(tellPos()) +
                        ")");
  }
  itsCurLength += sz1;
}

BlobIStream& BlobIStream::operator>>(bool& var) {
  char v;
  getBuf(&v, 1);
  var = v;
  return *this;
}
BlobIStream& BlobIStream::operator>>(char& var) {
  getBuf(&var, 1);
  return *this;
}
BlobIStream& BlobIStream::operator>>(int8_t& var) {
  getBuf(&var, 1);
  return *this;
}
BlobIStream& BlobIStream::operator>>(uint8_t& var) {
  getBuf(&var, 1);
  return *this;
}
BlobIStream& BlobIStream::operator>>(int16_t& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    var = common::dataConvert(itsDataFormat, var);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(uint16_t& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    var = common::dataConvert(itsDataFormat, var);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(int32_t& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    var = common::dataConvert(itsDataFormat, var);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(uint32_t& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    var = common::dataConvert(itsDataFormat, var);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(int64_t& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    var = common::dataConvert(itsDataFormat, var);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(uint64_t& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    var = common::dataConvert(itsDataFormat, var);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(float& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    common::dataConvertFloat(itsDataFormat, &var);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(double& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    common::dataConvertDouble(itsDataFormat, &var);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(std::complex<float>& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    common::dataConvertFloat(itsDataFormat, &var, 2);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(std::complex<double>& var) {
  getBuf(&var, sizeof(var));
  if (itsMustConvert) {
    common::dataConvertDouble(itsDataFormat, &var, 2);
  }
  return *this;
}
BlobIStream& BlobIStream::operator>>(std::string& var) {
  int64_t len;
  operator>>(len);
  var.resize(len);        // resize storage
  char* ptr = &(var[0]);  // get actual string
  getBuf(ptr, len);
  return *this;
}

void BlobIStream::get(bool* values, uint64_t nrval) {
  unsigned char buf[256];
  while (nrval > 0) {
    unsigned int nr = std::min(nrval, uint64_t(8 * 256));
    // Get and convert bits to bools.
    int nrb = (nr + 7) / 8;
    getBuf(buf, nrb);
    common::bitToBool(values, buf, nr);
    nrval -= nr;
    values += nr;
  }
}
void BlobIStream::get(char* values, uint64_t nrval) { getBuf(values, nrval); }
void BlobIStream::get(int8_t* values, uint64_t nrval) { getBuf(values, nrval); }
void BlobIStream::get(uint8_t* values, uint64_t nrval) {
  getBuf(values, nrval);
}
void BlobIStream::get(int16_t* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(int16_t));
  if (itsMustConvert) {
    common::dataConvert16(itsDataFormat, values, nrval);
  }
}
void BlobIStream::get(uint16_t* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(uint16_t));
  if (itsMustConvert) {
    common::dataConvert16(itsDataFormat, values, nrval);
  }
}
void BlobIStream::get(int32_t* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(int32_t));
  if (itsMustConvert) {
    common::dataConvert32(itsDataFormat, values, nrval);
  }
}
void BlobIStream::get(uint32_t* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(uint32_t));
  if (itsMustConvert) {
    common::dataConvert32(itsDataFormat, values, nrval);
  }
}
void BlobIStream::get(int64_t* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(int64_t));
  if (itsMustConvert) {
    common::dataConvert64(itsDataFormat, values, nrval);
  }
}
void BlobIStream::get(uint64_t* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(uint64_t));
  if (itsMustConvert) {
    common::dataConvert64(itsDataFormat, values, nrval);
  }
}
void BlobIStream::get(float* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(float));
  if (itsMustConvert) {
    common::dataConvertFloat(itsDataFormat, values, nrval);
  }
}
void BlobIStream::get(double* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(double));
  if (itsMustConvert) {
    common::dataConvertDouble(itsDataFormat, values, nrval);
  }
}
void BlobIStream::get(std::complex<float>* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(std::complex<float>));
  if (itsMustConvert) {
    common::dataConvertFloat(itsDataFormat, values, 2 * nrval);
  }
}
void BlobIStream::get(std::complex<double>* values, uint64_t nrval) {
  getBuf(values, nrval * sizeof(std::complex<double>));
  if (itsMustConvert) {
    common::dataConvertDouble(itsDataFormat, values, 2 * nrval);
  }
}
void BlobIStream::get(std::string* values, uint64_t nrval) {
  for (uint64_t i = 0; i < nrval; i++) {
    *this >> values[i];
  }
}

void BlobIStream::getBoolVec(std::vector<bool>& values, uint64_t sz) {
  values.resize(sz);
  bool buf[256];
  uint64_t inx = 0;
  while (sz > 0) {
    unsigned int nr = std::min(sz, uint64_t(256));
    // Get and convert bools to vector.
    get(buf, nr);
    for (unsigned int i = 0; i < nr; i++) {
      values[inx++] = buf[i];
    }
    sz -= nr;
  }
}

int64_t BlobIStream::getSpace(uint64_t nbytes) {
  checkGet();
  int64_t pos = tellPos();
  if (pos == -1) {
    throw BlobException(
        "BlobIStream::getSpace cannot be done; "
        "its BlobIBuffer is not seekable");
  }
  itsStream->setPos(pos + nbytes);
  itsCurLength += nbytes;
  return pos;
}

unsigned int BlobIStream::align(unsigned int n) {
  unsigned int nfill = 0;
  if (n > 1) {
    int64_t pos = tellPos();
    if (pos > 0) {
      nfill = pos % n;
    }
  }
  if (nfill > 0) {
    char fill;
    nfill = n - nfill;
    for (unsigned int i = 0; i < nfill; i++) {
      unsigned int sz1 = itsStream->get(&fill, 1);
      if (sz1 != 1) {
        throw BlobException("BlobIStream::align - could not read fill (pos=" +
                            std::to_string(tellPos()) + ")");
      }
      itsCurLength++;
    }
  }
  return nfill;
}

void BlobIStream::throwGet() const {
  throw BlobException("BlobIStream: getStart should be done first");
}

}  // namespace blob
}  // namespace dp3
