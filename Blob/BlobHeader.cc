// BlobHeader.tcc: Standard header for a blob
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Always #include <lofar_config.h> first!

#include "BlobHeader.h"

#include "../Common/DataFormat.h"

#include <cassert>

DP3::BlobHeader::BlobHeader(int version, unsigned int level)
    : itsLength(0),
      itsMagicValue(bobMagicValue()),
      itsVersion(version),
      itsDataFormat(DP3::dataFormat()),
      itsLevel(level),
      itsNameLength(0) {
  assert(version > -128 && version < 128);
  assert(level < 256);
}

void DP3::BlobHeader::setLocalDataFormat() {
  itsLength = DP3::dataConvert(getDataFormat(), itsLength);
  itsDataFormat = DP3::dataFormat();
}
