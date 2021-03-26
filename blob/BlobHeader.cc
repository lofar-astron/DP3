// BlobHeader.tcc: Standard header for a blob
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Always #include <lofar_config.h> first!

#include "BlobHeader.h"

#include "../common/DataFormat.h"

#include <cassert>

namespace dp3 {
namespace blob {

BlobHeader::BlobHeader(int version, unsigned int level)
    : itsLength(0),
      itsMagicValue(bobMagicValue()),
      itsVersion(version),
      itsDataFormat(common::dataFormat()),
      itsLevel(level),
      itsNameLength(0) {
  assert(version > -128 && version < 128);
  assert(level < 256);
}

void BlobHeader::setLocalDataFormat() {
  itsLength = dataConvert(getDataFormat(), itsLength);
  itsDataFormat = common::dataFormat();
}

}  // namespace blob
}  // namespace dp3
