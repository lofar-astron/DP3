// BlobException.h: Blob Exception class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBEXCEPTION_H
#define LOFAR_BLOB_BLOBEXCEPTION_H

#include <stdexcept>

namespace dp3 {
namespace blob {

/// \brief Blob Exception class.
typedef std::runtime_error BlobException;

}  // namespace blob
}  // namespace dp3

#endif
