// BlobException.h: Blob Exception class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_BLOB_BLOBEXCEPTION_H
#define LOFAR_BLOB_BLOBEXCEPTION_H

#include <stdexcept>

namespace DP3 {

/// \brief Blob Exception class.
typedef std::runtime_error BlobException;

}  // namespace DP3

#endif
