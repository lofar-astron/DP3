// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef FACET_IMAGE_H
#define FACET_IMAGE_H

#include "Facet.h"

#include <aocommon/uvector.h>

#include <vector>

class FacetImage {
 public:
  FacetImage();

  size_t Width() const { return _width; }

  size_t Height() const { return _height; }

  int OffsetX() const { return _offsetX; }

  int OffsetY() const { return _offsetY; }

  void Set(size_t width, size_t height, size_t spectralTerms, double value);

  void CopyFacetPart(const Facet& facet,
                     const std::vector<aocommon::UVector<double>>& inputs,
                     size_t inputWidth, size_t inputHeight, double padding,
                     bool makeSquare);

  void FillFacet(const Facet& facet, int colour);

  double* Data(size_t spectralTerm) { return _data[spectralTerm].data(); }

  aocommon::UVector<double> AcquireData(size_t spectralTerm) {
    return std::move(_data[spectralTerm]);
  }

  /// A vector with Nterms elements, each holding the image data.
  /// Each uvector holds the data for one spectral frequency term.
  std::vector<aocommon::UVector<double>> _data;
  size_t _width, _height;
  int _offsetX, _offsetY;
};

#endif
