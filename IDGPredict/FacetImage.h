// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

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
