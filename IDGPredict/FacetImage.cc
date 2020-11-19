// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "FacetImage.h"

namespace {

Facet clippedFacet(const Facet& input, int width, int height) {
  Facet clippedFacet(input);

  for (auto& v : clippedFacet) {
    if (v.x < 0) v.x = 0;
    if (v.y < 0) v.y = 0;
    if (v.x > width) v.x = width;
    if (v.y > height) v.y = height;
  }
  return clippedFacet;
}

}  // namespace

FacetImage::FacetImage() : _data(), _offsetX(0), _offsetY(0) {}

void FacetImage::Set(size_t width, size_t height, size_t spectralTerms,
                     double value) {
  _width = width;
  _height = height;
  _data.resize(spectralTerms);
  for (size_t i = 0; i != spectralTerms; ++i)
    _data[i].assign(width * height, value);
}

void FacetImage::CopyFacetPart(
    const Facet& facet, const std::vector<aocommon::UVector<double>>& inputs,
    size_t inputWidth, size_t inputHeight, double padding, bool makeSquare) {
  Facet cFacet = clippedFacet(facet, inputWidth, inputHeight);
  int x1, x2, y1, y2;
  cFacet.BoundingBox(x1, x2, y1, y2);
  size_t width = x2 - x1, height = y2 - y1;
  size_t paddedWidth = (size_t)ceil(width * padding);
  size_t paddedHeight = (size_t)ceil(height * padding);
  if (makeSquare) {
    paddedWidth = std::max(paddedWidth, paddedHeight);
    paddedHeight = paddedWidth;
  }
  // Make the width and height divisable by four.
  paddedWidth += (4 - (paddedWidth % 4)) % 4;
  paddedHeight += (4 - (paddedHeight % 4)) % 4;
  int padX = (paddedWidth - width) / 2;
  int padY = (paddedHeight - height) / 2;
  Set(paddedWidth, paddedHeight, inputs.size(), 0.0);
  _offsetX = x1 - padX;
  _offsetY = y1 - padY;
  for (size_t term = 0; term != inputs.size(); ++term) {
    for (int y = y1; y != y2; ++y) {
      int xi1, xi2;
      if (cFacet.HorizontalIntersections(y, xi1, xi2)) {
        if (xi1 < 0) xi1 = 0;
        if (xi2 < 0) xi2 = 0;
        if (xi2 - xi1 > int(width)) xi2 = xi1 + width;
        for (int x = xi1; x != xi2; ++x) {
          _data[term][(y - y1 + padY) * paddedWidth + x - x1 + padX] =
              inputs[term][y * inputWidth + x];
        }
      }
    }
  }
}

void FacetImage::FillFacet(const Facet& facet, int colour) {
  Facet cFacet = clippedFacet(facet, _width, _height);
  int x1, x2, y1, y2;
  cFacet.BoundingBox(x1, x2, y1, y2);
  size_t width = x2 - x1;

  for (size_t term = 0; term != _data.size(); ++term) {
    for (int y = y1; y != y2; ++y) {
      int xi1, xi2;
      if (cFacet.HorizontalIntersections(y, xi1, xi2)) {
        if (xi1 < 0) xi1 = 0;
        if (xi2 < 0) xi2 = 0;
        if (xi2 - xi1 > int(width)) xi2 = xi1 + width;
        for (int x = xi1; x != xi2; ++x) {
          _data[term][y * _width + x] = colour;
        }
      }
    }
  }
}
