// Box.cc:
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Box.h"

#include <iostream>

namespace dp3 {
namespace parmdb {

Box::Box(double x1, double x2, double y1, double y2, bool asStartEnd) {
  if (!asStartEnd) {
    x1 -= x2 * 0.5;
    x2 += x1;
    y1 -= y2 * 0.5;
    y2 += y1;
  }
  itsStart = Point(x1, y1);
  itsEnd = Point(x2, y2);
}

Box unite(const Box& lhs, const Box& rhs) {
  Point start(std::min(lhs.lowerX(), rhs.lowerX()),
              std::min(lhs.lowerY(), rhs.lowerY()));
  Point end(std::max(lhs.upperX(), rhs.upperX()),
            std::max(lhs.upperY(), rhs.upperY()));
  return Box(start, end);
}

Box intersect(const Box& lhs, const Box& rhs) {
  Point start(std::max(lhs.lowerX(), rhs.lowerX()),
              std::max(lhs.lowerY(), rhs.lowerY()));
  Point end(std::min(lhs.upperX(), rhs.upperX()),
            std::min(lhs.upperY(), rhs.upperY()));
  if (start.first < end.first && !casacore::near(start.first, end.first) &&
      start.second < end.second && !casacore::near(start.second, end.second)) {
    return Box(start, end);
  }
  return Box();
}

Box::Box(const std::vector<double>& values) {
  double stx = -1e30;
  double sty = -1e30;
  double enx = 1e30;
  double eny = 1e30;
  int sz = values.size();
  if (sz > 4) sz = 4;
  switch (sz) {
    case 4:
      eny = values[3];
      // fall through
    case 3:
      enx = values[2];
      // fall through
    case 2:
      sty = values[1];
      // fall through
    case 1:
      stx = values[0];
      break;
    default:
      break;
  }
  assert(stx <= enx);
  assert(sty <= eny);
  itsStart = Point(stx, sty);
  itsEnd = Point(enx, eny);
}

void Box::print() const {
  std::cout << itsStart.first << "\t" << itsStart.second << "\t" << itsEnd.first
            << "\t" << itsEnd.second << '\n';
}

}  // namespace parmdb
}  // namespace dp3
