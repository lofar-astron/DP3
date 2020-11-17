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

#ifndef FACET_H
#define FACET_H

#include <string>
#include <vector>
#include <cmath>

struct Vertex {
  Vertex(int _x, int _y) : x(_x), y(_y) {}
  int x, y;
};

class Facet {
 public:
  Facet() : _vertices(), _dirRA(0.0), _dirDec(0.0) {}

  typedef std::vector<Vertex>::iterator iterator;
  typedef std::vector<Vertex>::const_iterator const_iterator;

  iterator begin() { return _vertices.begin(); }
  iterator end() { return _vertices.end(); }
  const_iterator begin() const { return _vertices.begin(); }
  const_iterator end() const { return _vertices.end(); }

  void AddVertex(int x, int y) { _vertices.emplace_back(x, y); }

  void BoundingBox(int& x1, int& x2, int& y1, int& y2) const {
    if (_vertices.empty()) {
      x1 = 0;
      x2 = 0;
      y1 = 0;
      y2 = 0;
    } else {
      x1 = _vertices.front().x;
      x2 = _vertices.front().x;
      y1 = _vertices.front().y;
      y2 = _vertices.front().y;
      for (auto i = _vertices.begin() + 1; i != _vertices.end(); ++i) {
        x1 = std::min(x1, i->x);
        x2 = std::max(x2, i->x);
        y1 = std::min(y1, i->y);
        y2 = std::max(y2, i->y);
      }
    }
  }

  bool HorizontalIntersections(int yIntersect, int& x1, int& x2) const {
    size_t nInts = 0;
    x1 = 0;
    x2 = 0;
    for (size_t i = 0; i != _vertices.size(); ++i) {
      Vertex v1 = _vertices[i], v2 = _vertices[(i + 1) % _vertices.size()];
      if (v1.y > v2.y) std::swap(v1, v2);
      if (v1.y <= yIntersect && v2.y > yIntersect) {
        size_t x;
        if (v1.y == v2.y)
          x = std::min(v1.x, v2.x);
        else {
          double beta = double(v2.x - v1.x) / double(v2.y - v1.y);
          double xfl = v1.x + beta * (yIntersect - v1.y);
          x = round(xfl);
        }
        if (nInts == 0) {
          x1 = x;
          ++nInts;
        } else {
          x2 = x;
          if (x1 > x2) std::swap(x1, x2);
          return true;
        }
      }
    }
    return false;
  }

  double RA() const { return _dirRA; }
  double Dec() const { return _dirDec; }
  std::string Direction() const { return _direction; }

  void SetRA(double dirRA) { _dirRA = dirRA; }
  void SetDec(double dirDec) { _dirDec = dirDec; }
  void SetDirection(std::string direction) { _direction = direction; }

 private:
  std::vector<Vertex> _vertices;
  double _dirRA, _dirDec;
  std::string _direction;
};

#endif
