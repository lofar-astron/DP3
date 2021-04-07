// Grid.cc:
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Grid.h"

#include <algorithm>

namespace dp3 {
namespace parmdb {

GridRep::GridRep()
    : itsAxes{std::make_shared<RegularAxis>(), std::make_shared<RegularAxis>()},
      itsHash(0),
      itsIsDefault(true) {
  init();
}

GridRep::GridRep(Axis::ShPtr first, Axis::ShPtr second)
    : itsAxes{first, second}, itsHash(0), itsIsDefault(false) {
  init();
}

GridRep::GridRep(const std::vector<Box>& domains, bool unsorted)
    : itsIsDefault(false) {
  if (domains.empty()) {
    itsIsDefault = true;
    itsAxes[0] = std::make_shared<RegularAxis>();
    itsAxes[1] = std::make_shared<RegularAxis>();
  } else {
    if (unsorted) {
      std::vector<Box> sortDomains(domains);
      std::sort(sortDomains.begin(), sortDomains.end());
      setup(sortDomains);
    } else {
      setup(domains);
    }
  }
  init();
}

GridRep::GridRep(const std::vector<Grid>& grids, bool unsorted)
    : itsIsDefault(false) {
  if (grids.empty()) {
    itsIsDefault = true;
    itsAxes[0] = std::make_shared<RegularAxis>();
    itsAxes[1] = std::make_shared<RegularAxis>();
  } else {
    if (unsorted) {
      std::vector<Grid> sortGrids(grids);
      std::sort(sortGrids.begin(), sortGrids.end());
      setup(sortGrids);
    } else {
      setup(grids);
    }
  }
  init();
}

void GridRep::setup(const std::vector<Box>& domains) {
  double sx = domains[0].lowerX();
  double ex = domains[0].upperX();
  double dx = ex - sx;
  double sy = domains[0].lowerY();
  double ey = domains[0].upperY();
  double dy = ey - sy;
  // Determine nr of cells in X and test if the cells are regular or not.
  // They are regular if equal width and consecutive.
  unsigned int nx = 1;
  bool xregular = true;
  while (nx < domains.size()) {
    if (sy != domains[nx].lowerY()) {
      break;
    }
    if (!(casacore::near(ex, domains[nx].lowerX()) &&
          casacore::near(dx, domains[nx].upperX() - domains[nx].lowerX()))) {
      xregular = false;
    }
    ex = domains[nx].upperX();
    ++nx;
  }
  // Assure there is an integer nr of y domains.
  assert(domains.size() % nx == 0);
  unsigned int ny = domains.size() / nx;
  // Check if the x axis is the same for all y-s.
  std::vector<double> xaxisStart;
  std::vector<double> xaxisEnd;
  xaxisStart.reserve(nx);
  xaxisEnd.reserve(nx);
  for (unsigned int i = 0; i < nx; ++i) {
    xaxisStart.push_back(domains[i].lowerX());
    xaxisEnd.push_back(domains[i].upperX());
  }
  for (unsigned int j = 1; j < ny; ++j) {
    for (unsigned int i = 0; i < nx; ++i) {
      assert(casacore::near(xaxisStart[i], domains[j * nx + i].lowerX()));
      assert(casacore::near(xaxisEnd[i], domains[j * nx + i].upperX()));
    }
  }
  // Determine the start/end for Y and if it is regular.
  // Check if the y axis is the same for all x-s.
  std::vector<double> yaxisStart;
  std::vector<double> yaxisEnd;
  yaxisStart.reserve(ny);
  yaxisEnd.reserve(ny);
  bool yregular = true;
  ey = sy;
  for (unsigned int i = 0; i < ny; ++i) {
    unsigned int inx = i * nx;
    yaxisStart.push_back(domains[inx].lowerY());
    yaxisEnd.push_back(domains[inx].upperY());
    if (!(casacore::near(ey, domains[inx].lowerY()) &&
          casacore::near(dy, domains[inx].upperY() - domains[inx].lowerY()))) {
      yregular = false;
    }
    ey = domains[inx].upperY();
  }
  for (unsigned int j = 0; j < ny; ++j) {
    for (unsigned int i = 1; i < nx; ++i) {
      assert(casacore::near(yaxisStart[j], domains[j * nx + i].lowerY()));
      assert(casacore::near(yaxisEnd[j], domains[j * nx + i].upperY()));
    }
  }
  // Create the (ir)regular axis.
  // Note that OrderedAxis checks if the intervals are ordered.
  if (xregular) {
    itsAxes[0] = std::make_shared<RegularAxis>(xaxisStart[0], dx, nx);
  } else {
    itsAxes[0] = std::make_shared<OrderedAxis>(xaxisStart, xaxisEnd, true);
  }
  if (yregular) {
    itsAxes[1] = std::make_shared<RegularAxis>(yaxisStart[0], dy, ny);
  } else {
    itsAxes[1] = std::make_shared<OrderedAxis>(yaxisStart, yaxisEnd, true);
  }
}

void GridRep::setup(const std::vector<Grid>& grids) {
  // First form a grid from the bounding boxes.
  // This checks if it is regular and divides in x and y.
  std::vector<Box> domains;
  domains.reserve(grids.size());
  for (std::vector<Grid>::const_iterator iter = grids.begin();
       iter != grids.end(); ++iter) {
    domains.push_back(iter->getBoundingBox());
  }
  setup(domains);
  // Now combine the grid axes.
  unsigned int nx = itsAxes[0]->size();
  unsigned int ny = itsAxes[1]->size();
  itsAxes[0] = combineAxes(grids, 0, nx, 1);
  itsAxes[1] = combineAxes(grids, 1, ny, nx);
  // Check if the grids themselves are equal for all x and y.
  // Check if the x axis is the same for all y-s.
  Axis::ShPtr xaxis = grids[0].getAxis(0);
  for (unsigned int j = 1; j < ny; ++j) {
    for (unsigned int i = 0; i < nx; ++i) {
      assert(*grids[i].getAxis(0) == *grids[j * nx + i].getAxis(0));
    }
  }
  // Check if the y axis is the same for all x-s.
  for (unsigned int j = 0; j < ny; ++j) {
    for (unsigned int i = 1; i < nx; ++i) {
      assert(*grids[j * nx].getAxis(1) == *grids[j * nx + i].getAxis(1));
    }
  }
}

Axis::ShPtr GridRep::combineAxes(const std::vector<Grid>& grids,
                                 unsigned int axnr, unsigned int n,
                                 unsigned int step) const {
  // Nothing to be done if only one cell.
  const Axis::ShPtr& faxis = grids[0].getAxis(axnr);
  if (n == 1) {
    return faxis;
  }
  // Count total number of cells and check if fully regular.
  bool isRegular = faxis->isRegular();
  double width = faxis->width(0);
  double last = faxis->upper(faxis->size() - 1);
  unsigned int ncells = faxis->size();
  for (unsigned int i = 1; i < n; ++i) {
    const Axis::ShPtr& axis = grids[i * step].getAxis(axnr);
    ncells += axis->size();
    if (isRegular) {
      isRegular = axis->isRegular() && (casacore::near(width, axis->width(0)) &&
                                        casacore::near(last, axis->lower(0)));
      last = axis->upper(axis->size() - 1);
    }
  }
  // If regular, return as such.
  if (isRegular) {
    return std::make_shared<RegularAxis>(faxis->lower(0), width, ncells);
  }
  // Alas irregular, so create an ordered axis.
  std::vector<double> starts;
  std::vector<double> ends;
  starts.reserve(ncells);
  ends.reserve(ncells);
  for (unsigned int i = 0; i < n; ++i) {
    const Axis::ShPtr& axis = grids[i * step].getAxis(axnr);
    for (unsigned int j = 0; j < axis->size(); ++j) {
      starts.push_back(axis->lower(j));
      ends.push_back(axis->upper(j));
    }
  }
  return std::make_shared<OrderedAxis>(starts, ends, true);
}

void GridRep::init() {
  // Calculate the hash value as a set of individual domains.
  // Thus add up the start and end values of all cells.
  const Axis& x = *itsAxes[0];
  const Axis& y = *itsAxes[1];
  int64_t xval = 0;
  int64_t yval = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    xval += int64_t(x.lower(i)) + int64_t(x.upper(i));
  }
  for (unsigned int i = 0; i < y.size(); ++i) {
    yval += int64_t(y.lower(i)) + int64_t(y.upper(i));
  }
  itsHash = x.size() * yval + y.size() * xval;
}

Grid::Grid(const std::vector<Grid>& grids, bool unsorted) {
  // If only one entry, we can simply make a copy.
  if (grids.size() == 1) {
    this->operator=(grids[0]);
  } else {
    itsRep = std::make_shared<GridRep>(grids, unsorted);
  }
}

bool Grid::operator==(const Grid& that) const {
  if (&(*itsRep) == &(*that.itsRep)) return true;
  if (hash() != that.hash()) return false;
  if (getAxis(0) != that.getAxis(0)) return false;
  if (getAxis(1) != that.getAxis(1)) return false;
  return true;
}

bool Grid::checkIntervals(const Grid& that) const {
  // Check per axis.
  return (getAxis(0)->checkIntervals(*that.getAxis(0)) &&
          getAxis(1)->checkIntervals(*that.getAxis(1)));
}

int64_t Grid::hash(const std::vector<Grid>& grids) {
  double val = 0;
  for (unsigned int i = 0; i < grids.size(); ++i) {
    val += grids[i].hash();
  }
  return val;
}

int64_t Grid::hash(const std::vector<Box>& domains) {
  double val = 0;
  for (unsigned int i = 0; i < domains.size(); ++i) {
    const Box& box = domains[i];
    val += int64_t(box.lowerX()) + int64_t(box.upperX()) +
           int64_t(box.lowerY()) + int64_t(box.upperY());
  }
  return val;
}

Grid Grid::subset(const Box& domain) const {
  Location index;
  return subset(domain, index);
}

Grid Grid::subset(const Box& domain, Location& index) const {
  return Grid(
      getAxis(0)->subset(domain.lowerX(), domain.upperX(), index.first),
      getAxis(1)->subset(domain.lowerY(), domain.upperY(), index.second));
}

Grid Grid::subset(const Location& start, const Location& end) const {
  assert(start.first <= end.first && start.second <= end.second);
  return Grid(getAxis(0)->subset(start.first, end.first),
              getAxis(1)->subset(start.second, end.second));
}

void Grid::toDomains(std::vector<Box>& domains) const {
  const Axis& xaxis = *(getAxis(0));
  const Axis& yaxis = *(getAxis(1));
  unsigned int nrx = nx();
  unsigned int nry = ny();
  // Prefetch the lower and upper values for both axes.
  std::vector<double> sx(nrx), ex(nrx);
  std::vector<double> sy(nrx), ey(nrx);
  for (unsigned int i = 0; i < nrx; ++i) {
    sx[i] = xaxis.lower(i);
    ex[i] = xaxis.upper(i);
  }
  for (unsigned int i = 0; i < nry; ++i) {
    sy[i] = yaxis.lower(i);
    ey[i] = yaxis.upper(i);
  }
  // Create the domains and append them to the vector.
  domains.reserve(domains.size() + size());
  for (unsigned int iy = 0; iy < nry; ++iy) {
    for (unsigned int ix = 0; ix < nrx; ++ix) {
      domains.push_back(Box(Point(sx[ix], sy[iy]), Point(ex[ix], ey[iy])));
    }
  }
}

}  // namespace parmdb
}  // namespace dp3
