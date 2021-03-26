// ParmValue.cc: A class containing the values of a parameter
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParmValue.h"

#include <casacore/casa/Arrays/Matrix.h>

using namespace casacore;
using namespace std;

namespace dp3 {
namespace parmdb {

ParmValue::ParmValue(double value) : itsErrors(0), itsRowId(-1) {
  setScalar(value);
}

ParmValue::ParmValue(const ParmValue& that) : itsErrors(0) { copyOther(that); }

ParmValue& ParmValue::operator=(const ParmValue& that) {
  if (this != &that) {
    copyOther(that);
  }
  return *this;
}

ParmValue::~ParmValue() { delete itsErrors; }

void ParmValue::copyOther(const ParmValue& that) {
  itsGrid = that.itsGrid;
  itsRowId = that.itsRowId;
  itsValues.assign(that.itsValues);  // ensure a copy is made
  delete itsErrors;
  itsErrors = 0;
  if (that.itsErrors) {
    itsErrors = new Array<double>;
    *itsErrors = *that.itsErrors;
  }
}

void ParmValue::setScalar(double value) {
  itsValues.resize(IPosition(2, 1, 1));
  itsValues = value;
}

void ParmValue::setCoeff(const casacore::Array<double>& values) {
  itsValues.assign(values);
}

void ParmValue::setScalars(const Grid& grid,
                           const casacore::Array<double>& values) {
  assert(int(grid.nx()) == values.shape()[0] &&
         int(grid.ny()) == values.shape()[1]);
  itsValues.assign(values);
  itsGrid = grid;
}

void ParmValue::setErrors(const casacore::Array<double>& errors) {
  // Check that the errors have the same shape as the values.
  assert(errors.shape().isEqual(itsValues.shape()));
  // Make sure a copy is made of the errors.
  if (!itsErrors) {
    itsErrors = new Array<double>();
  }
  itsErrors->assign(errors);
}

bool ParmValue::rescale(double sx, double ex, double sy, double ey,
                        const Box& oldDomain) {
  // No need to rescale if polynomial did not have a domain
  // or if the axes have length 1 and others match.
  casacore::Matrix<double> coeff(getValues());
  if (oldDomain.empty() || coeff.size() == 1 ||
      (coeff.nrow() == 1 && sy == oldDomain.lowerY() &&
       ey == oldDomain.upperY()) ||
      (coeff.ncolumn() == 1 && sx == oldDomain.lowerX() &&
       ex == oldDomain.upperX())) {
    return false;
  }
  // Rescale from coeff like:  s2a=s2/s1 and o2a=(o2-o1)/s1
  // x1=(x-o1)/s1 and x2=(x-o2)/s2; hence x2=(x1-(o2-o1)/s1)/(s2/s1)
  // where o is start of box and s is width of box (1 is old, 2 is new).
  // Note this must be the same as in AxisMapping that also uses
  // lower and width when scaling a polynomial.
  double s1x = oldDomain.widthX();
  double s1y = oldDomain.widthY();
  itsValues =
      scale2(coeff, (sx - oldDomain.lowerX()) / s1x,
             (sy - oldDomain.lowerY()) / s1y, (ex - sx) / s1x, (ey - sy) / s1y);
  return true;
}

// Scale a 2D polynomial using a given offset and scale factor.
// Polynomial is:
//  f(x,y) = sigma (a[i,j] * x**i * y**j)
// When scaling x to (x-ox)/sx  and  y to (y-oy)/sy each term gets
//  a[i,j] * (x*sx+ox)**i * (y*sy+oy)**j
//  a[i,j] * sigma((i p) * x**p * sx**p * ox**(i-p)) *
//           sigma((j q) * x**q * sy**q * oy**(j-q))
// So for each i,j,p,q each new coeff(p,q) gets terms
//  a[i,j] * (i p) * sx**p * ox**(i-p)) * (j q) * sy**q * oy**(j-q))
// The factors sx**p and sy**q are independent of i,j and are applied later.
casacore::Matrix<double> ParmValue::scale2(
    const casacore::Matrix<double>& coeff, double offx, double offy,
    double scalex, double scaley) {
  // Fill the Pascal triangle (till order 10) if not done yet.
  static casacore::Matrix<double> pascal;
  if (pascal.empty()) {
    fillPascal(pascal, 10);
  }
  // Rescale by looping over all possible terms (see above).
  int nx = coeff.shape()[0];
  int ny = coeff.shape()[1];
  assert(nx < int(pascal.nrow()) && ny < int(pascal.nrow()));
  casacore::Matrix<double> scoeff(coeff.shape(), 0.);
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      double offpy = coeff(ix, iy);
      for (int jy = iy; jy >= 0; --jy) {
        double offpx = pascal(jy, iy) * offpy;
        for (int jx = ix; jx >= 0; --jx) {
          scoeff(jx, jy) += pascal(jx, ix) * offpx;
          offpx *= offx;
        }
        offpy *= offy;
      }
    }
  }
  double scy = 1;
  for (int iy = 0; iy < ny; ++iy) {
    double scx = scy;
    for (int ix = 0; ix < nx; ++ix) {
      scoeff(ix, iy) *= scx;
      scx *= scalex;
    }
    scy *= scaley;
  }
  return scoeff;
}

void ParmValue::fillPascal(casacore::Matrix<double>& pascal, int order) {
  // pascal(j,i) gives (i over j).
  pascal.resize(order, order);
  pascal = 0.;
  for (int i = 0; i < order; ++i) {
    pascal(0, i) = 1.;
    for (int j = 1; j <= i; ++j) {
      pascal(j, i) = pascal(j - 1, i - 1) + pascal(j, i - 1);
    }
  }
}

ParmValueSet::ParmValueSet(const ParmValue& defaultValue,
                           ParmValue::FunkletType type, double perturbation,
                           bool pertRel, const Box& scaleDomain)
    : itsType(type),
      itsPerturbation(perturbation),
      itsPertRel(pertRel),
      itsDefaultValue(defaultValue),
      itsScaleDomain(scaleDomain),
      itsDirty(false) {
  if (type == ParmValue::Scalar) {
    if (defaultValue.getValues().size() != 1)
      throw std::runtime_error(
          "Default value of funklet type SCALAR can have one value only");
  }
}

ParmValueSet::ParmValueSet(const Grid& domainGrid,
                           const std::vector<ParmValue::ShPtr>& values,
                           const ParmValue& defaultValue,
                           ParmValue::FunkletType type, double perturbation,
                           bool pertRel)
    : itsType(type),
      itsPerturbation(perturbation),
      itsPertRel(pertRel),
      itsDomainGrid(domainGrid),
      itsValues(values),
      itsDefaultValue(defaultValue),
      itsDirty(false) {
  assert(domainGrid.size() == values.size() && values.size() > 0);
  if (type == ParmValue::Scalar) {
    if (defaultValue.getValues().size() != 1)
      throw std::runtime_error(
          "Default value of funklet type SCALAR can have one value only");
    for (unsigned int i = 0; i < values.size(); ++i) {
      if (values[i]->getValues().size() != values[i]->getGrid().size())
        throw std::runtime_error(
            "ParmValues of funklet type SCALAR must contain scalar values");
    }
  }
}

ParmValueSet::ParmValueSet(const ParmValueSet& that) { operator=(that); }

ParmValueSet& ParmValueSet::operator=(const ParmValueSet& that) {
  if (this != &that) {
    itsType = that.itsType;
    itsPerturbation = that.itsPerturbation;
    itsPertRel = that.itsPertRel;
    itsSolvableMask.assign(that.itsSolvableMask);
    itsDomainGrid = that.itsDomainGrid;
    itsValues = that.itsValues;
    itsDefaultValue = that.itsDefaultValue;
    itsScaleDomain = that.itsScaleDomain;
    itsDirty = that.itsDirty;
  }
  return *this;
}

const ParmValue& ParmValueSet::getFirstParmValue() const {
  return itsValues.empty() ? itsDefaultValue : *itsValues[0];
}

void ParmValueSet::setSolveGrid(const Grid& solveGrid) {
  // If the grid is empty, we must add the entire solve grid.
  if (itsDomainGrid.isDefault()) {
    createValues(solveGrid);
  } else {
    // If the entire solve grid is part of the values, check if the
    // grid matches.
    if (itsDomainGrid.getBoundingBox().contains(solveGrid.getBoundingBox())) {
      checkGrid(solveGrid);
    } else {
      // Part of the solve grid does not exist, so they need to be added.
      addValues(solveGrid);
    }
  }
}

void ParmValueSet::createValues(const Grid& solveGrid) {
  assert(itsValues.empty());
  // If the ParmValue represents coefficients, copy it as often as needed.
  if (itsType != ParmValue::Scalar) {
    itsDomainGrid = solveGrid;
    const Axis& xaxis = *itsDomainGrid[0];
    const Axis& yaxis = *itsDomainGrid[1];
    unsigned int nrx = itsDomainGrid.nx();
    unsigned int nry = itsDomainGrid.ny();
    itsValues.reserve(nrx * nry);
    for (unsigned int iy = 0; iy < nry; ++iy) {
      for (unsigned int ix = 0; ix < nrx; ++ix) {
        ParmValue::ShPtr pval(new ParmValue(itsDefaultValue));
        itsValues.push_back(pval);
        if (!itsScaleDomain.empty()) {
          // Rescale from itsScaleDomain to new domain.
          pval->rescale(xaxis.lower(ix), xaxis.upper(ix), yaxis.lower(iy),
                        yaxis.upper(iy), itsScaleDomain);
        }
      }
    }
  } else {
    // Otherwise it is an array of scalar values, so form the array.
    Array<double> values(IPosition(2, solveGrid.nx(), solveGrid.ny()));
    // Set it to the default value.
    values = itsDefaultValue.getValues().data()[0];
    ParmValue::ShPtr newVal(new ParmValue());
    newVal->setScalars(solveGrid, values);
    itsValues.push_back(newVal);
    itsDomainGrid = Grid(vector<Box>(1, solveGrid.getBoundingBox()));
  }
}

void ParmValueSet::checkGrid(const Grid& solveGrid) {
  // Check if the solve grid intervals match the domain grid.
  // If the values represent coefficients, the domain grid is the final grid
  // which should match the solve grid.
  if (itsType != ParmValue::Scalar) {
    assert(itsDomainGrid.checkIntervals(solveGrid));
  } else {
    // Each ParmValue has its own grid which has to be checked.
    if (itsValues.size() == 1) {
      // Only one value, so its grid should match.
      assert(itsValues[0]->getGrid().checkIntervals(solveGrid));
    } else {
      // The domain grid is split, so check each part with the corresponding
      // subset of the solve grid.
      for (unsigned int i = 0; i < itsDomainGrid.size(); ++i) {
        assert(itsValues[i]->getGrid().checkIntervals(solveGrid));
      }
    }
  }
}

void ParmValueSet::addValues(const Grid& solveGrid) {
  // Add values and extend the domain grid.
  // If the values represent coefficients, the domain grid is the final grid
  // which should match the solve grid.
  if (itsType != ParmValue::Scalar) {
    addCoeffValues(solveGrid);
  } else {
    // The values is an array of scalars.
    // For now only a single array can be handled.
    // If there are multiple arrays, it is (too) hard to decide which one
    // gets extended or if a new array has to be added.
    assert(itsValues.size() == 1);
    ParmValue& value = *itsValues[0];
    int sx1, ex1, sx2, ex2, sy1, sy2, ey1, ey2;
    Axis::ShPtr xaxis = value.getGrid().getAxis(0)->combine(
        *solveGrid.getAxis(0), sx1, ex1, sx2, ex2);
    Axis::ShPtr yaxis = value.getGrid().getAxis(1)->combine(
        *solveGrid.getAxis(1), sy1, ey1, sy2, ey2);
    Grid newGrid(xaxis, yaxis);
    Array<double> newValues(IPosition(2, newGrid.nx(), newGrid.ny()));
    // Copy the old values.
    newValues(IPosition(2, sx1, sy1), IPosition(2, ex1 - 1, ey1 - 1)) =
        value.getValues();
    // Fill in the other values.
    // In the extreme case the old values are in the middle of the
    // new values, so all sides have to be filled.
    // The values before are filled with the first old value, the values
    // after with the last old value.
    // In this way we achieve that new solutions are initialized with
    // existing ones.
    // First copy the values before and after the x-part of the old values.
    for (int iy = sy1; iy < ey1; ++iy) {
      for (int ix = sx2; ix < sx1; ++ix) {
        newValues(IPosition(2, ix, iy)) = newValues(IPosition(2, sx1, iy));
      }
      for (int ix = ex1; ix < ex2; ++ix) {
        newValues(IPosition(2, ix, iy)) = newValues(IPosition(2, ex1 - 1, iy));
      }
    }
    // Now copy the values before and after the y-part of the old values.
    int nrx = newValues.shape()[0];
    for (int iy = sy2; iy < sy1; ++iy) {
      for (int ix = 0; ix < nrx; ++ix) {
        newValues(IPosition(2, ix, iy)) = newValues(IPosition(2, ix, sy1));
      }
    }
    for (int iy = ey1; iy < ey2; ++iy) {
      for (int ix = 0; ix < nrx; ++ix) {
        newValues(IPosition(2, ix, iy)) = newValues(IPosition(2, ix, ey1 - 1));
      }
    }
    value.setScalars(newGrid, newValues);
    itsDomainGrid = Grid(vector<Box>(1, newGrid.getBoundingBox()));
  }
}

void ParmValueSet::addCoeffValues(const Grid& solveGrid) {
  // Combine the domain grid and the solve grid to form the new domain grid.
  // The values sx and ex give for each axis the start and end of the old
  // axes in the new one.
  int sx1, ex1, sx2, ex2, sy1, sy2, ey1, ey2;
  Axis::ShPtr xaxis = itsDomainGrid.getAxis(0)->combine(*solveGrid.getAxis(0),
                                                        sx1, ex1, sx2, ex2);
  Axis::ShPtr yaxis = itsDomainGrid.getAxis(1)->combine(*solveGrid.getAxis(1),
                                                        sy1, ey1, sy2, ey2);
  Grid newGrid(xaxis, yaxis);
  // Now copy existing parm values and insert new ones as necessary.
  // Take care that the ParmValues are in the correct order
  // (i.e. in order of the cells in the new grid).
  int nx = xaxis->size();
  int ny = yaxis->size();
  vector<ParmValue::ShPtr> newValues(nx * ny);
  // Copy the old values.
  assert((ex1 - sx1) * (ey1 - sy1) == int(itsValues.size()));
  const ParmValue::ShPtr* oldValues = &(itsValues[0]);
  for (int iy = 0; iy < ey1 - sy1; ++iy) {
    for (int ix = 0; ix < ex1 - sx1; ++ix) {
      newValues[sx1 + ix + (sy1 + iy) * nx] = *oldValues++;
    }
  }
  // Fill in the other values.
  // In the extreme case the old values are in the middle of the
  // new values, so all sides have to be filled.
  // The values before are filled with the first old value, the values
  // after with the last old value.
  // In this way we achieve that new solutions are initialized with
  // existing ones.
  // First copy the values before and after the x-part of the old values.
  for (int iy = sy1; iy < ey1; ++iy) {
    for (int ix = sx2; ix < sx1; ++ix) {
      newValues[ix + iy * nx] = copyParmCoeff(newValues[sx1 + iy * nx]);
    }
    for (int ix = ex1; ix < ex2; ++ix) {
      newValues[ix + iy * nx] = copyParmCoeff(newValues[ex1 - 1 + iy * nx]);
    }
  }
  // Now copy the values before and after the y-part of the old values.
  for (int iy = sy2; iy < sy1; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      newValues[ix + iy * nx] = copyParmCoeff(newValues[ix + sy1 * nx]);
    }
  }
  for (int iy = ey1; iy < ey2; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      newValues[ix + iy * nx] = copyParmCoeff(newValues[ix + (ey1 - 1) * nx]);
    }
  }
  // Use new values and new grid.
  itsValues.swap(newValues);
  itsDomainGrid = newGrid;
}

ParmValue::ShPtr ParmValueSet::copyParmCoeff(const ParmValue::ShPtr& pval) {
  ParmValue::ShPtr newpval(new ParmValue(*pval));
  newpval->clearRowId();
  return newpval;
}

}  // namespace parmdb
}  // namespace dp3
