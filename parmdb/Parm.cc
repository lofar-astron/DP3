// ParmCache.cc: A class dealing with caching and handling ParmDB entries
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Parm.h"
#include "ParmCache.h"
#include "AxisMapping.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

using namespace casacore;
using namespace std;

namespace dp3 {
namespace parmdb {

Parm::Parm(ParmCache& cache, ParmId parmid)
    : itsCache(&cache), itsParmId(parmid) {}

Parm::Parm(ParmCache& cache, const string& name)
    : itsCache(&cache), itsParmId(cache.getParmSet().find(name)) {}

void Parm::setSolveGrid(const Grid& solveGrid) {
  itsCache->setSolveGrid(itsParmId, solveGrid);
  itsSolveGrid = solveGrid;
  // Solve grid is set, so calculate the perturbations.
  calcPerturbations();
}

unsigned int Parm::getCoeffSize(bool useMask) {
  const ParmValueSet& pvset = itsCache->getValueSet(itsParmId);
  // For a scalar array, only one coeff is used.
  if (pvset.getType() == ParmValue::Scalar) {
    return 1;
  }
  // Use the first ParmValue.
  const ParmValue& pv = pvset.getFirstParmValue();
  const Array<double>& values = pv.getValues();
  const Array<bool>& mask = pvset.getSolvableMask();
  if (!useMask || mask.size() == 0) {
    return values.size();
  }
  return ntrue(mask);
}

vector<double> Parm::getCoeff(const Grid::Location& where, bool useMask) {
  assert(!itsSolveGrid.isDefault());
  const ParmValueSet& pvset = itsCache->getValueSet(itsParmId);
  // Find the location in the ParmValueSet grid given the location in
  // the solve grid.
  unsigned int cellId = GridMapping::findCellId(
      itsCache->getAxisMappingCache(), where, itsSolveGrid, pvset.getGrid());
  const ParmValue& pv = pvset.getParmValue(cellId);
  if (pvset.getType() != ParmValue::Scalar) {
    return copyValues(pv.getValues(), pvset.getSolvableMask(), useMask);
  }
  // An array of scalar values; get the right one.
  cellId = GridMapping::findCellId(itsCache->getAxisMappingCache(), where,
                                   itsSolveGrid, pv.getGrid());
  return vector<double>(1, pv.getValues().data()[cellId]);
}

vector<double> Parm::getErrors(const Grid::Location& where, bool useMask) {
  const ParmValueSet& pvset = itsCache->getValueSet(itsParmId);
  // Find the location in the ParmValueSet grid given the location in
  // the solve grid.
  unsigned int cellId = GridMapping::findCellId(
      itsCache->getAxisMappingCache(), where, itsSolveGrid, pvset.getGrid());
  const ParmValue& pv = pvset.getParmValue(cellId);
  if (!pv.hasErrors()) {
    return vector<double>();
  }
  if (pvset.getType() != ParmValue::Scalar) {
    return copyValues(pv.getErrors(), pvset.getSolvableMask(), useMask);
  }
  // An array of scalar values; get the right one.
  cellId = GridMapping::findCellId(itsCache->getAxisMappingCache(), where,
                                   itsSolveGrid, pv.getGrid());
  return vector<double>(1, pv.getErrors().data()[cellId]);
}

vector<double> Parm::copyValues(const Array<double>& values,
                                const Array<Bool>& mask, bool useMask) {
  assert(values.contiguousStorage());
  if (!useMask || mask.size() == 0) {
    return vector<double>(values.cbegin(), values.cend());
  }
  assert(values.shape().isEqual(mask.shape()) && mask.contiguousStorage());
  vector<double> solvCoeff;
  solvCoeff.reserve(values.size());
  const double* valp = values.data();
  const bool* maskp = mask.data();
  for (unsigned int i = 0; i < values.size(); ++i) {
    if (maskp[i]) {
      solvCoeff.push_back(valp[i]);
    }
  }
  return solvCoeff;
}

void Parm::setCoeff(const Grid::Location& where, const double* newValues,
                    unsigned int nvalues, const double* newErrors,
                    bool useMask) {
  assert(!itsSolveGrid.isDefault());
  ParmValueSet& pvset = itsCache->getValueSet(itsParmId);
  pvset.setDirty();
  const Array<bool>& mask = pvset.getSolvableMask();
  unsigned int cellId = GridMapping::findCellId(
      itsCache->getAxisMappingCache(), where, itsSolveGrid, pvset.getGrid());
  bool cell0changed = (cellId == 0);
  ParmValue& pv = pvset.getParmValue(cellId);
  Array<double>& values(pv.getValues());
  if (newErrors && !pv.hasErrors()) {
    Array<double> err(values.shape());
    err = 0;
    pv.setErrors(err);
  }
  Array<double>& errors(pv.getErrors());
  assert(values.contiguousStorage());
  if (pvset.getType() != ParmValue::Scalar) {
    if (!useMask || mask.size() == 0) {
      // Coefficients without a mask; copy all.
      assert(nvalues == values.size());
      copy(newValues, newValues + nvalues, values.cbegin());
      if (newErrors) {
        assert(nvalues == errors.size());
        copy(newErrors, newErrors + nvalues, errors.cbegin());
      }
    } else {
      // Only copy values where mask is true.
      assert(values.shape().isEqual(mask.shape()) && mask.contiguousStorage());
      double* valp = values.data();
      double* errp = 0;
      if (newErrors) {
        errp = errors.data();
      }
      const bool* maskp = mask.data();
      unsigned int inx = 0;
      for (unsigned int i = 0; i < values.size(); ++i) {
        if (maskp[i]) {
          valp[i] = newValues[inx];
          if (newErrors) {
            errp[i] = newErrors[inx];
          }
          ++inx;
        }
      }
      assert(inx == nvalues);  // make sure everything is copied
    }
  } else {
    // A scalar array; copy to the correct cell.
    assert(nvalues == 1);
    cellId = GridMapping::findCellId(itsCache->getAxisMappingCache(), where,
                                     itsSolveGrid, pv.getGrid());
    values.data()[cellId] = newValues[0];
    if (newErrors) {
      errors.data()[cellId] = newErrors[0];
    }

    cell0changed = (cellId == 0);
  }
  // Coefficients have changed, so recalculate the perturbations.
  if (cell0changed) {
    calcPerturbations();
  }
}

void Parm::revertCoeff() {
  assert(!itsSolveGrid.isDefault());
  // Moet ik nog over nadenken.
  throw runtime_error("revertCoeff is not implemented yet");
  // Coefficients have changed, so recalculate the perturbations.
  calcPerturbations();
}

void Parm::calcPerturbations() {
  const ParmValueSet& pvset = itsCache->getValueSet(itsParmId);
  const ParmValue& pv = pvset.getFirstParmValue();
  if (pvset.getType() != ParmValue::Scalar) {
    itsPerturbations =
        copyValues(pv.getValues(), pvset.getSolvableMask(), true);
  } else {
    itsPerturbations.resize(1);
    itsPerturbations[0] = pv.getValues().data()[0];
  }
  double perturbation = pvset.getPerturbation();
  for (vector<double>::iterator iter = itsPerturbations.begin();
       iter != itsPerturbations.end(); ++iter) {
    if (pvset.getPertRel() && std::abs(*iter) > 1e-10) {
      *iter *= perturbation;
    } else {
      *iter = perturbation;
    }
  }
}

void Parm::getResult(vector<Array<double>>& result, const Grid& predictGrid,
                     bool perturb) {
  if (!perturb || itsPerturbations.empty()) {
    // No perturbed values need to be calculated.
    if (result.empty()) {
      result.resize(1);
    }
    getResult(result[0], predictGrid);
  } else {
    // Perturbed values need to be calculated. Make room for them.
    result.resize(itsPerturbations.size() + 1);
    // If no values found, return an empty array.
    ParmValueSet& pvset = itsCache->getValueSet(itsParmId);
    if (!pvset.empty()) {
      if (pvset.getType() != ParmValue::Scalar) {
        // It is a funklet, so evaluate it.
        getResultCoeff(&(result[0]), predictGrid, pvset, itsPerturbations,
                       itsCache->getAxisMappingCache());
      } else {
        // We have scalar values, thus only one perturbed value.
        // First get result and add perturbed value to it.
        assert(itsPerturbations.size() == 1);
        getResult(result[0], predictGrid);
        result[1].resize(result[0].shape());
        result[1] = result[0] + itsPerturbations[0];
      }
    }
  }
}

void Parm::getResult(Array<double>& result, const Grid& predictGrid,
                     bool emptyResult) {
  // Get the values.
  ParmValueSet& pvset = itsCache->getValueSet(itsParmId);
  if (emptyResult && pvset.empty()) {
    result.resize();
    return;
  }
  if (pvset.getType() != ParmValue::Scalar) {
    // It is a funklet, so evaluate it.
    getResultCoeff(&result, predictGrid, pvset, vector<double>(),
                   itsCache->getAxisMappingCache());
  } else if (pvset.getGrid().size() == 1) {
    // Optimize for the often occurring case of a single ParmValue object.
    const ParmValue& pval = pvset.getFirstParmValue();
    if (pval.getGrid().size() == 1) {
      // Only a single value, so size the array accordingly.
      result.resize(IPosition(2, 1, 1));
      result = pval.getValues();
    } else {
      // There are multiple values, so use the ParmValue's grid.
      getResultScalar(result, predictGrid, pval,
                      itsCache->getAxisMappingCache());
    }
  } else {
    // The hardest case; multiple ParmValues, possibly each with its own grid.
    getResultScalar(result, 0, predictGrid, pvset,
                    itsCache->getAxisMappingCache());
  }
}

void Parm::getResultCoeff(Array<double>* resultVec, const Grid& predictGrid,
                          const ParmValueSet& pvset,
                          const vector<double>& perturbations,
                          AxisMappingCache& axisMappingCache) {
  Array<double>& result = *resultVec;
  const Axis& paxisx = *predictGrid.getAxis(0);
  const Axis& paxisy = *predictGrid.getAxis(1);
  const Axis& daxisx = *pvset.getGrid().getAxis(0);
  const Axis& daxisy = *pvset.getGrid().getAxis(1);
  // Get the x and y axis mapping of predict grid to domain grid.
  const AxisMapping& mapx = axisMappingCache.get(paxisx, daxisx);
  const AxisMapping& mapy = axisMappingCache.get(paxisy, daxisy);
  int nrdx = daxisx.size();
  const double* cenx = mapx.getScaledCenters();
  const double* ceny = mapy.getScaledCenters();
  // First calculate the main result.
  // Size the array as needed and get an iterator for it.
  result.resize(IPosition(2, paxisx.size(), paxisy.size()));
  Array<double>::iterator resultIter = result.begin();
  const double* pvaly = ceny;
  // Loop over all cells of the predict y-axis.
  for (AxisMapping::const_iterator ity = mapy.begin(); ity != mapy.end();
       ++ity) {
    int inxy = *ity * nrdx;
    double valy = *pvaly++;
    const double* pvalx = cenx;
    // Loop over all cells of the predict x-axis.
    for (AxisMapping::const_iterator itx = mapx.begin(); itx != mapx.end();
         ++itx) {
      double valx = *pvalx++;
      // Get the coefficients.
      // If no ParmValues are found, take default one.
      const Array<double>* carr = &(pvset.getDefParmValue().getValues());
      if (*itx + inxy < int(pvset.size())) {
        carr = &(pvset.getParmValue(*itx + inxy).getValues());
      }
      const double* coeff = carr->data();
      int nrcx = carr->shape()[0];
      int nrcy = carr->shape()[1];
      // Calculate sigma(c[i,j] * x^i * y^j)
      double y = 1;
      double val = 0;
      for (int j = 0; j < nrcy; ++j) {
        double subval = 0;
        for (int i = nrcx - 1; i > 0; i--) {
          subval += coeff[i];
          subval *= valx;
        }
        subval += coeff[0];
        val += y * subval;
        y *= valy;
        coeff += nrcx;
      }
      *resultIter = val;
      ++resultIter;
    }
  }
  // Now calculate all perturbed values if needed.
  if (!perturbations.empty()) {
    vector<double> pertCoeff(perturbations.size(), 0.);
    for (unsigned int ip = 0; ip < perturbations.size(); ++ip) {
      pertCoeff[ip] = perturbations[ip];
      Array<double>& result = resultVec[ip + 1];
      // Size the array as needed and get an iterator for it.
      result.resize(IPosition(2, paxisx.size(), paxisy.size()));
      Array<double>::iterator resultIter = result.begin();
      const double* pvaly = ceny;
      // Loop over all cells of the predict y-axis.
      for (AxisMapping::const_iterator ity = mapy.begin(); ity != mapy.end();
           ++ity) {
        int inxy = *ity * nrdx;
        double valy = *pvaly++;
        const double* pvalx = cenx;
        // Loop over all cells of the predict x-axis.
        for (AxisMapping::const_iterator itx = mapx.begin(); itx != mapx.end();
             ++itx) {
          double valx = *pvalx++;
          // Get the coefficients.
          // If no ParmValues are found, take default one.
          const Array<double>* carr = &(pvset.getDefParmValue().getValues());
          if (*itx + inxy < int(pvset.size())) {
            carr = &(pvset.getParmValue(*itx + inxy).getValues());
          }
          assert(carr->size() == perturbations.size());
          const double* coeff = carr->data();
          const double* pcoeff = &(pertCoeff[0]);
          int nrcx = carr->shape()[0];
          int nrcy = carr->shape()[1];
          // Calculate sigma(c[i,j] * x^i * y^j)
          double y = 1;
          double val = 0;
          for (int j = 0; j < nrcy; ++j) {
            double subval = 0;
            for (int i = nrcx - 1; i > 0; i--) {
              subval += coeff[i] + pcoeff[i];
              subval *= valx;
            }
            subval += coeff[0] + pcoeff[0];
            val += y * subval;
            y *= valy;
            coeff += nrcx;
            pcoeff += nrcx;
          }
          *resultIter = val;
          ++resultIter;
        }
      }
      pertCoeff[ip] = 0.;
    }
  }
}

void Parm::getResultScalar(Array<double>& result, const Grid& predictGrid,
                           const ParmValue& pval,
                           AxisMappingCache& axisMappingCache) {
  const Axis& paxisx = *predictGrid.getAxis(0);
  const Axis& paxisy = *predictGrid.getAxis(1);
  const Axis& daxisx = *pval.getGrid().getAxis(0);
  const Axis& daxisy = *pval.getGrid().getAxis(1);
  // Get the x and y axis mapping of predict grid to domain grid.
  const AxisMapping& mapx = axisMappingCache.get(paxisx, daxisx);
  const AxisMapping& mapy = axisMappingCache.get(paxisy, daxisy);
  int nrdx = daxisx.size();
  const double* data = pval.getValues().data();
  // Size the array as needed and get an iterator for it.
  result.resize(IPosition(2, paxisx.size(), paxisy.size()));
  Array<double>::iterator resultIter = result.begin();
  // Loop over all cells of the predict y-axis.
  for (AxisMapping::const_iterator ity = mapy.begin(); ity != mapy.end();
       ++ity) {
    int inxy = *ity * nrdx;
    // Loop over all cells of the predict x-axis.
    for (AxisMapping::const_iterator itx = mapx.begin(); itx != mapx.end();
         ++itx) {
      *resultIter = data[*itx + inxy];
      ++resultIter;
    }
  }
}

void Parm::getResultScalar(Array<double>& result, Array<double>* errors,
                           const Grid& predictGrid, const ParmValueSet& pvset,
                           AxisMappingCache& axisMappingCache) {
  const Axis& paxisx = *predictGrid.getAxis(0);
  const Axis& paxisy = *predictGrid.getAxis(1);
  const Axis& saxisx = *pvset.getGrid().getAxis(0);
  const Axis& saxisy = *pvset.getGrid().getAxis(1);
  // Get the x and y axis mapping of predict grid to the set's domain grid.
  const AxisMapping& mapx = axisMappingCache.get(paxisx, saxisx);
  const AxisMapping& mapy = axisMappingCache.get(paxisy, saxisy);
  int nrsx = saxisx.size();
  // Size the array as needed and get a raw pointer to the result data.
  result.resize(IPosition(2, paxisx.size(), paxisy.size()));
  bool deleteRes;
  double* resData = result.getStorage(deleteRes);
  bool deleteErr;
  double* errData = 0;
  if (errors) {
    errors->resize(result.shape());
    *errors = -1;
    errData = errors->getStorage(deleteErr);
  }
  int nrx = result.shape()[0];
  // Loop through the cells of pvset's grid.
  // Each cell is a ParmValue with its own grid.
  // Fill a part of the result from the ParmValue.
  const vector<int>& bordersx = mapx.getBorders();
  const vector<int>& bordersy = mapy.getBorders();
  int sty = 0;
  for (unsigned int iy = 0; iy < bordersy.size(); ++iy) {
    int inxy = nrsx * mapy[sty];
    int stx = 0;
    for (unsigned int ix = 0; ix < bordersx.size(); ++ix) {
      const ParmValue& pval = pvset.getParmValue(mapx[stx] + inxy);
      fillArrayPV(resData, nrx, stx, sty, bordersx[ix], bordersy[iy],
                  pval.getValues().data(), pval, predictGrid);
      if (errors && pval.hasErrors()) {
        fillArrayPV(errData, nrx, stx, sty, bordersx[ix], bordersy[iy],
                    pval.getErrors().data(), pval, predictGrid);
      }
      stx = bordersx[ix];
    }
    sty = bordersy[iy];
  }
  result.putStorage(resData, deleteRes);
  if (errors) {
    errors->putStorage(errData, deleteErr);
  }
}

void Parm::fillArrayPV(double* resData, int nrx, int stx, int sty, int endx,
                       int endy, const double* data, const ParmValue& pval,
                       const Grid& predictGrid) {
  // Get the axes of predict grid and domain grid.
  const Axis& paxisx = *predictGrid.getAxis(0);
  const Axis& paxisy = *predictGrid.getAxis(1);
  const Axis& daxisx = *pval.getGrid().getAxis(0);
  const Axis& daxisy = *pval.getGrid().getAxis(1);
  int nrdx = daxisx.size();
  // Loop through all relevant cells of the predict grid.
  // Find the corresponding cell in the domaingrid and copy its value.
  int inxy = 0;
  int inxx = 0;
  for (int iy = sty; iy < endy; ++iy) {
    // Set result pointer to the beginning of this chunk.
    double* rData = resData + iy * nrx + stx;
    inxy = daxisy.locate(paxisy.center(iy), true, inxy);
    const double* pData = data + inxy * nrdx;
    for (int ix = stx; ix < endx; ++ix) {
      inxx = daxisx.locate(paxisx.center(ix), true, inxx);
      *rData++ = pData[inxx];
    }
  }
}

}  // namespace parmdb
}  // namespace dp3
