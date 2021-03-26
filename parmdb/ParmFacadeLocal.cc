// ParmFacadeLocal.cc: Object access the parameter database
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParmFacadeLocal.h"
#include "ParmDB.h"
#include "ParmSet.h"
#include "ParmCache.h"
#include "Parm.h"

#include "../common/StringTools.h"

#include <casacore/casa/Utilities/Regex.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <set>
#include <string>

using namespace std;
using namespace casacore;

// Create tParmFacade.in_mep with parmdb using:
//   create tablename='tParmFacade.in_mep'
//   add parm1 domain=[1,4,5,10],values=2
//   add parm2 domain=[1,4,5,10],values=[2,0.1],nx=2
//   add parm3 type='expression',expression='parm1*parm2'

namespace dp3 {
namespace parmdb {

ParmFacadeLocal::ParmFacadeLocal(const string& tableName, bool create)
    : itsPDB(ParmDBMeta("casa", tableName), create) {}

ParmFacadeLocal::~ParmFacadeLocal() {}

vector<double> ParmFacadeLocal::getRange(const string& parmNamePattern) const {
  string pp = parmNamePattern;
  if (pp.empty()) {
    pp = "*";
  }
  Box dom = itsPDB.getRange(pp);
  vector<double> res(4);
  res[0] = dom.lowerX();
  res[1] = dom.upperX();
  res[2] = dom.lowerY();
  res[3] = dom.upperY();
  return res;
}

// Get all parameter names in the table.
vector<string> ParmFacadeLocal::getNames(const string& parmNamePattern,
                                         Bool includeDefaults) const {
  string pp = parmNamePattern;
  if (pp.empty()) {
    pp = "*";
  }
  vector<string> names(itsPDB.getNames(pp));
  if (!includeDefaults) {
    return names;
  }
  // Get possible unused names of default values.
  ParmMap parmset;
  itsPDB.getDefValues(parmset, pp);
  std::set<std::string> nameSet(names.begin(), names.end());
  for (ParmMap::const_iterator iter = parmset.begin(); iter != parmset.end();
       ++iter) {
    if (nameSet.find(iter->first) == nameSet.end()) {
      names.push_back(iter->first);
    }
  }
  return names;
}

vector<string> ParmFacadeLocal::getDefNames(
    const string& parmNamePattern) const {
  string pp = parmNamePattern;
  if (pp.empty()) {
    pp = "*";
  }
  ParmMap parmset;
  itsPDB.getDefValues(parmset, pp);
  vector<string> names;
  names.reserve(parmset.size());
  for (ParmMap::const_iterator iter = parmset.begin(); iter != parmset.end();
       iter++) {
    names.push_back(iter->first);
  }
  return names;
}

Record ParmFacadeLocal::getDefValues(const string& parmNamePattern) const {
  string pp = parmNamePattern;
  if (pp.empty()) {
    pp = "*";
  }
  ParmMap parmset;
  itsPDB.getDefValues(parmset, pp);
  Record result;
  for (ParmMap::const_iterator iter = parmset.begin(); iter != parmset.end();
       iter++) {
    // Define name to default value.
    result.define(iter->first, iter->second.getDefParmValue().getValues());
  }
  return result;
}

void ParmFacadeLocal::addDefValues(const Record& rec, bool check) {
  itsPDB.lock();
  for (uInt i = 0; i < rec.nfields(); ++i) {
    addDefValue(rec.name(i), rec.subRecord(i), check);
  }
  itsPDB.unlock();
}

void ParmFacadeLocal::addDefValue(const string& parmName, const Record& rec,
                                  bool check) {
  // Create the default value.
  Array<double> value(rec.toArrayDouble("value"));
  ParmValue pval(value.data()[0]);
  // Get the type.
  int type = -1;
  if (rec.isDefined("type")) {
    type = getType(rec.asString("type"));
  }
  // Set coefficients if no scalar or multiple values.
  if (type > 0 || value.size() > 1) {
    pval.setCoeff(value);
    // In this case default type is Polc.
    if (type < 0) type = ParmValue::Polc;
  } else {
    type = ParmValue::Scalar;
  }
  double pert = 1e-6;
  if (rec.isDefined("perturbation")) {
    pert = rec.asDouble("perturbation");
  }
  bool pertRel = true;
  if (rec.isDefined("pertrel")) {
    pertRel = rec.asBool("pertrel");
  }
  // Create the default value.
  ParmValueSet pvset(pval, ParmValue::FunkletType(type), pert, pertRel);
  if (rec.isDefined("mask")) {
    Array<bool> mask(rec.toArrayBool("mask"));
    if (mask.size() > 0) {
      pvset.setSolvableMask(mask);
    }
  }
  itsPDB.putDefValue(parmName, pvset, check);
}

void ParmFacadeLocal::deleteDefValues(const string& parmNamePattern) {
  itsPDB.deleteDefValues(parmNamePattern);
}

Record ParmFacadeLocal::getValues(const string& parmNamePattern, double freqv1,
                                  double freqv2, double freqStep, double timev1,
                                  double timev2, double timeStep,
                                  bool asStartEnd, bool includeDefaults) {
  // Use default step values if needed.
  if (freqStep <= 0) {
    freqStep = itsPDB.getDefaultSteps()[0];
  }
  if (timeStep <= 0) {
    timeStep = itsPDB.getDefaultSteps()[1];
  }
  int nfreq, ntime;
  if (asStartEnd) {
    nfreq = std::max(1, int((freqv2 - freqv1) / freqStep + 0.5));
    ntime = std::max(1, int((timev2 - timev1) / timeStep + 0.5));
  } else {
    nfreq = std::max(1, int(freqv2 / freqStep + 0.5));
    ntime = std::max(1, int(timev2 / timeStep + 0.5));
  }
  // Create the predict grid.
  Axis::ShPtr axisx(new RegularAxis(freqv1, freqv2, nfreq, asStartEnd));
  Axis::ShPtr axisy(new RegularAxis(timev1, timev2, ntime, asStartEnd));
  return doGetValues(parmNamePattern, Grid(axisx, axisy), includeDefaults);
}

Record ParmFacadeLocal::getValues(const string& parmNamePattern,
                                  const vector<double>& freqv1,
                                  const vector<double>& freqv2,
                                  const vector<double>& timev1,
                                  const vector<double>& timev2, bool asStartEnd,
                                  bool includeDefaults) {
  // Create the predict grid.
  Axis::ShPtr axisx(new OrderedAxis(freqv1, freqv2, asStartEnd));
  Axis::ShPtr axisy(new OrderedAxis(timev1, timev2, asStartEnd));
  return doGetValues(parmNamePattern, Grid(axisx, axisy), includeDefaults);
}

void ParmFacadeLocal::clearTables() { itsPDB.clearTables(); }

void ParmFacadeLocal::flush(bool fsync) { itsPDB.flush(fsync); }

void ParmFacadeLocal::lock(bool lockForWrite) { itsPDB.lock(lockForWrite); }

void ParmFacadeLocal::unlock() { itsPDB.unlock(); }

vector<double> ParmFacadeLocal::getDefaultSteps() const {
  return itsPDB.getDefaultSteps();
}

void ParmFacadeLocal::setDefaultSteps(const vector<double>& steps) {
  itsPDB.setDefaultSteps(steps);
}

void ParmFacadeLocal::addValues(const Record& rec) {
  itsPDB.lock();
  for (uInt i = 0; i < rec.nfields(); ++i) {
    addValue(rec.name(i), rec.subRecord(i));
  }
  itsPDB.unlock();
}

void ParmFacadeLocal::addValue(const string& parmName, const Record& rec) {
  // Get the new values and shape. Make 2D if a single value.
  Array<double> values(rec.toArrayDouble("values"));
  if (values.size() == 1) {
    values.reference(values.reform(IPosition(2, 1, 1)));
  }
  // Check if values is 2D.
  if (values.ndim() != 2)
    throw std::runtime_error("Values of parameter " + parmName +
                             " must be a scalar or 2-dim array");
  const IPosition& nshape = values.shape();
  // Form the grid from the record fields.
  Grid grid(record2Grid(rec));
  Box domain(grid.getBoundingBox());
  // Read the values of parameter and domain from the ParmDB.
  ParmSet parmset;
  ParmId parmid = parmset.addParm(itsPDB, parmName);
  ParmCache cache(parmset, domain);
  ParmValueSet& pvset = cache.getValueSet(parmid);
  // Assure no value exists yet.
  if (pvset.size() != 0)
    throw std::runtime_error("Value for this parameter/domain " + parmName +
                             " already exists");
  // Check if the parm already exists (for another domain).
  // If so, only the values can be set (but not the meta info).
  // Check type and shape.
  bool isOldParm = cache.getParmSet().isInParmDB(parmid);
  // Set current values as the meta defaults.
  ParmValue defval(pvset.getFirstParmValue());
  int type = pvset.getType();
  double pert = pvset.getPerturbation();
  bool pertrel = pvset.getPertRel();
  Array<Bool> mask = pvset.getSolvableMask();
  IPosition shape = defval.getValues().shape();
  // Get possible new meta values.
  if (!isOldParm) {
    shape.resize(0);
    shape = nshape;
    if (rec.isDefined("type")) {
      type = getType(rec.asString("type"));
    }
    if (rec.isDefined("mask") && type != ParmValue::Scalar) {
      mask.reference(rec.toArrayBool("mask"));
    }
  } else {
    if (rec.isDefined("type")) {
      int newType = getType(rec.asString("type"));
      if (newType != type)
        throw std::runtime_error(
            "New type " + std::to_string(newType) + " of parameter " +
            parmName + " mismatches existing type " + std::to_string(type));
    }
  }
  // Check sizes.
  if (type == ParmValue::Scalar) {
    // Turn a single axis value into a regular axis.
    Axis::ShPtr freqAxis = grid[0];
    Axis::ShPtr timeAxis = grid[1];
    if (grid.nx() == 1 && nshape[0] > 1) {
      freqAxis = makeAxis(freqAxis->centers(), freqAxis->widths(), nshape[0]);
    }
    if (grid.ny() == 1 && nshape[1] > 1) {
      timeAxis = makeAxis(timeAxis->centers(), timeAxis->widths(), nshape[1]);
    }
    grid = Grid(freqAxis, timeAxis);
    if (int(grid.nx()) != nshape[0] || int(grid.ny()) != nshape[1])
      throw std::runtime_error(
          "Mismatch in shape of coeff and grid for scalar parameter " +
          parmName);
  } else {
    if (!nshape.isEqual(shape))
      throw std::runtime_error(
          "Non-scalar parameter " + parmName +
          " is used before; coeff shape cannot be changed");
    if (grid.size() != 1)
      throw std::runtime_error("Grid of non-scalar parameter " + parmName +
                               " should contain 1 element");
    if (!mask.empty()) {
      if (!nshape.isEqual(mask.shape()))
        throw std::runtime_error("Coeff and mask of non-scalar parameter " +
                                 parmName + " have different shapes");
    }
  }
  shape.resize(0);
  shape = nshape;
  // Create the parm value.
  ParmValue::ShPtr pval(new ParmValue);
  if (type == ParmValue::Scalar) {
    pval->setScalars(grid, values);
  } else {
    pval->setCoeff(values);
  }
  // Set the errors if given.
  if (rec.isDefined("errors")) {
    Array<double> errs = rec.toArrayDouble("errors");
    pval->setErrors(errs);
  }
  // Create the ParmValueSet.
  vector<ParmValue::ShPtr> pvalvec(1, pval);
  vector<Box> domains(1, domain);
  pvset = ParmValueSet(Grid(domains), pvalvec, defval,
                       ParmValue::FunkletType(type), pert, pertrel);
  pvset.setDirty();
  cache.flush();
}

void ParmFacadeLocal::deleteValues(const string& parmNamePattern, double freqv1,
                                   double freqv2, double timev1, double timev2,
                                   bool asStartEnd) {
  Box domain(freqv1, freqv2, timev1, timev2, asStartEnd);
  itsPDB.deleteValues(parmNamePattern, domain);
}

Record ParmFacadeLocal::doGetValues(const string& parmNamePattern,
                                    const Grid& predictGrid,
                                    Bool includeDefaults) {
  // Get all matching parm names.
  vector<string> names = getNames(parmNamePattern, includeDefaults);
  // The output is returned in a record.
  Record out;
  // Form the names to get.
  // The returned parmId should be the index.
  ParmSet parmSet;
  for (unsigned int i = 0; i < names.size(); ++i) {
    [[maybe_unused]] const ParmId index = parmSet.addParm(itsPDB, names[i]);
    assert(index == i);
  }
  const Axis& axisx = *predictGrid[0];
  const Axis& axisy = *predictGrid[1];
  unsigned int nfreq = axisx.size();
  unsigned int ntime = axisy.size();
  // Create and fill the cache for the given domain.
  Box domain(Point(axisx.lower(0), axisy.lower(0)),
             Point(axisx.upper(nfreq - 1), axisy.upper(ntime - 1)));
  ParmCache parmCache(parmSet, domain);
  // Now create the Parm object for each parm and get the values.
  Array<double> result;
  for (unsigned int i = 0; i < names.size(); ++i) {
    Parm parm(parmCache, i);
    parm.getResult(result, predictGrid, !includeDefaults);
    if (result.size() > 0) {
      // There is data in this domain.
      Record rec;
      // If the value is constant, the array has only one element.
      // In that case resize it to the full grid.
      if (result.nelements() == 1) {
        double dval = *result.data();
        result.resize(IPosition(2, axisx.size(), axisy.size()));
        result = dval;
      }
      rec.define("values", result);
      rec.define("freqs", Vector<double>(axisx.centers()));
      rec.define("times", Vector<double>(axisy.centers()));
      rec.define("freqwidths", Vector<double>(axisx.widths()));
      rec.define("timewidths", Vector<double>(axisy.widths()));
      out.defineRecord(names[i], rec);
    }
  }
  return out;
}

Record ParmFacadeLocal::getValuesGrid(const string& parmNamePattern,
                                      double freqv1, double freqv2,
                                      double timev1, double timev2,
                                      bool asStartEnd) {
  Box domain(freqv1, freqv2, timev1, timev2, asStartEnd);
  // Get all matching parm names.
  vector<string> names = getNames(parmNamePattern, false);
  // The output is returned in a record.
  Record out;
  // Form the names to get.
  // The returned parmId should be the index.
  ParmSet parmSet;
  for (unsigned int i = 0; i < names.size(); ++i) {
    [[maybe_unused]] const ParmId index = parmSet.addParm(itsPDB, names[i]);
    assert(index == i);
  }
  // Create and fill the cache for the given domain.
  ParmCache parmCache(parmSet, domain);
  // Now create the Parm object for each parm and get the values.
  Array<double> result;
  for (unsigned int i = 0; i < names.size(); ++i) {
    Grid grid(getGrid(parmCache.getValueSet(i), domain));
    if (!grid.isDefault()) {
      // There should be data in this domain.
      Parm parm(parmCache, i);
      parm.getResult(result, grid, true);
      if (result.size() > 0) {
        // There is data in this domain.
        Record rec;
        rec.define("values", result);
        rec.define("freqs", Vector<double>(grid[0]->centers()));
        rec.define("times", Vector<double>(grid[1]->centers()));
        rec.define("freqwidths", Vector<double>(grid[0]->widths()));
        rec.define("timewidths", Vector<double>(grid[1]->widths()));
        out.defineRecord(names[i], rec);
      }
    }
  }
  return out;
}

Grid ParmFacadeLocal::getGrid(const ParmValueSet& valueSet, const Box& domain) {
  Grid grid(valueSet.getGrid());
  if (valueSet.getType() == ParmValue::Scalar) {
    // For scalars the detailed grids have to be combined.
    vector<Grid> grids;
    grids.reserve(valueSet.size());
    for (unsigned int i = 0; i < valueSet.size(); ++i) {
      grids.push_back(valueSet.getParmValue(i).getGrid());
    }
    grid = Grid(grids, true);
  }
  return grid.subset(domain);
}

// Get coefficients, errors, and domains they belong to.
Record ParmFacadeLocal::getCoeff(const string& parmNamePattern, double freqv1,
                                 double freqv2, double timev1, double timev2,
                                 bool asStartEnd) {
  Box domain(freqv1, freqv2, timev1, timev2, asStartEnd);
  ParmMap result;
  itsPDB.getValues(result, parmNamePattern, domain);
  AxisMappingCache axesCache;
  Record rec;
  for (ParmMap::const_iterator iter = result.begin(); iter != result.end();
       ++iter) {
    const ParmValueSet& pvset = iter->second;
    Grid grid(getGrid(pvset, domain));
    if (!grid.isDefault()) {
      // Values found, so add them to the Record.
      if (pvset.getType() == ParmValue::Scalar) {
        // Scalars are put in 2D arrays.
        Array<double> result, errors;
        Parm::getResultScalar(result, &errors, grid, pvset, axesCache);
        const Axis& axisx = *grid[0];
        const Axis& axisy = *grid[1];
        Record vals;
        vals.define("values", result);
        vals.define("errors", errors);
        vals.define("freqs", Vector<double>(axisx.centers()));
        vals.define("times", Vector<double>(axisy.centers()));
        vals.define("freqwidths", Vector<double>(axisx.widths()));
        vals.define("timewidths", Vector<double>(axisy.widths()));
        rec.defineRecord(iter->first, vals);
      } else {
        // Funklets are put in 4D arrays.
        rec.defineRecord(iter->first, getFunkletCoeff(pvset));
      }
    }
  }
  return rec;
}

Record ParmFacadeLocal::getFunkletCoeff(const ParmValueSet& pvset) {
  // Get the grid of the funklets.
  const Grid& grid = pvset.getGrid();
  const Axis& axisx = *grid[0];
  const Axis& axisy = *grid[1];
  // Create a 4D array to hold the 2D coeff per funklet.
  // Its shape is formed by funklet shape and grid shape.
  IPosition shp = pvset.getParmValue(0).getValues().shape();
  shp.append(IPosition(2, axisx.size(), axisy.size()));
  Array<double> coeff(shp);
  Array<double> errors(shp);
  errors = -1;
  // Fill the arrays by iterating over them and each funklet.
  ArrayIterator<double> coeffIter(coeff, 2);
  ArrayIterator<double> errorIter(errors, 2);
  for (unsigned int i = 0; i < pvset.size(); ++i) {
    const ParmValue& pval = pvset.getParmValue(i);
    coeffIter.array() = pval.getValues();
    if (pval.hasErrors()) {
      errorIter.array() = pval.getErrors();
    }
    coeffIter.next();
    errorIter.next();
  }
  // Return the info in a Record.
  Record vals;
  vals.define("values", coeff);
  vals.define("errors", errors);
  vals.define("freqs", Vector<double>(axisx.centers()));
  vals.define("times", Vector<double>(axisy.centers()));
  vals.define("freqwidths", Vector<double>(axisx.widths()));
  vals.define("timewidths", Vector<double>(axisy.widths()));
  return vals;
}

Grid ParmFacadeLocal::record2Grid(const Record& rec) const {
  Array<double> freqs = rec.toArrayDouble("freqs");
  Array<double> freqw = rec.toArrayDouble("freqwidths");
  Array<double> times = rec.toArrayDouble("times");
  Array<double> timew = rec.toArrayDouble("timewidths");
  assert(freqs.size() == freqw.size());
  assert(times.size() == timew.size());
  return Grid(makeAxis(freqs, freqw, freqs.size()),
              makeAxis(times, timew, times.size()));
}

Axis::ShPtr ParmFacadeLocal::makeAxis(const Vector<double>& centers,
                                      const Vector<double>& widths,
                                      unsigned int n) const {
  if (centers.size() == 1) {
    return std::make_shared<RegularAxis>(centers[0] - 0.5 * widths[0],
                                         widths[0], n);
  }
  // Convert from center/width to start/end.
  vector<double> low, upp;
  (centers - 0.5 * widths).tovector(low);
  (centers + 0.5 * widths).tovector(upp);
  return Axis::makeAxis(low, upp);
}

int ParmFacadeLocal::getType(const string& str) const {
  String strc(str);
  strc.downcase();
  if (strc == "scalar") {
    return ParmValue::Scalar;
  } else if (strc == "polc") {
    return ParmValue::Polc;
  } else if (strc == "polclog") {
    return ParmValue::PolcLog;
  }
  throw std::runtime_error(strc + " is an unknown funklet type");
}

}  // namespace parmdb
}  // namespace dp3
