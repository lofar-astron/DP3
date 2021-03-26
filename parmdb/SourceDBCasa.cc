// SourceDBCasa.cc: Class for a Casa table holding sources and their parameters
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SourceDBCasa.h"
#include "ParmMap.h"

#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/ExprNodeSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableLocker.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/casa/Containers/RecordField.h>

using casacore::ArrayColumn;
using casacore::ArrayColumnDesc;
using casacore::Block;
using casacore::FileLocker;
using casacore::Regex;
using casacore::ScalarColumn;
using casacore::ScalarColumnDesc;
using casacore::SetupNewTable;
using casacore::Table;
using casacore::TableDesc;
using casacore::TableIterator;
using casacore::TableLock;
using casacore::TableLocker;

namespace dp3 {
namespace parmdb {

SourceDBCasa::SourceDBCasa(const ParmDBMeta& pdm, bool forceNew)
    : SourceDBRep(pdm, forceNew), itsSetsFilled(false), itsRowNr(1, 0) {
  string tableName = pdm.getTableName() + "/SOURCES";
  // Create the table if needed or if it does not exist yet.
  if (forceNew || !Table::isReadable(tableName)) {
    createTables(pdm.getTableName());
  }
  // Open the main table.
  itsSourceTable = Table(tableName, TableLock::UserLocking);
  // Open the names table.
  itsPatchTable = itsSourceTable.keywordSet().asTable("PATCHES");
}

SourceDBCasa::~SourceDBCasa() {}

void SourceDBCasa::lock(bool lockForWrite) {
  // Lock all tables involved to avoid unwanted flushes.
  getParmDB().lock(lockForWrite);
  itsSourceTable.lock(lockForWrite);
  itsPatchTable.lock(lockForWrite);
}

void SourceDBCasa::unlock() {
  itsPatchTable.unlock();
  itsSourceTable.unlock();
  getParmDB().unlock();
}

void SourceDBCasa::createTables(const string& tableName) {
  // Create description of SOURCES.
  TableDesc td("Local Sky Model Sources", TableDesc::Scratch);
  td.comment() =
      casacore::String("Table containing the sources in the Local Sky Model");
  td.addColumn(ScalarColumnDesc<casacore::String>("SOURCENAME"));
  td.addColumn(ScalarColumnDesc<unsigned int>("PATCHID"));
  td.addColumn(ScalarColumnDesc<int>("SOURCETYPE"));
  td.addColumn(ScalarColumnDesc<casacore::String>("REFTYPE"));
  td.addColumn(ScalarColumnDesc<unsigned int>("SPINX_NTERMS"));
  td.addColumn(ScalarColumnDesc<bool>("LOG_SI"));
  td.addColumn(ScalarColumnDesc<double>("SPINX_REFFREQ"));
  td.addColumn(ScalarColumnDesc<bool>("USE_ROTMEAS"));
  td.addColumn(ScalarColumnDesc<double>("SHAPELET_ISCALE"));
  td.addColumn(ScalarColumnDesc<double>("SHAPELET_QSCALE"));
  td.addColumn(ScalarColumnDesc<double>("SHAPELET_USCALE"));
  td.addColumn(ScalarColumnDesc<double>("SHAPELET_VSCALE"));
  td.addColumn(ArrayColumnDesc<double>("SHAPELET_ICOEFF"));
  td.addColumn(ArrayColumnDesc<double>("SHAPELET_QCOEFF"));
  td.addColumn(ArrayColumnDesc<double>("SHAPELET_UCOEFF"));
  td.addColumn(ArrayColumnDesc<double>("SHAPELET_VCOEFF"));
  // Create description of PATCHES.
  TableDesc tdpat("Local Sky Model patches", TableDesc::Scratch);
  tdpat.comment() =
      casacore::String("Table containing the patches in the Local Sky Model");
  tdpat.addColumn(ScalarColumnDesc<casacore::String>("PATCHNAME"));
  tdpat.addColumn(ScalarColumnDesc<unsigned int>("CATEGORY"));
  tdpat.addColumn(ScalarColumnDesc<double>("APPARENT_BRIGHTNESS"));
  tdpat.addColumn(ScalarColumnDesc<double>("RA"));
  tdpat.addColumn(ScalarColumnDesc<double>("DEC"));
  // Create the tables.
  string tabNameSrc = tableName + "/SOURCES";
  SetupNewTable newtab(tabNameSrc, td, Table::New);
  SetupNewTable newpattab(tabNameSrc + "/PATCHES", tdpat, Table::New);
  Table tab(newtab);
  Table pattab(newpattab);
  // PATCHES is subtable of SOURCES.
  tab.rwKeywordSet().defineTable("PATCHES", pattab);
  // Set type info.
  tab.tableInfo().setType("LSM");
  tab.tableInfo().readmeAddLine("Sources in the Local Sky Model");
  pattab.tableInfo().setType("LSMpatches");
  pattab.tableInfo().readmeAddLine("Patches in the Local Sky Model");
  // Make SOURCES subtable of the ParmDB.
  Table ptab(tableName, Table::Update);
  TableLocker locker(ptab, FileLocker::Write);
  ptab.rwKeywordSet().defineTable("SOURCES", tab);
}

void SourceDBCasa::clearTables() {
  {
    TableLocker locker(itsSourceTable, FileLocker::Write);
    casacore::Vector<common::rownr_t> rows = itsSourceTable.rowNumbers();
    itsSourceTable.removeRow(rows);
  }
  {
    TableLocker locker(itsPatchTable, FileLocker::Write);
    casacore::Vector<common::rownr_t> rows = itsPatchTable.rowNumbers();
    itsPatchTable.removeRow(rows);
  }
}

void SourceDBCasa::checkDuplicates() {
  TableLocker lockerp(itsPatchTable, FileLocker::Read);
  Table tabp = itsPatchTable.sort(
      "PATCHNAME", casacore::Sort::Ascending,
      casacore::Sort::QuickSort + casacore::Sort::NoDuplicates);
  if (tabp.nrow() != itsPatchTable.nrow())
    throw std::runtime_error(
        "The PATCHES table has " +
        std::to_string(itsPatchTable.nrow() - tabp.nrow()) +
        " duplicate patch names");
  TableLocker lockers(itsSourceTable, FileLocker::Read);
  Table tabs = itsSourceTable.sort(
      "SOURCENAME", casacore::Sort::Ascending,
      casacore::Sort::QuickSort + casacore::Sort::NoDuplicates);
  if (tabs.nrow() != itsSourceTable.nrow())
    throw std::runtime_error(
        "The SOURCES table has " +
        std::to_string(itsSourceTable.nrow() - tabs.nrow()) +
        " duplicate source names");
}

std::vector<string> SourceDBCasa::findDuplicates(Table& table,
                                                 const string& columnName) {
  TableLocker locker(table, FileLocker::Read);
  TableIterator iter(table, columnName);
  std::vector<string> result;
  while (!iter.pastEnd()) {
    if (iter.table().nrow() > 1) {
      result.push_back(
          ScalarColumn<casacore::String>(iter.table(), columnName)(0));
    }
    ++iter;
  }
  return result;
}

std::vector<string> SourceDBCasa::findDuplicatePatches() {
  return findDuplicates(itsPatchTable, "PATCHNAME");
}

std::vector<string> SourceDBCasa::findDuplicateSources() {
  return findDuplicates(itsSourceTable, "SOURCENAME");
}

void SourceDBCasa::fillSets() {
  if (!itsSetsFilled) {
    TableLocker plocker(itsPatchTable, FileLocker::Read);
    ScalarColumn<casacore::String> patchCol(itsPatchTable, "PATCHNAME");
    itsPatchSet.clear();
    for (common::rownr_t i = 0; i < itsPatchTable.nrow(); ++i) {
      itsPatchSet.insert(patchCol(i));
    }
    TableLocker slocker(itsSourceTable, FileLocker::Read);
    ScalarColumn<casacore::String> sourceCol(itsSourceTable, "SOURCENAME");
    itsSourceSet.clear();
    for (common::rownr_t i = 0; i < itsSourceTable.nrow(); ++i) {
      itsSourceSet.insert(sourceCol(i));
    }
    itsSetsFilled = true;
  }
}

bool SourceDBCasa::patchExists(const string& patchName) {
  if (!itsSetsFilled) {
    fillSets();
  }
  return itsPatchSet.find(patchName) != itsPatchSet.end();
}

bool SourceDBCasa::sourceExists(const string& sourceName) {
  if (!itsSetsFilled) {
    fillSets();
  }
  return itsSourceSet.find(sourceName) != itsSourceSet.end();
}

unsigned int SourceDBCasa::addPatch(const string& patchName, int catType,
                                    double apparentBrightness, double ra,
                                    double dec, bool check) {
  itsPatchTable.reopenRW();
  TableLocker locker(itsPatchTable, FileLocker::Write);
  // See if already existing.
  if (check) {
    if (patchExists(patchName))
      throw std::runtime_error("Patch " + patchName + " already exists");
  }
  itsPatchSet.insert(patchName);
  common::rownr_t rownr = itsPatchTable.nrow();
  itsPatchTable.addRow();
  ScalarColumn<casacore::String> nameCol(itsPatchTable, "PATCHNAME");
  ScalarColumn<unsigned int> catCol(itsPatchTable, "CATEGORY");
  nameCol.put(rownr, patchName);
  catCol.put(rownr, catType);
  writePatch(apparentBrightness, ra, dec, rownr);
  return rownr;
}

void SourceDBCasa::updatePatch(unsigned int patchId, double apparentBrightness,
                               double ra, double dec) {
  // Note: patchId is the rownr in the PATCHES subtable.
  writePatch(apparentBrightness, ra, dec, patchId);
}

void SourceDBCasa::writePatch(double apparentBrightness, double ra, double dec,
                              unsigned int rownr) {
  ScalarColumn<double> brCol(itsPatchTable, "APPARENT_BRIGHTNESS");
  ScalarColumn<double> raCol(itsPatchTable, "RA");
  ScalarColumn<double> decCol(itsPatchTable, "DEC");
  brCol.put(rownr, apparentBrightness);
  raCol.put(rownr, ra);
  decCol.put(rownr, dec);
}

void SourceDBCasa::addSource(const SourceInfo& sourceInfo,
                             const string& patchName,
                             const ParmMap& defaultParameters, double ra,
                             double dec, bool check) {
  common::rownr_t patchId;
  {
    // Find the patch.
    TableLocker locker(itsPatchTable, FileLocker::Read);
    Table table = itsPatchTable(itsPatchTable.col("PATCHNAME") ==
                                casacore::String(patchName));
    if (table.nrow() != 1)
      throw std::runtime_error("Patch " + patchName + " does not exist");
    patchId = table.rowNumbers()[0];
  }
  itsSourceTable.reopenRW();
  TableLocker locker(itsSourceTable, FileLocker::Write);
  if (check) {
    if (sourceExists(sourceInfo.getName()))
      throw std::runtime_error("Source " + sourceInfo.getName() +
                               " already exists");
  }
  itsSourceSet.insert(sourceInfo.getName());
  addSrc(sourceInfo, patchId, defaultParameters, ra, dec);
}

void SourceDBCasa::addSource(const SourceData& source, bool check) {
  ParmMap parms;
  source.getParms(parms);
  addSource(source.getInfo(), source.getPatchName(), parms, 0, 0, check);
}

void SourceDBCasa::addSource(const SourceInfo& sourceInfo,
                             const string& patchName, int catType,
                             double apparentBrightness,
                             const ParmMap& defaultParameters, double ra,
                             double dec, bool check) {
  itsPatchTable.reopenRW();
  itsSourceTable.reopenRW();
  TableLocker lockerp(itsPatchTable, FileLocker::Write);
  TableLocker lockers(itsSourceTable, FileLocker::Write);
  if (check) {
    if (patchExists(patchName))
      throw std::runtime_error("Patch " + patchName + " already exists");
    if (sourceExists(sourceInfo.getName()))
      throw std::runtime_error("Source " + sourceInfo.getName() +
                               " already exists");
  }
  itsPatchSet.insert(patchName);
  itsSourceSet.insert(sourceInfo.getName());
  common::rownr_t patchId =
      addPatch(patchName, catType, apparentBrightness, ra, dec, false);
  addSrc(sourceInfo, patchId, defaultParameters, ra, dec);
}

void SourceDBCasa::addSrc(const SourceInfo& sourceInfo, unsigned int patchId,
                          const ParmMap& defaultParameters, double ra,
                          double dec) {
  // Okay, add it to the source table.
  ScalarColumn<casacore::String> nameCol(itsSourceTable, "SOURCENAME");
  ScalarColumn<unsigned int> idCol(itsSourceTable, "PATCHID");
  ScalarColumn<int> typeCol(itsSourceTable, "SOURCETYPE");
  ScalarColumn<casacore::String> reftCol(itsSourceTable, "REFTYPE");
  ScalarColumn<unsigned int> spinxCol(itsSourceTable, "SPINX_NTERMS");
  ScalarColumn<bool> logSICol(itsSourceTable, "LOG_SI");
  ScalarColumn<double> sirefCol(itsSourceTable, "SPINX_REFFREQ");
  ScalarColumn<bool> usermCol(itsSourceTable, "USE_ROTMEAS");
  ScalarColumn<double> iscalCol(itsSourceTable, "SHAPELET_ISCALE");
  ScalarColumn<double> qscalCol(itsSourceTable, "SHAPELET_QSCALE");
  ScalarColumn<double> uscalCol(itsSourceTable, "SHAPELET_USCALE");
  ScalarColumn<double> vscalCol(itsSourceTable, "SHAPELET_VSCALE");
  ArrayColumn<double> icoefCol(itsSourceTable, "SHAPELET_ICOEFF");
  ArrayColumn<double> qcoefCol(itsSourceTable, "SHAPELET_QCOEFF");
  ArrayColumn<double> ucoefCol(itsSourceTable, "SHAPELET_UCOEFF");
  ArrayColumn<double> vcoefCol(itsSourceTable, "SHAPELET_VCOEFF");
  common::rownr_t rownr = itsSourceTable.nrow();
  itsSourceTable.addRow();
  nameCol.put(rownr, sourceInfo.getName());
  idCol.put(rownr, patchId);
  typeCol.put(rownr, sourceInfo.getType());
  reftCol.put(rownr, sourceInfo.getRefType());
  spinxCol.put(rownr, sourceInfo.getNSpectralTerms());
  logSICol.put(rownr, sourceInfo.getHasLogarithmicSI());
  sirefCol.put(rownr, sourceInfo.getSpectralTermsRefFreq());
  usermCol.put(rownr, sourceInfo.getUseRotationMeasure());
  iscalCol.put(rownr, sourceInfo.getShapeletScaleI());
  qscalCol.put(rownr, sourceInfo.getShapeletScaleQ());
  uscalCol.put(rownr, sourceInfo.getShapeletScaleU());
  vscalCol.put(rownr, sourceInfo.getShapeletScaleV());
  if (sourceInfo.getType() == SourceInfo::SHAPELET) {
    if (sourceInfo.getShapeletCoeffI().empty())
      throw std::runtime_error("No coefficients defined for shapelet source " +
                               sourceInfo.getName());
    icoefCol.put(rownr, sourceInfo.getShapeletCoeffI());
    qcoefCol.put(rownr, sourceInfo.getShapeletCoeffQ());
    ucoefCol.put(rownr, sourceInfo.getShapeletCoeffU());
    vcoefCol.put(rownr, sourceInfo.getShapeletCoeffV());
  }
  // Now add the default parameters to the ParmDB DEFAULTVALUES table.
  bool foundRa = false;
  bool foundDec = false;
  string suffix = ':' + sourceInfo.getName();
  for (ParmMap::const_iterator iter = defaultParameters.begin();
       iter != defaultParameters.end(); ++iter) {
    // Add source name suffix if not part of name yet.
    string name = iter->first;
    if (name.size() <= suffix.size() ||
        name.substr(name.size() - suffix.size()) != suffix) {
      name = name + suffix;
    }
    getParmDB().putDefValue(name, iter->second, false);
    if (name.substr(0, 3) == "Ra:") {
      foundRa = true;
      ra = iter->second.getFirstParmValue().getValues().data()[0];
    } else if (name.substr(0, 4) == "Dec:") {
      foundDec = true;
      dec = iter->second.getFirstParmValue().getValues().data()[0];
    }
  }
  // If Ra or Dec given and not in parameters, put it.
  // Use absolute perturbations for them.
  if (!foundRa && ra != -1e9) {
    ParmValue pval(ra);
    getParmDB().putDefValue("Ra:" + sourceInfo.getName(),
                            ParmValueSet(pval, ParmValue::Scalar, 1e-6, false),
                            false);
  }
  if (!foundDec && dec != -1e9) {
    ParmValue pval(dec);
    getParmDB().putDefValue("Dec:" + sourceInfo.getName(),
                            ParmValueSet(pval, ParmValue::Scalar, 1e-6, false),
                            false);
  }
}

void SourceDBCasa::deleteSources(const string& sourceNamePattern) {
  Table table = itsSourceTable;
  table.reopenRW();
  TableLocker locker(table, FileLocker::Write);
  // Get the selection from the patch table.
  Regex regex(Regex::fromPattern(sourceNamePattern));
  table = table(table.col("SOURCENAME") == regex);
  // Delete all rows found.
  itsSourceTable.removeRow(table.rowNumbers());
  // A patch will never be removed from the PATCH table, otherwise the
  // PATCHID keys (which are row numbers) do not match anymore.
  // Delete the sources from the ParmDB tables.
  string parmNamePattern = "*:" + sourceNamePattern;
  getParmDB().deleteDefValues(parmNamePattern);
  getParmDB().deleteValues(parmNamePattern,
                           Box(Point(-1e30, -1e30), Point(1e30, 1e30)));
}

Table SourceDBCasa::selectPatches(int category, const string& pattern,
                                  double minBrightness,
                                  double maxBrightness) const {
  Table table = itsPatchTable;
  if (category >= 0) {
    table = table(table.col("CATEGORY") == category);
  }
  if (!pattern.empty() && pattern != "*") {
    Regex regex(Regex::fromPattern(pattern));
    table = table(table.col("PATCHNAME") == regex);
  }
  if (minBrightness >= 0) {
    table = table(table.col("APPARENT_BRIGHTNESS") >= minBrightness);
  }
  if (maxBrightness >= 0) {
    table = table(table.col("APPARENT_BRIGHTNESS") <= maxBrightness);
  }
  return table;
}

std::vector<string> SourceDBCasa::getPatches(int category,
                                             const string& pattern,
                                             double minBrightness,
                                             double maxBrightness) {
  TableLocker locker(itsPatchTable, FileLocker::Read);
  Table table(selectPatches(category, pattern, minBrightness, maxBrightness));
  Block<casacore::String> keys(3);
  Block<int> orders(3);
  keys[0] = "CATEGORY";
  keys[1] = "APPARENT_BRIGHTNESS";
  keys[2] = "PATCHNAME";
  orders[0] = casacore::Sort::Ascending;
  orders[1] = casacore::Sort::Descending;
  orders[2] = casacore::Sort::Ascending;
  table = table.sort(keys, orders);
  casacore::Vector<casacore::String> nm(
      ScalarColumn<casacore::String>(table, "PATCHNAME").getColumn());
  return std::vector<string>(nm.cbegin(), nm.cend());
}

std::vector<PatchInfo> SourceDBCasa::getPatchInfo(int category,
                                                  const string& pattern,
                                                  double minBrightness,
                                                  double maxBrightness) {
  TableLocker locker(itsPatchTable, FileLocker::Read);
  Table table(selectPatches(category, pattern, minBrightness, maxBrightness));
  casacore::Vector<casacore::String> nm(
      ScalarColumn<casacore::String>(table, "PATCHNAME").getColumn());
  casacore::Vector<double> ra(ScalarColumn<double>(table, "RA").getColumn());
  casacore::Vector<double> dc(ScalarColumn<double>(table, "DEC").getColumn());
  casacore::Vector<unsigned int> ca(
      ScalarColumn<unsigned int>(table, "CATEGORY").getColumn());
  casacore::Vector<double> ab(
      ScalarColumn<double>(table, "APPARENT_BRIGHTNESS").getColumn());
  std::vector<PatchInfo> vec;
  vec.reserve(nm.size());
  for (size_t i = 0; i < nm.size(); ++i) {
    vec.push_back(PatchInfo(nm[i], ra[i], dc[i], ca[i], ab[i]));
  }
  return vec;
}

std::vector<SourceInfo> SourceDBCasa::getPatchSources(const string& patchName) {
  TableLocker lockerp(itsPatchTable, FileLocker::Read);
  TableLocker lockers(itsSourceTable, FileLocker::Read);
  Table table = itsPatchTable(itsPatchTable.col("PATCHNAME") ==
                              casacore::String(patchName));
  if (table.nrow() == 0) {
    return std::vector<SourceInfo>();
  }
  if (table.nrow() != 1)
    throw std::runtime_error("Patch name " + patchName +
                             " multiply defined in " +
                             std::string(itsPatchTable.tableName()));
  common::rownr_t patchid = table.rowNumbers()[0];
  table = itsSourceTable(itsSourceTable.col("PATCHID") == patchid);
  return readSources(table);
}

std::vector<SourceData> SourceDBCasa::getPatchSourceData(
    const string& patchName) {
  std::vector<SourceInfo> info = getPatchSources(patchName);
  std::vector<SourceData> result;
  result.reserve(info.size());
  for (std::vector<SourceInfo>::const_iterator iter = info.begin();
       iter != info.end(); ++iter) {
    ParmMap pmap;
    getParmDB().getDefValues(pmap, "*:" + iter->getName());
    SourceData srcData(*iter, patchName, 0, 0);
    srcData.setParms(pmap);
    result.push_back(srcData);
  }
  return result;
}

SourceInfo SourceDBCasa::getSource(const string& sourceName) {
  TableLocker lockers(itsSourceTable, FileLocker::Read);
  Table table = itsSourceTable(itsSourceTable.col("SOURCENAME") ==
                               casacore::String(sourceName));
  if (table.nrow() == 0)
    throw std::runtime_error("Source name " + sourceName + " not found in " +
                             std::string(itsSourceTable.tableName()));
  if (table.nrow() != 1)
    throw std::runtime_error("Source name " + sourceName +
                             " multiply defined in " +
                             std::string(itsSourceTable.tableName()));
  return readSources(table)[0];
}

std::vector<SourceInfo> SourceDBCasa::getSources(const string& pattern) {
  TableLocker locker(itsSourceTable, FileLocker::Read);
  // Get the selection from the patch table.
  Regex regex(Regex::fromPattern(pattern));
  Table table = itsSourceTable(itsSourceTable.col("SOURCENAME") == regex);
  return readSources(table);
}

std::vector<SourceInfo> SourceDBCasa::readSources(const Table& table) {
  casacore::Vector<casacore::String> nm(
      ScalarColumn<casacore::String>(table, "SOURCENAME").getColumn());
  casacore::Vector<int> tp(ScalarColumn<int>(table, "SOURCETYPE").getColumn());
  // Default RefType is J2000 (for backward compatibility).
  casacore::Vector<casacore::String> rt(tp.size(), "J2000");
  if (table.tableDesc().isColumn("REFTYPE")) {
    ScalarColumn<casacore::String>(table, "REFTYPE").getColumn(rt);
  }
  std::vector<SourceInfo> res;
  res.reserve(nm.size());
  if (table.tableDesc().isColumn("SPINX_NTERMS")) {
    casacore::Vector<unsigned int> sd(
        ScalarColumn<unsigned int>(table, "SPINX_NTERMS").getColumn());
    casacore::Vector<double> sr(
        ScalarColumn<double>(table, "SPINX_REFFREQ").getColumn());
    casacore::Vector<bool> rm(
        ScalarColumn<bool>(table, "USE_ROTMEAS").getColumn());
    ScalarColumn<double> iscalCol(table, "SHAPELET_ISCALE");
    ScalarColumn<double> qscalCol(table, "SHAPELET_QSCALE");
    ScalarColumn<double> uscalCol(table, "SHAPELET_USCALE");
    ScalarColumn<double> vscalCol(table, "SHAPELET_VSCALE");
    ArrayColumn<double> icoefCol(table, "SHAPELET_ICOEFF");
    ArrayColumn<double> qcoefCol(table, "SHAPELET_QCOEFF");
    ArrayColumn<double> ucoefCol(table, "SHAPELET_UCOEFF");
    ArrayColumn<double> vcoefCol(table, "SHAPELET_VCOEFF");
    ScalarColumn<bool> logSICol;
    bool hasLogSICol = table.tableDesc().isColumn("LOG_SI");
    if (hasLogSICol) {
      logSICol = ScalarColumn<bool>(table, "LOG_SI");
    }
    for (common::rownr_t i = 0; i < nm.size(); ++i) {
      SourceInfo::Type type = SourceInfo::Type((tp[i]));
      bool useLogSI = true;
      if (hasLogSICol) {
        useLogSI = logSICol(i);
      }
      res.push_back(
          SourceInfo(nm[i], type, rt[i], useLogSI, sd[i], sr[i], rm[i]));
      if (type == SourceInfo::SHAPELET) {
        if (!icoefCol.isDefined(i))
          throw std::runtime_error(
              "No coefficients defined for "
              " shapelet source " +
              nm[i]);
        res[i].setShapeletScale(iscalCol(i), qscalCol(i), uscalCol(i),
                                vscalCol(i));
        res[i].setShapeletCoeff(icoefCol(i), qcoefCol(i), ucoefCol(i),
                                vcoefCol(i));
      }
    }
  } else {
    // Columns SPINX_NTERMS, SPINX_REFFREQ, and USE_ROTMEAS were added later,
    // so be backward compatible.
    // In this case get degree and reffreq from associated parmdb.
    for (common::rownr_t i = 0; i < nm.size(); ++i) {
      ParmMap parmd;
      int degree = -1;
      double refFreq = 0.;
      getParmDB().getDefValues(parmd, nm[i] + ":SpectralIndexDegree");
      if (parmd.size() == 1) {
        degree = *(parmd.begin()->second.getDefParmValue().getValues().data());
        ParmMap parmf;
        getParmDB().getDefValues(parmf, nm[i] + ":SpectralIndexDegree");
        if (parmf.size() == 1) {
          refFreq =
              *(parmf.begin()->second.getDefParmValue().getValues().data());
        }
      }
      res.push_back(
          SourceInfo(nm[i], SourceInfo::Type(tp[i]), rt[i], true, degree + 1,
                     refFreq));  // In old format, useLogSI = true always
    }
  }
  return res;
}

bool SourceDBCasa::atEnd() { return itsRowNr[0] >= itsSourceTable.nrow(); }

void SourceDBCasa::rewind() { itsRowNr[0] = 0; }

void SourceDBCasa::getNextSource(SourceData& src) {
  TableLocker slocker(itsSourceTable, FileLocker::Read);
  TableLocker plocker(itsPatchTable, FileLocker::Read);
  // Read the main info for the current row.
  src.setInfo(readSources(itsSourceTable(itsRowNr))[0]);
  ScalarColumn<casacore::String> patchNameCol(itsPatchTable, "PATCHNAME");
  ScalarColumn<unsigned int> patchIdCol(itsSourceTable, "PATCHID");
  src.setPatchName(patchNameCol(patchIdCol(itsRowNr[0])));
  // Read the other SourceData info from the default Parm values.
  const string& srcName = src.getInfo().getName();
  // Fetch position.
  src.setRa(getDefaultParmValue("Ra:" + srcName));
  src.setDec(getDefaultParmValue("Dec:" + srcName));
  src.setI(getDefaultParmValue("I:" + srcName));
  src.setV(getDefaultParmValue("V:" + srcName));
  src.setQ(getDefaultParmValue("Q:" + srcName));
  src.setU(getDefaultParmValue("U:" + srcName));
  if (src.getInfo().getType() == SourceInfo::GAUSSIAN) {
    src.setOrientation(getDefaultParmValue("Orientation:" + srcName));
    src.setMajorAxis(getDefaultParmValue("MajorAxis:" + srcName));
    src.setMinorAxis(getDefaultParmValue("MinorAxis:" + srcName));
  } else {
    src.setOrientation(0);
    src.setMajorAxis(0);
    src.setMinorAxis(0);
  }
  // Fetch spectral index attributes (if applicable).
  size_t nTerms = src.getInfo().getNSpectralTerms();
  std::vector<double> terms;
  if (nTerms > 0) {
    terms.reserve(nTerms);
    for (size_t i = 0; i < nTerms; ++i) {
      std::ostringstream oss;
      oss << "SpectralIndex:" << i << ":" << srcName;
      terms.push_back(getDefaultParmValue(oss.str()));
    }
  }
  src.setSpectralTerms(terms);
  // Fetch rotation measure attributes (if applicable).
  if (src.getInfo().getUseRotationMeasure()) {
    src.setPolarizedFraction(
        getDefaultParmValue("PolarizedFraction:" + srcName));
    src.setPolarizationAngle(
        getDefaultParmValue("PolarizationAngle:" + srcName));
    src.setRotationMeasure(getDefaultParmValue("RotationMeasure:" + srcName));
  } else {
    src.setPolarizedFraction(0);
    src.setPolarizationAngle(0);
    src.setRotationMeasure(0);
  }
  itsRowNr[0]++;
}

double SourceDBCasa::getDefaultParmValue(const string& name) {
  ParmValueSet valueSet = getParmDB().getDefValue(name, ParmValue());
  assert(valueSet.empty() && valueSet.getType() == ParmValue::Scalar);
  const casacore::Array<double>& values =
      valueSet.getDefParmValue().getValues();
  assert(values.size() == 1);
  return values.data()[0];
}

}  // namespace parmdb
}  // namespace dp3
