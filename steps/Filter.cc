// Filter.cc: DPPP step to filter out baselines and channels
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Filter.h"

#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"
#include "../base/DPLogger.h"
#include "../base/Exceptions.h"

#include "../common/ParameterSet.h"

#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/casa/Containers/Record.h>

#include <cassert>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using casacore::IPosition;
using casacore::Matrix;
using casacore::Record;
using casacore::RecordGram;
using casacore::TableExprNode;

namespace dp3 {
namespace steps {

Filter::Filter(InputStep* input, const common::ParameterSet& parset,
               const string& prefix)
    : itsInput(input),
      itsName(prefix),
      itsStartChanStr(parset.getString(prefix + "startchan", "0")),
      itsNrChanStr(parset.getString(prefix + "nchan", "0")),
      itsRemoveAnt(parset.getBool(prefix + "remove", false)),
      itsBaselines(parset, prefix),
      itsDoSelect(false) {}

Filter::Filter(InputStep* input, const base::BaselineSelection& baselines)
    : itsInput(input),
      itsStartChanStr("0"),
      itsNrChanStr("0"),
      itsRemoveAnt(false),
      itsBaselines(baselines),
      itsDoSelect(false) {}

Filter::~Filter() {}

void Filter::updateInfo(const base::DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
  info().setWriteFlags();
  if (itsRemoveAnt) {
    info().setMetaChanged();
  }
  // Parse the chan expressions.
  // Nr of channels can be used as 'nchan' in the expressions.
  Record rec;
  rec.define("nchan", infoIn.nchan());
  TableExprNode node1(RecordGram::parse(rec, itsStartChanStr));
  TableExprNode node2(RecordGram::parse(rec, itsNrChanStr));
  // nchan=0 means until the last channel.
  double result;
  node1.get(rec, result);
  itsStartChan = (unsigned int)(result + 0.001);
  node2.get(rec, result);
  unsigned int nrChan = (unsigned int)(result + 0.0001);
  unsigned int nAllChan = getInfo().nchan();
  if (itsStartChan >= nAllChan)
    throw Exception("startchan " + std::to_string(itsStartChan) +
                    " exceeds nr of available channels (" +
                    std::to_string(nAllChan) + ')');
  unsigned int maxNrChan = nAllChan - itsStartChan;
  if (nrChan == 0) {
    nrChan = maxNrChan;
  } else {
    nrChan = std::min(nrChan, maxNrChan);
  }
  itsDoSelect = itsStartChan > 0 || nrChan < maxNrChan;
  // Handle possible baseline selection.
  if (itsBaselines.hasSelection()) {
    Matrix<bool> selbl(itsBaselines.apply(infoIn));
    const casacore::Vector<int>& ant1 = getInfo().getAnt1();
    const casacore::Vector<int>& ant2 = getInfo().getAnt2();
    itsSelBL.reserve(ant1.size());
    for (unsigned int i = 0; i < ant1.size(); ++i) {
      if (selbl(ant1[i], ant2[i])) {
        itsSelBL.push_back(i);
      }
    }
    if (itsSelBL.size() < ant1.size()) {
      itsDoSelect = true;
    }
  }
  if (itsDoSelect || itsRemoveAnt) {
    // Update the DPInfo object.
    info().update(itsStartChan, nrChan, itsSelBL, itsRemoveAnt);
    if (itsDoSelect) {
      // Shape the arrays in the buffer.
      IPosition shape(3, infoIn.ncorr(), nrChan, getInfo().nbaselines());
      itsBuf.getData().resize(shape);
      itsBuf.getFlags().resize(shape);
      itsBuf.getWeights().resize(shape);
      if (!itsSelBL.empty()) {
        itsBuf.getUVW().resize(IPosition(2, 3, shape[2]));
      }
    }
  }
}

void Filter::show(std::ostream& os) const {
  os << "Filter " << itsName << '\n';
  os << "  startchan:      " << itsStartChan << "  (" << itsStartChanStr << ')'
     << '\n';
  os << "  nchan:          " << getInfo().nchan() << "  (" << itsNrChanStr
     << ')' << '\n';
  itsBaselines.show(os);
  os << "  remove:         " << itsRemoveAnt << '\n';
}

void Filter::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " Filter " << itsName << '\n';
}

bool Filter::process(const DPBuffer& buf) {
  itsTimer.start();
  if (!itsDoSelect) {
    itsBuf.referenceFilled(buf);
    itsTimer.stop();
    getNextStep()->process(buf);
    return true;
  }
  // Get the various data arrays.
  itsBufTmp.referenceFilled(buf);
  const casacore::Array<casacore::Complex>& data = buf.getData();
  const casacore::Array<bool>& flags = buf.getFlags();
  const casacore::Array<float>& weights =
      itsInput->fetchWeights(buf, itsBufTmp, itsTimer);
  const casacore::Array<double>& uvws =
      itsInput->fetchUVW(buf, itsBufTmp, itsTimer);
  const casacore::Array<bool>& frFlags =
      itsInput->fetchFullResFlags(buf, itsBufTmp, itsTimer);
  // Size fullResFlags if not done yet.
  int frfAvg = frFlags.shape()[0] / data.shape()[1];
  if (itsBuf.getFullResFlags().empty()) {
    IPosition frfShp = frFlags.shape();
    frfShp[0] = getInfo().nchan() * frfAvg;
    frfShp[2] = getInfo().nbaselines();
    itsBuf.getFullResFlags().resize(frfShp);
  }
  // Form the blc and trc for the channel selection.
  IPosition first(3, 0);
  IPosition last(data.shape() - 1);
  first[1] = itsStartChan;
  last[1] = itsStartChan + getInfo().nchan() - 1;
  IPosition frfFirst(3, 0);
  IPosition frfLast(frFlags.shape() - 1);
  frfFirst[0] = first[1] * frfAvg;
  frfLast[0] = (last[1] + 1) * frfAvg - 1;
  // Copy the data into the output buffer.
  if (itsSelBL.empty()) {
    // No baseline selection; copy all data for given channels to
    // make them contiguous.
    // UVW can be referenced, because not dependent on channel.
    itsBuf.getData().assign(data(first, last));
    itsBuf.getFlags().assign(flags(first, last));
    itsBuf.getWeights().assign(weights(first, last));
    itsBuf.getFullResFlags().assign(frFlags(frfFirst, frfLast));
    itsBuf.setUVW(buf.getUVW());
    itsBuf.setRowNrs(buf.getRowNrs());
  } else {
    casacore::Vector<common::rownr_t> rowNrs;
    if (!buf.getRowNrs().empty()) {
      rowNrs.resize(getInfo().nbaselines());
    }
    // Copy the data of the selected baselines and channels.
    casacore::Complex* toData = itsBuf.getData().data();
    bool* toFlag = itsBuf.getFlags().data();
    float* toWeight = itsBuf.getWeights().data();
    double* toUVW = itsBuf.getUVW().data();
    bool* toFrf = itsBuf.getFullResFlags().data();
    size_t off = data.shape()[0] * first[1];  // offset of first channel
    const casacore::Complex* frData = data.data() + off;
    const bool* frFlag = flags.data() + off;
    const float* frWeight = weights.data() + off;
    const double* frUVW = uvws.data();
    int ndfr = data.shape()[0] * data.shape()[1];
    int ndto = itsBuf.getData().shape()[0] * itsBuf.getData().shape()[1];
    int nffr = frFlags.shape()[0];
    int nfto = itsBuf.getFullResFlags().shape()[0];
    for (size_t i = 0; i < itsSelBL.size(); ++i) {
      if (!buf.getRowNrs().empty()) {
        rowNrs[i] = buf.getRowNrs()[itsSelBL[i]];
      }
      casacore::objcopy(toData, frData + itsSelBL[i] * ndfr, ndto);
      toData += ndto;
      casacore::objcopy(toFlag, frFlag + itsSelBL[i] * ndfr, ndto);
      toFlag += ndto;
      casacore::objcopy(toWeight, frWeight + itsSelBL[i] * ndfr, ndto);
      toWeight += ndto;
      casacore::objcopy(toUVW, frUVW + itsSelBL[i] * 3, 3);
      toUVW += 3;
      // Copy FullResFlags for all times.
      const bool* frFrf = (frFlags.data() + frfFirst[0] +
                           itsSelBL[i] * nffr * frFlags.shape()[1]);
      for (size_t j = 0; j <= size_t(frfLast[1]); ++j) {
        casacore::objcopy(toFrf, frFrf, nfto);
        toFrf += nfto;
        frFrf += nffr;
      }
    }
    itsBuf.setRowNrs(rowNrs);
  }
  itsBuf.setTime(buf.getTime());
  itsBuf.setExposure(buf.getExposure());
  itsTimer.stop();
  getNextStep()->process(itsBuf);
  return true;
}

void Filter::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

void Filter::addToMS(const string& msName) {
  getPrevStep()->addToMS(msName);
  if (!itsRemoveAnt) {
    return;
  }
  // See if and which stations have been removed.
  casacore::Table antTab(msName + "/ANTENNA", casacore::Table::Update);
  casacore::Table selTab =
      antTab(!antTab.col("NAME").in(info().antennaNames()));
  if (selTab.nrow() == 0) {
    return;
  }
  // Remove these rows from the ANTENNA table.
  // Note that stations of baselines that have been filtered out before,
  // will also be removed.
  casacore::Vector<common::rownr_t> removedAnt = selTab.rowNumbers();
  casacore::Vector<int> antMap = createIdMap(antTab.nrow(), removedAnt);
  antTab.removeRow(removedAnt);
  // Remove and renumber the stations in other subtables.
  casacore::Table ms(msName);
  common::rownr_t nr;
  renumberSubTable(ms, "FEED", "ANTENNA_ID", removedAnt, antMap, nr);
  renumberSubTable(ms, "POINTING", "ANTENNA_ID", removedAnt, antMap, nr);
  renumberSubTable(ms, "SYSCAL", "ANTENNA_ID", removedAnt, antMap, nr);
  renumberSubTable(ms, "QUALITY_BASELINE_STATISTIC", "ANTENNA1", removedAnt,
                   antMap, nr);
  renumberSubTable(ms, "QUALITY_BASELINE_STATISTIC", "ANTENNA2", removedAnt,
                   antMap, nr);
  // Finally remove and renumber in the beam tables.
  common::rownr_t nrAntFldId;
  casacore::Vector<common::rownr_t> remAntFldId = renumberSubTable(
      ms, "LOFAR_ANTENNA_FIELD", "ANTENNA_ID", removedAnt, antMap, nrAntFldId);
  if (!remAntFldId.empty()) {
    casacore::Vector<int> antFldIdMap = createIdMap(nrAntFldId, remAntFldId);
    renumberSubTable(ms, "LOFAR_ELEMENT_FAILURE", "ANTENNA_FIELD_ID",
                     remAntFldId, antFldIdMap, nr);
  }
}

casacore::Vector<int> Filter::createIdMap(
    common::rownr_t nrId,
    const casacore::Vector<common::rownr_t>& removedIds) const {
  // Create the mapping from old to new id.
  casacore::Vector<int> idMap(nrId);
  indgen(idMap);  // fill with 0,1,2,...
  int nrrem = 0;
  for (common::rownr_t i = 0; i < removedIds.size(); ++i) {
    idMap[removedIds[i]] = -1;
    nrrem++;
    if (i < removedIds.size() - 1) {
      for (common::rownr_t j = removedIds[i] + 1; j < removedIds[i + 1]; ++j) {
        idMap[j] -= nrrem;
      }
    }
  }
  for (common::rownr_t j = removedIds[removedIds.size() - 1] + 1;
       j < idMap.size(); ++j) {
    idMap[j] -= nrrem;
  }
  return idMap;
}

casacore::Vector<common::rownr_t> Filter::renumberSubTable(
    const casacore::Table& ms, const casacore::String& name,
    const casacore::String& colName,
    const casacore::Vector<common::rownr_t>& removedAnt,
    const casacore::Vector<int>& antMap, common::rownr_t& nrId) const {
  // Exit if no such subtable.
  if (!ms.keywordSet().isDefined(name)) {
    return casacore::Vector<common::rownr_t>();
  }
  // Remove the rows of the removed stations.
  casacore::Table subTab(ms.tableName() + '/' + name, casacore::Table::Update);
  nrId = subTab.nrow();
  casacore::Table selTab =
      subTab(subTab.col(colName).in(TableExprNode(removedAnt)));
  subTab.removeRow(selTab.rowNumbers());
  // Renumber the rest.
  casacore::ScalarColumn<int> antCol(subTab, colName);
  casacore::Vector<int> antIds = antCol.getColumn();
  for (unsigned int i = 0; i < antIds.size(); ++i) {
    int newId = antMap[antIds[i]];
    assert(newId >= 0);
    antIds[i] = newId;
  }
  antCol.putColumn(antIds);
  return selTab.rowNumbers();
}

}  // namespace steps
}  // namespace dp3
