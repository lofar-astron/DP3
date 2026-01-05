// Filter.cc: DPPP step to filter out baselines and channels
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Filter.h"

#include <cassert>
#include <cstddef>
#include <tuple>

#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/casa/Containers/Record.h>

#include <xtensor/xview.hpp>

#include "base/DPBuffer.h"
#include "base/DPInfo.h"

#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"
#include "MsReader.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using casacore::Matrix;
using casacore::Record;
using casacore::RecordGram;
using casacore::TableExprNode;

namespace dp3 {
namespace steps {

Filter::Filter(const common::ParameterSet& parset, const std::string& prefix)
    : itsName(prefix),
      itsStartChanStr(parset.getString(prefix + "startchan", "0")),
      itsNrChanStr(parset.getString(prefix + "nchan", "0")),
      itsRemoveAnt(parset.getBool(prefix + "remove", false)),
      itsBaselines(parset, prefix),
      itsDoSelect(false) {}

Filter::Filter(const base::BaselineSelection& baselines)
    : itsStartChanStr("0"),
      itsNrChanStr("0"),
      itsRemoveAnt(false),
      itsBaselines(baselines),
      itsDoSelect(false) {}

Filter::~Filter() {}

common::Fields Filter::getRequiredFields() const {
  return getProvidedFields();  // A Filter requires the fields it will change.
}

common::Fields Filter::getProvidedFields() const {
  // Check if there is any setting that enables filtering. We cannot use
  // itsDoSelect here, since updateInfo() sets it and getRequiredFields() may be
  // called before updateInfo().
  common::Fields fields;
  if (itsStartChanStr != "0" || itsNrChanStr != "0" || itsRemoveAnt ||
      itsBaselines.hasSelection()) {
    fields |= kDataField | kFlagsField | kWeightsField;

    if (itsRemoveAnt || itsBaselines.hasSelection()) {
      fields |= kUvwField;
    }
  }
  return fields;
}

void Filter::updateInfo(const base::DPInfo& infoIn) {
  Step::updateInfo(infoIn);

  unsigned int nrChan;
  std::tie(itsStartChan, nrChan) = MsReader::ParseChannelSelection(
      itsStartChanStr, itsNrChanStr, getInfoOut().nchan());

  itsDoSelect = itsStartChan > 0 || nrChan < getInfoOut().nchan();

  // Handle possible baseline selection.
  if (itsBaselines.hasSelection()) {
    Matrix<bool> selbl(itsBaselines.apply(infoIn));
    const std::vector<int>& ant1 = getInfoOut().getAnt1();
    const std::vector<int>& ant2 = getInfoOut().getAnt2();
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

  // Update the DPInfo object.
  GetWritableInfoOut().SelectChannels(itsStartChan, nrChan);
  if (!itsSelBL.empty()) GetWritableInfoOut().SelectBaselines(itsSelBL);
  if (itsRemoveAnt) GetWritableInfoOut().RemoveUnusedAntennas();
}

void Filter::show(std::ostream& os) const {
  os << "Filter " << itsName << '\n';
  os << "  startchan:      " << itsStartChan << "  (" << itsStartChanStr << ')'
     << '\n';
  os << "  nchan:          " << getInfoOut().nchan() << "  (" << itsNrChanStr
     << ')' << '\n';
  itsBaselines.show(os);
  os << "  remove:         " << itsRemoveAnt << '\n';
}

void Filter::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " Filter " << itsName << '\n';
}

bool Filter::process(std::unique_ptr<DPBuffer> buffer) {
  itsTimer.start();
  if (!itsDoSelect) {
    itsTimer.stop();
    getNextStep()->process(std::move(buffer));
    return true;
  }

  // Create the new buffer and reshape it based on the sizes set in updateInfo
  std::unique_ptr<DPBuffer> filter_buffer =
      std::make_unique<DPBuffer>(buffer->GetTime(), buffer->GetExposure());
  const std::array<std::size_t, 3> filter_shape{
      getInfoOut().nbaselines(), getInfoOut().nchan(), getInfoOut().ncorr()};
  filter_buffer->GetData().resize(filter_shape);
  filter_buffer->GetFlags().resize(filter_shape);
  filter_buffer->GetWeights().resize(filter_shape);
  filter_buffer->GetUvw().resize({getInfoOut().nbaselines(), 3});

  // Select the filtered data and copy it into the new filter_buffer
  const DPBuffer::DataType& data = buffer->GetData();
  const DPBuffer::FlagsType& flags = buffer->GetFlags();
  const DPBuffer::WeightsType& weights = buffer->GetWeights();
  const DPBuffer::UvwType& uvws = buffer->GetUvw();
  if (itsSelBL.empty()) {
    // filtering on channel not baseline; copy all data for given channels to
    // make them contiguous.
    // UVW is a straight copy with no filtering because it is not dependent on
    // channel.
    filter_buffer->GetData() =
        xt::view(data, xt::all(),
                 xt::range(itsStartChan, itsStartChan + getInfoOut().nchan()),
                 xt::all());
    filter_buffer->GetFlags() =
        xt::view(flags, xt::all(),
                 xt::range(itsStartChan, itsStartChan + getInfoOut().nchan()),
                 xt::all());
    filter_buffer->GetWeights() =
        xt::view(weights, xt::all(),
                 xt::range(itsStartChan, itsStartChan + getInfoOut().nchan()),
                 xt::all());
    filter_buffer->GetUvw() = xt::view(uvws, xt::all(), xt::all());
    filter_buffer->SetRowNumbers(buffer->GetRowNumbers());
  } else {
    // Filtering on baseline and/or channel; copy all data for selected
    // baselines and channels to make them contiguous. UVW needs to be filtered
    // as well as it is dependent on baselines.
    casacore::Vector<common::rownr_t> rowNrs;
    if (!buffer->GetRowNumbers().empty()) {
      rowNrs.resize(getInfoOut().nbaselines());
    }
    // Copy the data of the selected baselines and channels.
    std::complex<float>* toData = filter_buffer->GetData().data();
    bool* toFlag = filter_buffer->GetFlags().data();
    float* toWeight = filter_buffer->GetWeights().data();
    double* toUVW = filter_buffer->GetUvw().data();
    std::size_t off = data.shape(2) * itsStartChan;  // offset of first channel
    const std::complex<float>* frData = data.data() + off;
    const bool* frFlag = flags.data() + off;
    const float* frWeight = weights.data() + off;
    const double* frUVW = uvws.data();
    int ndfr = data.shape(2) * data.shape(1);
    int ndto =
        filter_buffer->GetData().shape(2) * filter_buffer->GetData().shape(1);
    for (std::size_t i = 0; i < itsSelBL.size(); ++i) {
      if (!buffer->GetRowNumbers().empty()) {
        rowNrs[i] = buffer->GetRowNumbers()[itsSelBL[i]];
      }
      std::copy_n(frData + itsSelBL[i] * ndfr, ndto, toData);
      toData += ndto;
      std::copy_n(frFlag + itsSelBL[i] * ndfr, ndto, toFlag);
      toFlag += ndto;
      std::copy_n(frWeight + itsSelBL[i] * ndfr, ndto, toWeight);
      toWeight += ndto;
      std::copy_n(frUVW + itsSelBL[i] * 3, 3, toUVW);
      toUVW += 3;
    }
    filter_buffer->SetRowNumbers(rowNrs);
  }
  itsTimer.stop();
  getNextStep()->process(std::move(filter_buffer));
  return true;
}

void Filter::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

void Filter::addToMS(const std::string& msName) {
  Step::addToMS(msName);
  if (!itsRemoveAnt) {
    return;
  }
  // See if and which stations have been removed.
  casacore::Table antTab(msName + "/ANTENNA", casacore::Table::Update);
  casacore::Table selTab = antTab(!antTab.col("NAME").in(
      casacore::Vector<casacore::String>(getInfoOut().antennaNames())));
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
