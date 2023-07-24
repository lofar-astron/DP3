// StationAdder.cc: DPPP step class to add stations as a superstation
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "StationAdder.h"

#include <iomanip>
#include <iostream>

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableRow.h>
#include <casacore/casa/BasicMath/Functors.h>
#include <casacore/casa/Utilities/LinearSearch.h>
#include <casacore/casa/Utilities/Regex.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "../base/DPLogger.h"
#include "../base/FlagCounter.h"

using casacore::ArrayColumn;
using casacore::MPosition;
using casacore::MVPosition;
using casacore::Regex;
using casacore::ScalarColumn;
using casacore::ScalarMeasColumn;
using casacore::Table;
using casacore::TableRow;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

StationAdder::StationAdder(const common::ParameterSet& parset,
                           const string& prefix)
    : itsName(prefix),
      itsStatRec(parset.getRecord(prefix + "stations")),
      itsMinNPoint(parset.getUint(prefix + "minpoints", 1)),
      itsMakeAutoCorr(parset.getBool(prefix + "autocorr", false)),
      itsSumAutoCorr(parset.getBool(prefix + "sumauto", true)),
      itsDoAverage(parset.getBool(prefix + "average", true)),
      itsUseWeight(parset.getBool(prefix + "useweights", true)),
      itsUVWCalc(),
      itsTimer() {}

StationAdder::~StationAdder() {}

std::vector<int> StationAdder::getMatchingStations(
    const std::vector<std::string>& antennaNames,
    const std::vector<string>& patterns) {
  casacore::Vector<bool> used(antennaNames.size(), false);
  for (std::vector<string>::const_iterator iter = patterns.begin();
       iter != patterns.end(); ++iter) {
    int n = 0;
    if (iter->size() > 1 && ((*iter)[0] == '!' || (*iter)[0] == '^')) {
      Regex regex(Regex::fromPattern(iter->substr(1)));
      for (unsigned int i = 0; i < antennaNames.size(); ++i) {
        if (casacore::String(antennaNames[i]).matches(regex)) {
          used[i] = false;
          n++;
        }
      }
    } else {
      Regex regex(Regex::fromPattern(*iter));
      for (unsigned int i = 0; i < antennaNames.size(); ++i) {
        if (casacore::String(antennaNames[i]).matches(regex)) {
          used[i] = true;
          n++;
        }
      }
    }
    if (n == 0) {
      DPLOG_WARN_STR("StationAdder: no matching stations found for pattern "
                     << *iter);
    }
  }
  std::vector<int> parts;
  parts.reserve(12);  // Usually up to 12 stations are used
  for (unsigned int i = 0; i < used.size(); ++i) {
    if (used[i]) {
      parts.push_back(i);
    }
  }
  return parts;
}

void StationAdder::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  info().setMetaChanged();
  // Check the superstation definition(s).
  // They are specified as a ParameterRecord like:
  //    stations = {new1:[s1,s2,s3], new2:[s4,s5,s6]}
  // where s1, etc. can be glob patterns.
  std::vector<std::string> antennaNames = infoIn.antennaNames();
  std::vector<double> antennaDiam = infoIn.antennaDiam();
  std::vector<MPosition> antennaPos = info().antennaPos();
  // For each existing station, give id of new superstation it is used in.
  std::vector<int> newStations(antennaNames.size());
  std::fill(newStations.begin(), newStations.end(), -1);
  std::vector<string> newNames;  // Names of new superstations
  std::vector<double> newDiam;
  std::vector<MPosition> newPoss;
  for (common::ParameterRecord::const_iterator iter = itsStatRec.begin();
       iter != itsStatRec.end(); ++iter) {
    if (std::find(antennaNames.begin(), antennaNames.end(), iter->first) !=
            antennaNames.end() ||
        std::find(newNames.begin(), newNames.end(), iter->first) !=
            newNames.end()) {
      throw std::runtime_error("StationAdder: new station name " + iter->first +
                               " already exists");
    }
    // Get the ids of the stations forming the new superstation.
    // Expand possible .. in the parameter value.
    std::vector<int> parts = getMatchingStations(
        antennaNames, iter->second.expand().getStringVector());
    if (parts.empty()) {
      continue;
    }
    MVPosition newPosition;
    // Check if the stations exist and not used for other superstations.
    // Add their ITRF positions.
    for (unsigned int i = 0; i < parts.size(); ++i) {
      int inx = parts[i];
      if (newStations[inx] >= 0)
        throw std::runtime_error("Station " + antennaNames[inx] +
                                 " is used in multiple superstations");
      newStations[inx] = newNames.size();
      newPosition +=
          MPosition::Convert(antennaPos[inx], MPosition::ITRF)().getValue();
    }
    // The superstation position is the average of its parts.
    newNames.push_back(iter->first);
    newPosition *= 1.0 / parts.size();
    newPoss.push_back(MPosition(newPosition, MPosition::ITRF));
    itsParts.push_back(casacore::Vector<int>(parts));
    // Set the diameter of the new station by determining the
    // maximum distance to the center.
    double maxdist = 0;
    for (unsigned int i = 0; i < parts.size(); ++i) {
      int inx = parts[i];
      MVPosition mvdiff =
          newPosition -
          MPosition::Convert(antennaPos[inx], MPosition::ITRF)().getValue();
      const casacore::Vector<double>& diff = mvdiff.getValue();
      double dist = std::sqrt(std::accumulate(diff.cbegin(), diff.cend(), 0.0,
                                              casacore::SumSqr<double>()));
      // Add the radius of the station used.
      maxdist = std::max(maxdist, dist + 0.5 * antennaDiam[inx]);
    }
    newDiam.push_back(2 * maxdist);
  }
  // Add the new stations to the info's vectors.
  std::vector<int> ant1(info().getAnt1());
  std::vector<int> ant2(info().getAnt2());
  unsigned int nrold = antennaNames.size();
  unsigned int nrnew = nrold + newNames.size();
  antennaNames.resize(nrnew);
  antennaDiam.resize(nrnew);
  antennaPos.reserve(nrnew);
  for (unsigned int i = 0; i < newNames.size(); ++i) {
    antennaNames[nrold + i] = newNames[i];
    antennaDiam[nrold + i] = newDiam[i];
    antennaPos.push_back(newPoss[i]);
  }
  // Now determine the new baselines.
  // A new baseline is formed if ANTENNA1 or ANTENNA2 is in the list
  // of stations, but not both.
  // However, for auto-correlations they can be the same.
  std::vector<int> newbl(nrnew);
  // Loop over the superstations.
  // Note that by making this the outer loop, the baselines between
  // superstations are also formed.
  // At the end the new baselines are added to itsAnt1 and itsAnt2.
  // itsBufRows contains for each new baseline the rownrs in the DPBuffer
  // to be added for the new baseline. If rownr<0, the conjugate has to be
  // added (1 is added to rownr, otherwise 0 is ambiguous).
  // Note that a rownr can be the rownr of a new baseline.
  for (unsigned int j = 0; j < itsParts.size(); ++j) {
    std::fill(newbl.begin(), newbl.end(), -1);
    std::vector<int> newAnt1;
    std::vector<int> newAnt2;
    // Loop through all baselines and find out if a baseline should
    // be used for a superstation.
    for (unsigned int i = 0; i < ant1.size(); ++i) {
      bool havea1 = linearSearch1(itsParts[j], ant1[i]) >= 0;
      bool havea2 = linearSearch1(itsParts[j], ant2[i]) >= 0;
      int ant = nrold + j;
      int take = 0;
      if (havea1) {
        // If both stations are in same superstation, only use them
        // if autocorrelations have to be made.
        // Taking auto or cross depends on summing mode.
        if (havea2) {
          take = itsMakeAutoCorr && itsSumAutoCorr == (ant1[i] == ant2[i]);
        } else {
          ant = ant2[i];
          take = -1;  // conjugate has to be added
        }
      } else if (havea2) {
        ant = ant1[i];
        take = 1;
      }
      if (take != 0) {
        // We have a baseline for the superstation.
        // Get its index; create it if not used before.
        int blinx = newbl[ant];
        if (blinx < 0) {
          blinx = newbl[ant] = itsBufRows.size();
          itsBufRows.push_back(std::vector<int>());
          newAnt1.push_back(ant);
          newAnt2.push_back(nrold + j);
        }
        itsBufRows[blinx].push_back(take * int(i + 1));
      }
    }
    // Copy the new baselines for this superstation to the baseline list.
    // Give a warning if nothing found.
    if (newAnt1.empty()) {
      DPLOG_WARN_STR("StationAdder: no baseline found for superstation");
    } else {
      unsigned int oldsz = ant1.size();
      ant1.resize(oldsz + newAnt1.size(), true);
      ant2.resize(oldsz + newAnt1.size(), true);
      for (unsigned int i = 0; i < newAnt1.size(); ++i) {
        ant1[oldsz + i] = newAnt1[i];
        ant2[oldsz + i] = newAnt2[i];
      }
    }
  }
  // Set the new info.
  info().setAntennas(antennaNames, antennaDiam, antennaPos, ant1, ant2);
  // Setup the UVW calculator (for new baselines).
  itsUVWCalc = std::make_unique<base::UVWCalculator>(
      infoIn.phaseCenter(), infoIn.arrayPos(), antennaPos);
}

void StationAdder::show(std::ostream& os) const {
  os << "StationAdder " << itsName << '\n';
  os << "  stations:       " << itsStatRec << '\n';
  // Show the stations used for the new stations.
  unsigned int nold = getInfo().antennaNames().size() - itsParts.size();
  for (unsigned int i = 0; i < itsParts.size(); ++i) {
    os << "      " << getInfo().antennaNames()[nold + i] << ": [";
    for (unsigned int j = 0; j < itsParts[i].size(); ++j) {
      if (j > 0) os << ", ";
      os << getInfo().antennaNames()[itsParts[i][j]];
    }
    os << ']' << '\n';
  }
  os << "  minpoints:      " << itsMinNPoint << '\n';
  os << "  autocorr:       " << std::boolalpha << itsMakeAutoCorr << '\n';
  os << "  sumauto:        " << std::boolalpha << itsSumAutoCorr << '\n';
  os << "  average:        " << std::boolalpha << itsDoAverage << '\n';
  os << "  useweights:     " << std::boolalpha << itsUseWeight << '\n';
}

void StationAdder::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " StationAdder " << itsName << '\n';
}

bool StationAdder::process(std::unique_ptr<base::DPBuffer> buffer) {
  itsTimer.start();

  // Resize of data buffers can be destructive, therefore we must:
  // 1. Store the data locally.
  // TODO(AST-1330): Avoid copying the data, e.g., using 'swap'.
  const xt::xtensor<std::complex<float>, 3> input_data = buffer->GetData();
  const xt::xtensor<bool, 3> input_flags = buffer->GetFlags();
  const xt::xtensor<float, 3> input_weights = buffer->GetWeights();
  const xt::xtensor<double, 2> input_uvws = buffer->GetUvw();
  // 2. Resize the buffer to the new sizes that were set in our info object
  const std::array<std::size_t, 3> new_shape{
      getInfo().nbaselines(), getInfo().nchan(), getInfo().ncorr()};
  buffer->GetData().resize(new_shape);
  buffer->GetFlags().resize(new_shape);
  buffer->GetWeights().resize(new_shape);
  buffer->GetUvw().resize({getInfo().nbaselines(), 3});
  // 3. Copy the data back into the resized buffer; only the existing baselines
  // at the start will be filled the additional new baselines at the end will be
  // empty.
  std::copy_n(input_data.data(), input_data.size(), buffer->GetData().data());
  std::copy_n(input_flags.data(), input_flags.size(),
              buffer->GetFlags().data());
  std::copy_n(input_weights.data(), input_weights.size(),
              buffer->GetWeights().data());
  std::copy_n(input_uvws.data(), input_uvws.size(), buffer->GetUvw().data());

  // Now calculate the data pointers of the new baselines.
  unsigned int nrOldBL = input_data.shape()[0];
  unsigned int nrcc = input_data.shape()[2] * input_data.shape()[1];
  std::complex<float>* dataPtr = buffer->GetData().data() + input_data.size();
  bool* flagPtr = buffer->GetFlags().data() + input_data.size();
  float* wghtPtr = buffer->GetWeights().data() + input_data.size();
  double* uvwPtr = buffer->GetUvw().data() + input_uvws.size();
  std::vector<unsigned int> npoints(nrcc);
  std::vector<std::complex<float>> dataFlg(nrcc);
  std::vector<float> wghtFlg(nrcc);
  // Loop over all new baselines.
  for (unsigned int i = 0; i < itsBufRows.size(); ++i) {
    // Clear the data for the new baseline.
    for (unsigned int k = 0; k < nrcc; ++k) {
      dataPtr[k] = std::complex<float>();
      wghtPtr[k] = 0.;
      npoints[k] = 0;
      dataFlg[k] = std::complex<float>();
      wghtFlg[k] = 0.;
    }

    for (unsigned int k = 0; k < 3; ++k) {
      uvwPtr[k] = 0.;
    }
    double uvwWghtSum = 0.;

    // Sum the baselines forming the new baselines.
    for (unsigned int j = 0; j < itsBufRows[i].size(); ++j) {
      // Get the baseline number to use.
      // A negative one means using the conjugate.
      int blnr = itsBufRows[i][j];
      bool useConj = false;
      if (blnr < 0) {
        blnr = -blnr;
        useConj = true;
      }
      blnr--;  // decrement because blnr+1 is stored in itsBufRows
      // Get pointers to the input baseline data.
      const std::complex<float>* inDataPtr = &buffer->GetData()(blnr, 0, 0);
      const bool* inFlagPtr = &buffer->GetFlags()(blnr, 0, 0);
      const float* inWghtPtr = &buffer->GetWeights()(blnr, 0, 0);
      const double* inUvwPtr = &buffer->GetUvw()(blnr, 0);

      // Add the data, uvw, and weights if not flagged.
      // Write 4 loops to avoid having to test inside the loop.
      // Count the flagged points separately, so it can be used
      // if too many points are flagged.
      if (useConj) {
        if (itsUseWeight) {
          for (unsigned int k = 0; k < nrcc; ++k) {
            if (inFlagPtr[k]) {
              dataFlg[k] += conj(inDataPtr[k]) * inWghtPtr[k];
              wghtFlg[k] += inWghtPtr[k];
            } else {
              npoints[k]++;
              dataPtr[k] += conj(inDataPtr[k]) * inWghtPtr[k];
              wghtPtr[k] += inWghtPtr[k];
              for (int ui = 0; ui < 3; ++ui) {
                uvwPtr[ui] -= inUvwPtr[ui] * inWghtPtr[k];
              }
              uvwWghtSum += inWghtPtr[k];
            }
          }
        } else {
          for (unsigned int k = 0; k < nrcc; ++k) {
            if (inFlagPtr[k]) {
              dataFlg[k] += conj(inDataPtr[k]);
              wghtFlg[k] += 1.;
            } else {
              npoints[k]++;
              dataPtr[k] += conj(inDataPtr[k]);
              wghtPtr[k] += 1.;
              for (int ui = 0; ui < 3; ++ui) {
                uvwPtr[ui] -= inUvwPtr[ui];
              }
              uvwWghtSum += 1;
            }
          }
        }
      } else {
        if (itsUseWeight) {
          for (unsigned int k = 0; k < nrcc; ++k) {
            if (inFlagPtr[k]) {
              dataFlg[k] += inDataPtr[k] * inWghtPtr[k];
              wghtFlg[k] += inWghtPtr[k];
            } else {
              npoints[k]++;
              dataPtr[k] += inDataPtr[k] * inWghtPtr[k];
              wghtPtr[k] += inWghtPtr[k];
              for (int ui = 0; ui < 3; ++ui) {
                uvwPtr[ui] += inUvwPtr[ui] * inWghtPtr[k];
              }
              uvwWghtSum += inWghtPtr[k];
            }
          }
        } else {
          for (unsigned int k = 0; k < nrcc; ++k) {
            if (inFlagPtr[k]) {
              dataFlg[k] += inDataPtr[k];
              wghtFlg[k] += 1.;
            } else {
              npoints[k]++;
              dataPtr[k] += inDataPtr[k];
              wghtPtr[k] += 1.;
              for (int ui = 0; ui < 3; ++ui) {
                uvwPtr[ui] += inUvwPtr[ui];
              }
              uvwWghtSum += 1;
            }
          }
        }
      }
    }
    // Set the resulting flags. Average if needed.
    // Set flag if too few unflagged data points; use flagged data too.
    for (unsigned int k = 0; k < nrcc; ++k) {
      if (wghtPtr[k] == 0 || npoints[k] < itsMinNPoint) {
        flagPtr[k] = true;
        dataPtr[k] += dataFlg[k];
        wghtPtr[k] += wghtFlg[k];
      } else {
        flagPtr[k] = false;
      }
      if (itsDoAverage) {
        dataPtr[k] /= wghtPtr[k];
      }
    }

    // Average or calculate the UVW coordinate of the new station.
    if (itsDoAverage && uvwWghtSum != 0) {
      for (int ui = 0; ui < 3; ++ui) {
        uvwPtr[ui] /= uvwWghtSum;
      }
    } else {
      unsigned int blnr = nrOldBL + i;
      const std::array<double, 3> uvws =
          itsUVWCalc->getUVW(getInfo().getAnt1()[blnr],
                             getInfo().getAnt2()[blnr], buffer->getTime());
      uvwPtr[0] = uvws[0];
      uvwPtr[1] = uvws[1];
      uvwPtr[2] = uvws[2];
    }

    dataPtr += nrcc;
    flagPtr += nrcc;
    wghtPtr += nrcc;
    uvwPtr += 3;
  }
  itsTimer.stop();
  getNextStep()->process(std::move(buffer));
  return true;
}

void StationAdder::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

void StationAdder::addToMS(const string& msName) {
  Step::addToMS(msName);
  // Add the new stations to the ANTENNA subtable.
  Table antTab(msName + "/ANTENNA", Table::Update);
  ScalarColumn<casacore::String> nameCol(antTab, "NAME");
  ScalarColumn<casacore::String> typeCol(antTab, "TYPE");
  ScalarColumn<casacore::String> mountCol(antTab, "MOUNT");
  ArrayColumn<double> offsetCol(antTab, "OFFSET");
  ScalarColumn<double> diamCol(antTab, "DISH_DIAMETER");
  ScalarColumn<bool> flagCol(antTab, "FLAG_ROW");
  ScalarColumn<casacore::String> statCol;
  ScalarColumn<int> stidCol;
  ScalarMeasColumn<MPosition> phrefCol;
  ScalarMeasColumn<MPosition> posCol(antTab, "POSITION");
  // LOFAR columns are optional.
  if (antTab.tableDesc().isColumn("STATION")) {
    statCol.attach(antTab, "STATION");
  }
  if (antTab.tableDesc().isColumn("LOFAR_STATION_ID")) {
    stidCol.attach(antTab, "LOFAR_STATION_ID");
  }
  if (antTab.tableDesc().isColumn("LOFAR_PHASE_REFERENCE")) {
    phrefCol.attach(antTab, "LOFAR_PHASE_REFERENCE");
  }
  unsigned int origNant = antTab.nrow();
  // Take common info from the first row.
  casacore::String type, mount, stat;
  if (origNant > 0) {
    type = typeCol(0);
    mount = mountCol(0);
    if (!statCol.isNull()) {
      stat = statCol(0);
    }
  }
  casacore::Vector<double> offset(3, 0.0);
  // Put the data for each new antenna.
  for (unsigned int i = origNant; i < getInfo().antennaNames().size(); ++i) {
    antTab.addRow();
    nameCol.put(i, getInfo().antennaNames()[i]);
    typeCol.put(i, type);
    mountCol.put(i, mount);
    offsetCol.put(i, offset);
    diamCol.put(i, getInfo().antennaDiam()[i]);
    flagCol.put(i, false);
    posCol.put(i, getInfo().antennaPos()[i]);
    if (!statCol.isNull()) {
      statCol.put(i, stat);
    }
    if (!stidCol.isNull()) {
      stidCol.put(i, -1);
    }
    if (!phrefCol.isNull()) {
      phrefCol.put(i, getInfo().antennaPos()[i]);
    }
  }
  // For each new station, add a row to the FEED subtable.
  // It is a copy of the first row, except for the station id.
  Table feedTab(msName + "/FEED", Table::Update);
  TableRow feedRow(feedTab);
  ScalarColumn<int> antidCol(feedTab, "ANTENNA_ID");
  for (unsigned int i = origNant; i < getInfo().antennaNames().size(); ++i) {
    size_t rownr = feedTab.nrow();
    feedTab.addRow();
    feedRow.put(rownr, feedRow.get(0));
    antidCol.put(rownr, i);
  }
  // Now update the BeamInfo tables if they are present.
  Table ms(msName);
  if (ms.keywordSet().isDefined("LOFAR_ANTENNA_FIELD")) {
    updateBeamInfo(msName, origNant, antTab);
  }
}

void StationAdder::updateBeamInfo(const string& msName, unsigned int origNant,
                                  Table& antTab) {
  // Update the LOFAR_ANTENNA_FIELD table.
  // Copy all rows where ANTENNA_ID matches one of the stations used in
  // a superstation.
  Table afTab(msName + "/LOFAR_ANTENNA_FIELD", Table::Update);
  TableRow afRow(afTab);
  ScalarColumn<int> antIdCol(afTab, "ANTENNA_ID");
  for (unsigned int i = 0; i < afTab.nrow(); ++i) {
    for (unsigned int j = 0; j < itsParts.size(); ++j) {
      int inx = linearSearch1(itsParts[j], antIdCol(i));
      if (inx >= 0) {
        int rownr = afTab.nrow();
        afTab.addRow();
        afRow.put(rownr, afRow.get(i));
        antIdCol.put(rownr, j + origNant);
      }
    }
  }
  // Flush the table.
  afTab.flush();
  // Update the LOFAR_STATION table.
  Table statTab(msName + "/LOFAR_STATION", Table::Update);
  ScalarColumn<casacore::String> nameCol(statTab, "NAME");
  ScalarColumn<int> clockCol(statTab, "CLOCK_ID");
  ScalarColumn<bool> flagCol(statTab, "FLAG_ROW");
  // LOFAR_STATION_ID in the ANTENNA table needs to be updated.
  ScalarColumn<int> stidCol(antTab, "LOFAR_STATION_ID");
  // Get station-ids for each antenna-id.
  // Get clock-ids for each station-id.
  casacore::Vector<int> statIds(stidCol.getColumn());
  casacore::Vector<int> clockIds(clockCol.getColumn());
  int nextClockId = max(clockIds);
  // Loop over all new antennae.
  for (unsigned int i = 0; i < itsParts.size(); ++i) {
    if (itsParts[i].empty()) {
      break;
    }
    // Do all antennae of a new antenna share the same clock?
    // If so, use that clock-id, otherwise make a new one.
    int cid = clockIds[statIds[itsParts[i][0]]];
    for (unsigned int j = 1; j < itsParts[i].size(); ++j) {
      if (clockIds[statIds[itsParts[i][j]]] != cid) {
        cid = ++nextClockId;
        break;
      }
    }
    // Update LOFAR_STATION_ID for this new antenna.
    int rownr = statTab.nrow();
    stidCol.put(origNant + i, rownr);
    // Add new station.
    statTab.addRow();
    nameCol.put(rownr, getInfo().antennaNames()[origNant + i]);
    clockCol.put(rownr, cid);
    flagCol.put(rownr, false);
  }
}

}  // namespace steps
}  // namespace dp3
