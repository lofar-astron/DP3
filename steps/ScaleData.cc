// ScaleData.cc: DPPP step class for freq-dependent scaling of the data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "ScaleData.h"

#include "../base/DPBuffer.h"
#include "../base/BDABuffer.h"
#include "../base/DPInfo.h"
#include "../base/Exceptions.h"

#include "../common/ParameterSet.h"
#include "../common/ParameterValue.h"
#include "../common/StreamUtil.h"

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Functors.h>
#include <casacore/casa/Utilities/Regex.h>

#include <cassert>
#include <iostream>

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

using casacore::IPosition;
using casacore::Regex;

namespace dp3 {
namespace steps {

ScaleData::ScaleData(InputStep*, const common::ParameterSet& parset,
                     const string& prefix, MsType inputType)
    : itsName(prefix),
      itsInputType(inputType),
      itsScaleSizeGiven(false),
      itsScaleSize(false),
      itsStationExp(parset.getStringVector(prefix + "stations",
                                           std::vector<std::string>())),
      itsCoeffStr(parset.getStringVector(prefix + "coeffs",
                                         std::vector<std::string>())) {
  if (itsStationExp.size() != itsCoeffStr.size())
    throw Exception("ScaleData parameters stations and coeffs differ in size");
  // Determine if scaling for size is explicitly given.
  if (parset.isDefined(prefix + "scalesize")) {
    itsScaleSizeGiven = true;
    itsScaleSize = parset.getBool(prefix + "scalesize");
  }
}

ScaleData::~ScaleData() {}

void ScaleData::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
  // Find out if the observation has LBA or HBA data.
  // Add the default factors to itsCoeffStr as being valid for all stations.
  // In that way they will be used if a station matches no others.
  std::string antennaSet(infoIn.antennaSet());
  itsStationExp.push_back("*");
  std::vector<double> defCoeff(5);
  unsigned int nNominal = 0;
  if (antennaSet.substr(0, 3) == "LBA") {
    nNominal = 48;
    itsCoeffStr.push_back(
        "[1.02306755e+06, -7.31426342e+04,"
        " 2.05537660e+03, -2.61310245e+01,"
        " 1.26031118e-01]");
  } else {
    nNominal = 24;
    itsCoeffStr.push_back(
        "[2.35166277e+05, -6.20793100e+03,"
        " 6.22685124e+01, -2.78418826e-01,"
        " 4.67920578e-04]");
  }
  // Get the frequencies.
  const casacore::Vector<double>& freqs = infoIn.chanFreqs();
  // Convert the coefficients to scale factors per freq per station regex.
  std::vector<std::vector<double>> scaleVec(itsStationExp.size());
  std::vector<Regex> stationRegex(itsStationExp.size());
  for (unsigned int i = 0; i < scaleVec.size(); ++i) {
    // Convert the station string to a proper Regex object.
    stationRegex[i] = Regex(Regex::fromPattern(itsStationExp[i]));
    // Convert coefficients from string to double.
    common::ParameterValue coeffPar(itsCoeffStr[i]);
    std::vector<double> coeff(coeffPar.getDoubleVector());
    if (coeff.size() <= 0)
      throw Exception("A ScaleData coeffs vector is empty");
    std::vector<double>& scales = scaleVec[i];
    scales.reserve(freqs.size());
    // Evaluate the polynomial for each frequency giving the scale factors.
    for (unsigned int j = 0; j < freqs.size(); ++j) {
      double fact = coeff[coeff.size() - 1];
      for (unsigned int k = coeff.size() - 1; k > 0; --k) {
        fact *= freqs[j] / 1e6;  // use freq in MHz
        fact += coeff[k - 1];
      }
      scales.push_back(fact);
    }
  }
  // If needed, find the nr of tiles/dipoles used for each station and
  // fill the size scale factors.
  unsigned int nant = infoIn.antennaNames().size();
  std::vector<double> extraFactors(nant, 1.);
  if (itsScaleSize || !itsScaleSizeGiven) {
    fillSizeScaleFactors(nNominal, extraFactors);
    if (extraFactors.size() != nant)
      throw Exception(
          "Maybe stations have been added before doing the scaling; "
          "that should not be done");
  }
  // Find the scale factors for each station.
  // The first matching regex is used.
  // The nr of tiles in use gives an extra scale factor.
  itsStationFactors.reserve(nant);
  for (unsigned int i = 0; i < nant; ++i) {
    for (unsigned int j = 0; j < stationRegex.size(); ++j) {
      if (infoIn.antennaNames()[i].matches(stationRegex[j])) {
        itsStationFactors.push_back(scaleVec[j]);
        // If needed, scale with the nr of dipoles/tiles actually used.
        // Do that if explicitly told so or if default coeffs are used.
        if (itsScaleSize ||
            (!itsScaleSizeGiven && j == stationRegex.size() - 1)) {
          for (unsigned int k = 0; k < itsStationFactors[i].size(); ++k) {
            itsStationFactors[i][k] *= extraFactors[i];
          }
        }
        break;
      }
    }
  }
  // Now calculate the factors per baseline,freq,pol.
  unsigned int nb = infoIn.nbaselines();
  unsigned int nf = freqs.size();
  unsigned int nc = infoIn.ncorr();
  itsFactors.resize(nc, nf, nb);
  double* factPtr = itsFactors.data();
  for (unsigned int i = 0; i < nb; ++i) {
    const std::vector<double>& f1 = itsStationFactors[infoIn.getAnt1()[i]];
    const std::vector<double>& f2 = itsStationFactors[infoIn.getAnt2()[i]];
    for (unsigned int j = 0; j < nf; ++j) {
      double fact = sqrt(f1[j] * f2[j]);
      for (unsigned int k = 0; k < nc; ++k) {
        *factPtr++ = fact;
      }
    }
  }
}

void ScaleData::show(std::ostream& os) const {
  os << "ScaleData " << itsName << '\n';
  os << "  stations:       " << itsStationExp << '\n';
  os << "  coeffs:         " << itsCoeffStr << '\n';
  os << "  scalesize       ";
  if (itsScaleSizeGiven) {
    os << itsScaleSize;
  } else if (itsCoeffStr.size() == 1) {
    os << true;
  } else {
    os << true << " for stations using default coeffs, otherwise " << false;
  }
  os << '\n';
  os << "  Scale factors per station/frequency:" << '\n';
  for (unsigned int i = 0; i < itsStationFactors.size(); ++i) {
    os << "   " << getInfo().antennaNames()[i] << ' ' << itsStationFactors[i]
       << '\n';
  }
}

void ScaleData::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " ScaleData " << itsName << '\n';
}

bool ScaleData::process(const DPBuffer& buf) {
  itsTimer.start();
  // Apply the scale factors.
  DPBuffer bufNew(buf);
  const IPosition shp = itsFactors.shape();
  assert(buf.getData().shape() == shp);
  // Multiply the data and factors giving a new data array.
  casacore::Array<casacore::Complex> data(shp);
  arrayContTransform(
      static_cast<const casacore::Array<casacore::Complex>&>(buf.getData()),
      static_cast<const casacore::Array<double>&>(itsFactors), data,
      casacore::Multiplies<casacore::Complex, double, casacore::Complex>());
  bufNew.setData(data);
  itsTimer.stop();
  getNextStep()->process(bufNew);
  return true;
}

bool ScaleData::process(std::unique_ptr<BDABuffer> bda_buffer) {
  itsTimer.start();

  // Apply the scale factors.
  std::vector<BDABuffer::Row> rows = bda_buffer->GetRows();
  for (std::size_t row_nr = 0; row_nr < rows.size(); ++row_nr) {
    const casacore::Array<double>& factors =
        itsFactors[rows[row_nr].baseline_nr];
    // Verify vectors are the same size
    assert(rows[row_nr].n_correlations * rows[row_nr].n_channels ==
           factors.size());

    std::complex<float>* data = bda_buffer->GetData(row_nr);
    for (const double& factor : factors) {
      *data *= factor;
      ++data;
    }
  }

  itsTimer.stop();
  getNextStep()->process(std::move(bda_buffer));
  return true;
}

void ScaleData::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

void ScaleData::fillSizeScaleFactors(unsigned int nNominal,
                                     std::vector<double>& fact) {
  casacore::Table ms(getInfo().msName());
  if (!ms.keywordSet().isDefined("LOFAR_ANTENNA_FIELD"))
    throw Exception(
        "ScaleData: subtable LOFAR_ANTENNA_FIELD is missing, but "
        "is needed unless scalesize=false is given");
  casacore::Table tab(ms.keywordSet().asTable("LOFAR_ANTENNA_FIELD"));
  // Get nr of antennae from the table to be sure it matches the
  // contents of LOFAR_ANTENNA_FIELD. Later it is checked if it matches
  // the actual nr of antennae.
  unsigned int nant = ms.keywordSet().asTable("ANTENNA").nrow();
  fact.resize(nant);
  for (unsigned int i = 0; i < nant; ++i) {
    fact[i] = 0;
  }
  // Count the nr of used tiles (for which ELEMENT_FLAG is false).
  // A station can have multiple fields (e.g. both ears for HBA_JOINED).
  casacore::ScalarColumn<int> antId(tab, "ANTENNA_ID");
  casacore::ArrayColumn<bool> elemFlag(tab, "ELEMENT_FLAG");
  for (unsigned int i = 0; i < tab.nrow(); ++i) {
    fact[antId(i)] += 0.5 * nfalse(elemFlag(i));  // X and Y are separate
  }
  // Determine the scale factor.
  for (unsigned int i = 0; i < nant; ++i) {
    fact[i] = nNominal / fact[i];
  }
}

}  // namespace steps
}  // namespace dp3
