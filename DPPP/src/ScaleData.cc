//# ScaleData.cc: DPPP step class for freq-dependent scaling of the data
//# Copyright (C) 2013
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: ScaleData.cc 23223 2012-12-07 14:09:42Z schoenmakers $
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/ScaleData.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/ParameterValue.h>
#include <Common/StreamUtil.h>
#include <Common/LofarLogger.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Utilities/Regex.h>

#include <iostream>


using namespace casa;

namespace LOFAR {
  namespace DPPP {

    ScaleData::ScaleData (DPInput*,
                        const ParameterSet& parset,
                        const string& prefix)
      : itsName           (prefix),
        itsScaleSizeGiven (False),
        itsScaleSize      (False),
        itsStationExp     (parset.getStringVector(prefix+"stations",
                                                  vector<string>())),
        itsCoeffStr       (parset.getStringVector(prefix+"coeffs",
                                                  vector<string>()))
    {
      ASSERTSTR (itsStationExp.size() == itsCoeffStr.size(),
                 "ScaleData parameters stations and coeffs differ in size");
      // Determine if scaling for size is explicitly given.
      if (parset.isDefined (prefix+"scalesize")) {
        itsScaleSizeGiven = True;
        itsScaleSize      = parset.getBool (prefix+"scalesize");
      }
    }

    ScaleData::~ScaleData()
    {}

    void ScaleData::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setNeedWrite();
      // Find out if the observation has LBA or HBA data.
      // Add the default factors to itsCoeffStr as being valid for all stations.
      // In that way they will be used if a station matches no others.
      string antennaSet (infoIn.antennaSet());
      itsStationExp.push_back ("*");
      vector<double> defCoeff(5);
      uint nNominal = 0;
      if (antennaSet.substr(0,3) == "LBA") {
        nNominal = 48;
        itsCoeffStr.push_back ("[1.02306755e+06, -7.31426342e+04,"
                               " 2.05537660e+03, -2.61310245e+01,"
                               " 1.26031118e-01]");
      } else {
        nNominal = 24;
        itsCoeffStr.push_back ("[2.35166277e+05, -6.20793100e+03,"
                               " 6.22685124e+01, -2.78418826e-01,"
                               " 4.67920578e-04]");
      }
      // Get the frequencies.
      const Vector<double>& freqs = infoIn.chanFreqs();
      // Convert the coefficients to scale factors per freq per station regex.
      vector<vector<double> > scaleVec(itsStationExp.size());
      vector<Regex> stationRegex(itsStationExp.size());
      for (uint i=0; i<scaleVec.size(); ++i) {
        // Convert the station string to a proper Regex object.
        stationRegex[i] = Regex(Regex::fromPattern(itsStationExp[i]));
        // Convert coefficients from string to double.
	ParameterValue coeffPar(itsCoeffStr[i]);
	vector<double> coeff (coeffPar.getDoubleVector());
        ASSERTSTR (coeff.size() > 0, "A ScaleData coeffs vector is empty");
        vector<double>& scales = scaleVec[i];
        scales.reserve (freqs.size());
        // Evaluate the polynomial for each frequency giving the scale factors.
        for (uint j=0; j<freqs.size(); ++j) {
          double fact = coeff[coeff.size() -1];
          for (uint k=coeff.size()-1; k>0; --k) {
            fact *= freqs[j] / 1e6;     // use freq in MHz
            fact += coeff[k-1];
          }
          scales.push_back (fact);
        }
      }
      // If needed, find the nr of tiles/dipoles used for each station and
      // fill the size scale factors.
      uint nant = infoIn.antennaNames().size();
      vector<double> extraFactors(nant, 1.);
      if (itsScaleSize  ||  !itsScaleSizeGiven) {
        fillSizeScaleFactors (nNominal, extraFactors);
	ASSERTSTR (extraFactors.size() == nant,
		   "Maybe stations have been added before doing the scaling; "
		   "that should not be done");
      }
      // Find the scale factors for each station.
      // The first matching regex is used.
      // The nr of tiles in use gives an extra scale factor.
      itsStationFactors.reserve (nant);
      for (uint i=0; i<nant; ++i) {
        for (uint j=0; j<stationRegex.size(); ++j) {
          if (infoIn.antennaNames()[i].matches (stationRegex[j])) {
            itsStationFactors.push_back (scaleVec[j]);
            // If needed, scale with the nr of dipoles/tiles actually used.
            // Do that if explicitly told so or if default coeffs are used.
            if (itsScaleSize  ||
                (!itsScaleSizeGiven  &&  j == stationRegex.size()-1)) {
              for (uint k=0; k<itsStationFactors[i].size(); ++k) {
                itsStationFactors[i][k] *= extraFactors[i];
              }
            }
            break;
          }
        }
      }
      ASSERT (itsStationFactors.size() == nant);
      // Now calculate the factors per baseline,freq,pol.
      uint nb = infoIn.nbaselines();
      uint nf = freqs.size();
      uint nc = infoIn.ncorr();
      itsFactors.resize (nc, nf, nb);
      double* factPtr = itsFactors.data();
      for (uint i=0; i<nb; ++i) {
        const vector<double>& f1 = itsStationFactors[infoIn.getAnt1()[i]];
        const vector<double>& f2 = itsStationFactors[infoIn.getAnt2()[i]];
        for (uint j=0; j<nf; ++j) {
          double fact = sqrt(f1[j] * f2[j]);
          for (uint k=0; k<nc; ++k) {
            *factPtr++ = fact;
          }
        }
      }
    }

    void ScaleData::show (std::ostream& os) const
    {
      os << "ScaleData " << itsName << std::endl;
      os << "  stations:       " << itsStationExp << std::endl;
      os << "  coeffs:         " << itsCoeffStr << std::endl;
      os << "  scalesize       ";
      if (itsScaleSizeGiven) {
        os << itsScaleSize;
      } else if (itsCoeffStr.size() == 1) {
        os << True;
      } else {
        os << True << " for stations using default coeffs, otherwise "
           << False;
      }
      os << endl;
      os << "  Scale factors per station/frequency:" << endl;
      for (uint i=0; i<itsStationFactors.size(); ++i) {
        os << "   " << getInfo().antennaNames()[i] << ' '
           << itsStationFactors[i] << endl;
      }
    }

    void ScaleData::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " ScaleData " << itsName << endl;
    }

    bool ScaleData::process (const DPBuffer& buf)
    {
      itsTimer.start();
      // Apply the scale factors.
      DPBuffer bufNew(buf);
      const IPosition shp = itsFactors.shape();
      ASSERT (buf.getData().shape() == shp);
      // Multiply the data and factors giving a new data array.
      Array<Complex> data(shp);
      arrayContTransform (static_cast<const Array<Complex>&>(buf.getData()),
                          static_cast<const Array<double>&>(itsFactors),
                          data,
                          casa::Multiplies<Complex,double,Complex>());
      bufNew.setData (data);
      itsTimer.stop();
      getNextStep()->process (bufNew);
      return true;
    }

    void ScaleData::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }

    void ScaleData::fillSizeScaleFactors (uint nNominal, vector<double>& fact)
    {
      Table ms(getInfo().msName());
      ASSERTSTR (ms.keywordSet().isDefined ("LOFAR_ANTENNA_FIELD"),
                 "ScaleData: subtable LOFAR_ANTENNA_FIELD is missing, but "
                 "is needed unless scalesize=false is given");
      Table tab(ms.keywordSet().asTable ("LOFAR_ANTENNA_FIELD"));
      // Get nr of antennae from the table to be sure it matches the
      // contents of LOFAR_ANTENNA_FIELD. Later it is checked if it matches
      // the actual nr of antennae.
      uint nant = ms.keywordSet().asTable("ANTENNA").nrow();
      fact.resize (nant);
      for (uint i=0; i<nant; ++i) {
        fact[i] = 0;
      }
      // Count the nr of used tiles (for which ELEMENT_FLAG is false).
      // A station can have multiple fields (e.g. both ears for HBA_JOINED).
      ROScalarColumn<Int> antId (tab, "ANTENNA_ID");
      ROArrayColumn<Bool> elemFlag (tab, "ELEMENT_FLAG");
      for (uint i=0; i<tab.nrow(); ++i) {
        fact[antId(i)] += 0.5*nfalse(elemFlag(i));  // X and Y are separate
      }
      // Determine the scale factor.
      for (uint i=0; i<nant; ++i) {
        fact[i] = nNominal / fact[i];
      }
    }

  } //# end namespace
}
