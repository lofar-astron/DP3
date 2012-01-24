//# BBSExpr.cc: Create the expression tree for BBS
//#
//# Copyright (C) 2012
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
//# $Id$

#include <lofar_config.h>
#include <DPPP/BBSExpr.h>

#include <Common/LofarLogger.h>
#include <Common/lofar_iostream.h>
#include <Common/lofar_iomanip.h>

#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MeasConvert.h>

using namespace casa;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    BBSExpr::BBSExpr(const DPInput& input, const DPInfo& info,
                     const string& sourceName)
      : itsBaselineMask (True),
        itsRefFreq      (0),
        itsTimeInterval (input.timeInterval())
    {
      // Open ParmDB and SourceDB.
      try {
        itsSourceDB = boost::shared_ptr<SourceDB>
          (new SourceDB(ParmDBMeta("casa", "sky")));
        ParmManager::instance().initCategory(SKY, itsSourceDB->getParmDB());
      } catch (Exception &e) {
        THROW(Exception, "Failed to open sky model parameter database: "
              << "sky");
      }
      try {
        ParmManager::instance().initCategory(INSTRUMENT,
                                             ParmDB(ParmDBMeta("casa",
                                                               "instrument")));
      } catch (Exception &e) {
        THROW(Exception, "Failed to open instrument model parameter database: "
              << "instrument");
      }
      // Create Instrument instance using information present in DPInput.
      size_t nStations = input.antennaNames().size();
      vector<Station::Ptr> stations;
      stations.reserve(nStations);
      for (size_t i = 0; i < nStations; ++i) {
        // Get station name and ITRF position.
        casa::MPosition position = MPosition::Convert(input.antennaPos()[i],
                                                      MPosition::ITRF)();
        // Store station information.
        stations.push_back(Station::Ptr(new Station(input.antennaNames()(i),
                                                    position)));
      }
      MPosition position = MPosition::Convert(input.arrayPos(),
                                              MPosition::ITRF)();
      itsInstrument = Instrument::Ptr
        (new Instrument("LOFAR", position, stations.begin(), stations.end()));
      // Get directions and make sure they are in J2000.
      // Note that the phase center has to be taken from the DPInfo object,
      // because it can be used from a PhaseShift object.
      /// What to do if phasecenter is moving (e.g. SUN)?
      itsPhaseReference = MDirection::Convert(info.phaseCenter(),
                                              MDirection::J2000)();
      MDirection delayReference = MDirection::Convert(input.delayCenter(),
                                                      MDirection::J2000)();
      MDirection tileBeamDir    = MDirection::Convert(input.tileBeamDir(),
                                                      MDirection::J2000)();
      // Construct frequency axis (needs channel width information).
      double chanWidth = input.chanWidths()[0];
      ASSERT (allEQ (input.chanWidths(), chanWidth));
      itsFreqAxis = Axis::ShPtr(new RegularAxis
                                (input.chanFreqs()(0) - chanWidth*0.5,
                                 chanWidth, input.chanFreqs().size()));
      // Construct measurement expression.
      // Ignore auto-correlations.
      ASSERT(input.getAnt1().size() == input.getAnt2().size());
      for(size_t i=0; i<input.getAnt1().size(); ++i) {
        int ant1 = input.getAnt1()[i];
        int ant2 = input.getAnt2()[i];
        itsBaselines.append(baseline_t(ant1, ant2));
        if (ant1 == ant2) {
          itsBaselineMask.clear (ant1, ant2);
        }
      }
      ASSERT(input.ncorr() == 4);
      itsCorrelations.append(Correlation::XX);
      itsCorrelations.append(Correlation::XY);
      itsCorrelations.append(Correlation::YX);
      itsCorrelations.append(Correlation::YY);
      ModelConfig config;
      config.setDirectionalGain();
      config.setCache();
      config.setSources (vector<string>(1, sourceName));
      /// config.setBeam();
      /// beaminfo lezen/toevoegen aan itsInstrument
      // Model of this source plus directional gain for this source. How?
      BufferMap bufferMap;
      try {
        itsModel = MeasurementExprLOFAR::Ptr(new MeasurementExprLOFAR
                                             (*itsSourceDB, bufferMap,
                                              config, itsInstrument,
                                              itsBaselines, itsRefFreq,
                                              itsPhaseReference,
                                              delayReference,
                                              tileBeamDir));
      } catch (Exception &x) {
        THROW(Exception, "Unable to construct the BBS model expression; " +
              string(x.what()));
      }
    }

    BBSExpr::~BBSExpr()
    {}

  } //# namespace DPPP
} //# namespace LOFAR
