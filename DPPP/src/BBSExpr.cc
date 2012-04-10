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
#include <DPPP/EstimateNDPPP.h>

#include <BBSKernel/MeasurementAIPS.h>

#include <Common/LofarLogger.h>
#include <Common/lofar_iostream.h>
#include <Common/lofar_iomanip.h>
#include <Common/StreamUtil.h>

#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MeasConvert.h>

using namespace casa;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    BBSExpr::BBSExpr (const DPInput& input, const string& skyName,
                      const string& instrumentName, double elevationCutoff)
      : itsBaselineMask    (true),
        itsCorrelationMask (true),
        itsElevationCutoff (elevationCutoff)
    {
      // Open ParmDB and SourceDB.
      try {
        itsSourceDB = boost::shared_ptr<SourceDB>
          (new SourceDB(ParmDBMeta("casa", skyName)));
        ParmManager::instance().initCategory(SKY, itsSourceDB->getParmDB());
      } catch (Exception &e) {
        THROW(Exception, "Failed to open sky model parameter database: "
          << skyName);
      }

      try {
        ParmManager::instance().initCategory(INSTRUMENT,
                                             ParmDB(ParmDBMeta("casa",
                                                               instrumentName)));
      } catch (Exception &e) {
        THROW(Exception, "Failed to open instrument model parameter database: "
          << instrumentName);
      }

      // Create Instrument instance using the information present in DPInput.
      initInstrument(input);

      // Store reference directions.
      MDirection itsDelayRef = MDirection::Convert(input.delayCenter(),
        MDirection::J2000)();
      MDirection itsTileRef  = MDirection::Convert(input.tileBeamDir(),
        MDirection::J2000)();

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

      // Use all correlations.
      ASSERT(input.ncorr() == 4);
      itsCorrelations.append(Correlation::XX);
      itsCorrelations.append(Correlation::XY);
      itsCorrelations.append(Correlation::YX);
      itsCorrelations.append(Correlation::YY);
    }

    BBSExpr::~BBSExpr()
    {
    }

    void BBSExpr::setOptions (const SolverOptions& lsqOptions)
    {
      // Initialize parameter estimation options.
      itsOptions = EstimateOptions(EstimateOptions::COMPLEX,
                                   EstimateOptions::L2,
                                   false,
                                   0,
                                   false,
                                   ~flag_t(0),
                                   flag_t(4),
                                   lsqOptions);
    }

    void BBSExpr::addModel (const string &source, const MDirection &phaseRef)
    {
      // NB: The phase reference needs to be taken from the DPInfo object,
      // because it is not necessarily equal to the phase center of the
      // observation: A PhaseShift step may have shifted the phase center
      // somewhere else.
      MDirection phaseRefJ2000 = MDirection::Convert(phaseRef,
        MDirection::J2000)();

      // Define the model configuration.
      ModelConfig config;
      config.setSources (vector<string>(1, source));
      config.setDirectionalGain();
      config.setCache();
      if (itsElevationCutoff > 0) {
        config.setElevationCutConfig (ElevationCutConfig(itsElevationCutoff));
      }
      // TODO: Find a better way to handle the reference frequency. Here we
      // use a bogus reference frequency, which does not matter since we are
      // not using the beam model. But it's still a hack, of course.
      const double refFreq = 60.0e6;

      // Create the model expression.
      itsModels.push_back (MeasurementExpr::Ptr
        (new MeasurementExprLOFAR (*itsSourceDB, BufferMap(), config,
                                   itsInstrument, itsBaselines, refFreq,
                                   phaseRefJ2000, itsDelayRef, itsTileRef)));

      // Determine parameters to estimate and store for later reference.
      vector<string> incl(1, "DirectionalGain:*"), excl;
      ParmGroup parms = ParmManager::instance().makeSubset(incl, excl,
        itsModels.back()->parms());
      itsModelParms.push_back(parms);

      // Update list of unique parameters to estimate (models can in principle
      // share parameters).
      itsParms.insert(parms.begin(), parms.end());
    }

    void BBSExpr::estimate (vector<vector<DPBuffer> >& buffers,
                            const Grid& visGrid, const Grid& solveGrid,
                            const vector<Array<DComplex> >& factors)
    {
      // Set parameter domain.
      ParmManager::instance().setDomain(solveGrid.getBoundingBox());

      // Force computation of partial derivatives of the model with respect to
      // the parameters of interest.
      setSolvables();

      // Estimate parameter values.
      DPPP::estimate (buffers, itsModels, factors, itsBaselines,
                      itsCorrelations, itsBaselineMask, itsCorrelationMask,
                      visGrid, solveGrid, itsOptions);

      // Flush solutions to disk.
      ParmManager::instance().flush();
    }

    void BBSExpr::subtract (vector<DPBuffer>& buffer,
                            const Grid& visGrid,
                            const vector<Array<DComplex> >& factors,
                            uint target,
                            uint nSources)
    {
      // Ensure computation of partial derivatives is switched off.
      clearSolvables();

      // Subtract selected sources from the data.
      LOG_DEBUG_STR("nSources: " << nSources << ",  target: " << target);
      DPPP::subtract (buffer, itsModels, factors, itsBaselines,
                      itsCorrelations, itsBaselineMask, itsCorrelationMask,
                      visGrid, target, nSources);
    }

    void BBSExpr::clearSolvables()
    {
      for (uint i=0; i<itsModels.size(); ++i) {
        itsModels[i]->clearSolvables();
      }
    }

    void BBSExpr::setSolvables()
    {
      for (uint i=0; i<itsModels.size(); ++i) {
        itsModels[i]->setSolvables (itsModelParms[i]);
      }
    }

    void BBSExpr::initInstrument(const DPInput &input)
    {
      const size_t nStations = input.antennaNames().size();

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
      itsInstrument = Instrument::Ptr(new Instrument("LOFAR", position,
        stations.begin(), stations.end()));
    }

  } //# namespace DPPP
} //# namespace LOFAR
