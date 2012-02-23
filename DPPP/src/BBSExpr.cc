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
                      const string& instrumentName)
      : itsBaselineMask    (true),
        itsCorrelationMask (true)
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
      // Create Instrument instance using information present in DPInput.
      itsInstrument = MeasurementAIPS(input.msName()).instrument();
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
      // Define the model configuration.
      /// itsConfig.setDirectionalGain();
      itsConfig.setGain();
      BeamConfig beamConfig(BeamConfig::DEFAULT, false);
      itsConfig.setBeamConfig(beamConfig);
      itsConfig.setCache();
    }

    BBSExpr::~BBSExpr()
    {}

    void BBSExpr::addModel (const DPInput& input, const DPInfo& info,
                            const string& sourceName, double refFreq)
    {
      // Model of this source plus directional gain for this source. How?
      // Get directions and make sure they are in J2000.
      // Note that the phase center has to be taken from the DPInfo object,
      // because it can be used from a PhaseShift object.
      /// What to do if phasecenter is moving (e.g. SUN)?
      MDirection phaseReference = MDirection::Convert(info.phaseCenter(),
                                                      MDirection::J2000)();
      MDirection delayReference = MDirection::Convert(input.delayCenter(),
                                                      MDirection::J2000)();
      MDirection tileBeamDir    = MDirection::Convert(input.tileBeamDir(),
                                                      MDirection::J2000)();
      vector<string> incl, excl;
      incl.push_back ("Gain.*");
      /// incl.push_back ("DirectionalGain.*");
      try {
        itsConfig.setSources (vector<string>(1, sourceName));
        itsModels.push_back (MeasurementExpr::Ptr
          (new MeasurementExprLOFAR (*itsSourceDB, BufferMap(),
                                     itsConfig, itsInstrument,
                                     itsBaselines, refFreq,
                                     phaseReference,
                                     delayReference,
                                     tileBeamDir)));
       } catch (Exception &x) {
        THROW(Exception, "Unable to construct the BBS model expression; " +
              string(x.what()));
      }
      ParmGroup parms = ParmManager::instance().makeSubset
        (incl, excl, itsModels.back()->parms());
      itsModelParms.push_back(parms);
      itsParms.insert(parms.begin(), parms.end());
      itsSources.push_back (sourceName);
      // Set solver options.
      SolverOptions lsqOptions;
      lsqOptions.maxIter = 40;
      lsqOptions.epsValue = 1e-9;
      lsqOptions.epsDerivative = 1e-9;
      lsqOptions.colFactor = 1e-9;
      lsqOptions.lmFactor = 1.0;
      lsqOptions.balancedEq = false;
      lsqOptions.useSVD = true;
      itsOptions = EstimateOptions(EstimateOptions::COMPLEX,
                                   EstimateOptions::L2, false, 1,
                                   false, ~flag_t(0), flag_t(4),
                                   lsqOptions);
    }

    void BBSExpr::estimate (vector<vector<DPBuffer> >& buffers,
                            const Grid& visGrid, const Grid& solveGrid,
                            const vector<Array<DComplex> >& factors)
    {
      // Set parameter domain.
      ParmManager::instance().setDomain(solveGrid.getBoundingBox());
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
                            uint nsources)
    {
      clearSolvables();
      LOG_DEBUG_STR("nsources: " << nsources << ",  target: " << target);
      DPPP::subtract (buffer, itsModels, factors, itsBaselines,
                      itsCorrelations, itsBaselineMask, itsCorrelationMask,
                      visGrid, target, nsources);
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

  } //# namespace DPPP
} //# namespace LOFAR
