//# BBSExpr.h: Create the expression tree for BBS
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

#ifndef DPPP_BBSEXPR_H
#define DPPP_BBSEXPR_H

// \file
// Predict visibilities given a source model

#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/DPInput.h>

#include <BBSKernel/MeasurementExprLOFAR.h>
#include <BBSKernel/CorrelationMask.h>
#include <BBSKernel/ParmManager.h>
#include <BBSKernel/VisBuffer.h>
#include <BBSKernel/Estimate.h>
#include <ParmDB/SourceDB.h>
#include <ParmDB/Grid.h>

#include <Common/lofar_vector.h>
#include <Common/lofar_string.h>

namespace LOFAR {
  namespace DPPP {

    // \addtogroup NDPPP
    // @{

    class BBSExpr
    {
    public:
      // Construct the expression for the given sky and instrument
      // parameter databases.
      BBSExpr (const DPInput& input, const string& skyName,
               const string& instrumentName);

      ~BBSExpr();

      // Set the options.
      void setOptions (const BBS::SolverOptions&);

      // Add a model expression for the given source and phase reference.
      void addModel (const string &source, const casa::MDirection &phaseRef);

      // Estimate the model parameters.
      void estimate (vector<vector<DPBuffer> >& buffers,
                     const BBS::Grid& visGrid, const BBS::Grid& solveGrid,
                     const vector<casa::Array<casa::DComplex> >& factors);

      // Subtract the sources.
      void subtract (vector<DPBuffer>& buffer, const BBS::Grid& visGrid,
                     const vector<casa::Array<casa::DComplex> >& factors,
                     uint target, uint nSources);

    private:
      // For now, forbid copy construction and assignment.
      BBSExpr(const BBSExpr& other);
      BBSExpr& operator= (const BBSExpr& other);

      // Clear the solvables in the model expressions.
      void clearSolvables();

      // Set the solvables in the model expressions to the gains.
      void setSolvables();

      // Gather information about the instrument from the input meta-data.
      void initInstrument(const DPInput &input);

      //# Data members
      boost::shared_ptr<BBS::SourceDB>  itsSourceDB;
      vector<BBS::MeasurementExpr::Ptr> itsModels;
      BBS::Instrument::Ptr              itsInstrument;
      casa::MDirection                  itsDelayRef;
      casa::MDirection                  itsTileRef;
      BBS::BaselineSeq                  itsBaselines;
      BBS::CorrelationSeq               itsCorrelations;
      BBS::BaselineMask                 itsBaselineMask;
      BBS::CorrelationMask              itsCorrelationMask;
      BBS::EstimateOptions              itsOptions;
      BBS::ParmGroup                    itsParms;
      vector<BBS::ParmGroup>            itsModelParms;
    };

// @}

  } //# namespace DPPP
} //# namespace LOFAR

#endif
