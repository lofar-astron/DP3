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
#include <ParmDB/SourceDB.h>

#include <Common/lofar_vector.h>
#include <Common/lofar_string.h>

namespace LOFAR {
  namespace DPPP {

    // \addtogroup NDPPP
    // @{

    class BBSExpr
    {
    public:
      // Define the shared pointer for this type.
      typedef shared_ptr<BBSExpr> ShPtr;

      // Construct the expression for the given source.
      BBSExpr(const DPInput&, const DPInfo&, const string& sourceName);

      ~BBSExpr();

      // Get the resulting expression.
      BBS::MeasurementExprLOFAR::Ptr getModel()
        { return itsModel; }

      // Get the frequency axis.
      const BBS::Axis::ShPtr& getFreqAxis() const
        { return itsFreqAxis; }

      // Get the baseline mask.
      const BBS::BaselineMask& getBaselineMask() const
        { return itsBaselineMask; }

    private:
      // For now, forbid copy construction and assignment.
      BBSExpr(const BBSExpr& other);
      BBSExpr& operator= (const BBSExpr& other);

      //# Data members
      boost::shared_ptr<BBS::SourceDB> itsSourceDB;
      BBS::MeasurementExprLOFAR::Ptr   itsModel;
      BBS::ParmGroup                   itsParms;
      BBS::BaselineMask                itsBaselineMask;
      casa::MDirection                 itsPhaseReference;
      double                           itsRefFreq;
      double                           itsTimeInterval;
      BBS::Instrument::Ptr             itsInstrument;
      BBS::BaselineSeq                 itsBaselines;
      BBS::CorrelationSeq              itsCorrelations;
      BBS::Axis::ShPtr                 itsFreqAxis;
    };

// @}

  } //# namespace DPPP
} //# namespace LOFAR

#endif
