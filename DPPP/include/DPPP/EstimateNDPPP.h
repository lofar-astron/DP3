//# EstimateNDPPP.h: NDPPP specific variant of BBS estimation routines.
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

#ifndef DPPP_ESTIMATENDPPP_H
#define DPPP_ESTIMATENDPPP_H

// \file
// NDPPP specific variant of BBS estimation routines.

#include <DPPP/DPBuffer.h>
#include <BBSKernel/MeasurementExprLOFAR.h>
#include <BBSKernel/VisDimensions.h>
#include <BBSKernel/BaselineMask.h>
#include <BBSKernel/CorrelationMask.h>
#include <ParmDB/Grid.h>
#include <BBSKernel/Estimate.h>
#include <Common/lofar_vector.h>

namespace LOFAR
{
namespace DPPP
{

// \addtogroup NDPPP
// @{

// Estimate values for the solvable parameters of \c models that best fit the
// observed visibilities in \c buffers. Each entry in \c buffers is assumed to
// be a set of visibilities phase shifted to a particular direction (and
// optionally averaged). The model in \c models at the corresponding index is
// assumed to be the model for that direction. The mixing coefficient array
// \c coeff quantifies the influence of each direction on each other direction.
// When estimating parameters this interdependence is taken into account.
//
// \param[in]   buffers
// A \c vector of visibility buffers. The outer \c vector contains separate
// directions, the inner \c vector contains separate timeslots.
// \param[in]   models
// A \c vector of measurement expressions for the separate directions.
// \param[in]   baselines
// An ordered list of baselines for which visibilities are available in \c
// buffers.
// \param[in]   correlations
// An ordered list of correlations for which visibilities are available in \c
// buffers.
// \param[in]   baselineMask
// A mask that select the baselines to process (a subset of \c baselines).
// \param[in]   correlationMask
// A mask that select the correlations to process (a subset of \c correlations).
// \param[in]   visGrid
// The frequency x time grid associated with the visibility data in \c buffers.
// \param[in]   solGrid
// The frequency x time grid of solution cells. A separate solution will be
// computed for each cell. Should be a partition of \c visGrid.
// \param[in]   options
// Options that control various details of the estimation algorithm.
void estimate(const vector<vector<DPPP::DPBuffer> > &buffers,
    const vector<BBS::MeasurementExpr::Ptr> &models,
    const vector<casa::Array<casa::DComplex> > &coeff,
    const BBS::BaselineSeq &baselines,
    const BBS::CorrelationSeq &correlations,
    const BBS::BaselineMask &baselineMask,
    const BBS::CorrelationMask &correlationMask,
    const BBS::Grid &visGrid,
    const BBS::Grid &solGrid,
    const BBS::EstimateOptions &options = BBS::EstimateOptions());

void subtract(vector<DPPP::DPBuffer> &buffer,
    const vector<BBS::MeasurementExpr::Ptr> &models,
    const vector<casa::Array<casa::DComplex> > &coeff,
    const BBS::BaselineSeq &baselines,
    const BBS::CorrelationSeq &correlations,
    const BBS::BaselineMask &baselineMask,
    const BBS::CorrelationMask &correlationMask,
    const BBS::Grid &visGrid,
    unsigned int target,
    const vector<unsigned int> &directions);

// @}

} //# namespace BBS
} //# namespace LOFAR

#endif
