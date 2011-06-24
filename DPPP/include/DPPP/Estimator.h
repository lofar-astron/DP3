//# Estimator.h: Estimate complex gains introduced by the station electronics.
//#
//# Copyright (C) 2011
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

#ifndef DPPP_ESTIMATOR_H
#define DPPP_ESTIMATOR_H

// \file
// Estimate complex gains introduced by the station electronics.

#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/DPInput.h>
#include <DPPP/DPStep.h>
#include <BBSKernel/Estimate.h>
#include <BBSKernel/MeasurementExprLOFAR.h>
#include <BBSKernel/ParmManager.h>
#include <BBSKernel/VisBuffer.h>
#include <ParmDB/SourceDB.h>
#include <ParmDB/ParmDBLog.h>
#include <Common/lofar_vector.h>
#include <Common/lofar_string.h>

namespace LOFAR
{
namespace DPPP
{

// \addtogroup NDPPP
// @{

class Estimator: public DPStep
{
public:
    Estimator(DPInput *input, const ParSet &parset, const string &prefix);
    virtual ~Estimator();

    // Process the data.
    // It keeps the data.
    // When processed, it invokes the process function of the next step.
    virtual bool process(const DPBuffer &buffer);

    // Finish the processing of this step and subsequent steps.
    virtual void finish();

    // Update the general info.
    virtual void updateInfo(DPInfo &info);

    // Show the step parameters.
    virtual void show(std::ostream &out) const;

    // Show the timings.
    virtual void showTimings(std::ostream &os, double duration) const;

private:
    // For now, forbid copy construction and assignment.
    Estimator(const Estimator &other);
    Estimator &operator=(const Estimator &other);

    static const double NTIMES = 100;

    BBS::VisBuffer::Ptr convert(size_t nTimeSlots) const;
    BBS::Grid makeSolGrid(const BBS::Axis::ShPtr &timeAxis) const;

//    void copy(DPBuffer &out, size_t t);

    DPInput                         *itsInput;
    BBS::SourceDB                   *itsSourceDB;
    BBS::ParmDBLog                  *itsLog;
    BBS::MeasurementExprLOFAR::Ptr  itsModel;
    BBS::EstimateOptions            itsOptions;
    BBS::ParmGroup                  itsParms;
    BBS::BaselineMask               itsBaselineMask;
    BBS::CorrelationMask            itsCorrelationMask;
    unsigned int                    itsWindowSize;
    unsigned int                    itsCellSize;
    unsigned int                    itsCellCount;
    unsigned int                    itsTimeCount;
    casa::MDirection                itsPhaseReference;
    casa::MDirection                itsDelayReference;
    double                          itsReferenceFreq;
    BBS::Instrument::Ptr            itsInstrument;
    BBS::BaselineSeq                itsBaselines;
    BBS::CorrelationSeq             itsCorrelations;
    BBS::Axis::ShPtr                itsFreqAxis;
    vector<double>                  itsTimes;
    vector<double>                  itsTimeWidths;
    double                          itsTimeInterval;
    vector<DPBuffer>                itsDPBuffers;
    NSTimer                         itsTimer;
    NSTimer                         itsEvalTimer;
    NSTimer                         itsFlushTimer;
};

// @}

} //# namespace DPPP
} //# namespace LOFAR

#endif
