//# Estimator.cc: Estimate complex gains introduced by the station electronics.
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

#include <lofar_config.h>
#include <DPPP/Estimator.h>
#include <DPPP/ParSet.h>
//#include <BBSKernel/MeasurementAIPS.h>
#include <Common/ParameterSet.h>
#include <Common/LofarLogger.h>
#include <Common/lofar_iostream.h>
#include <Common/lofar_iomanip.h>

#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MeasConvert.h>

using namespace casa;
using namespace LOFAR::BBS;

namespace LOFAR
{
namespace DPPP
{

Estimator::Estimator(DPInput *input, const ParSet &parset, const string &prefix)
    :   itsInput(input),
        itsSourceDB(0),
        itsLog(0),
        itsTimeCount(0),
        itsReferenceFreq(0.0)
{
    ASSERT(input);

    itsWindowSize = parset.getUint(prefix + "windowsize", 100);
    itsCellSize = parset.getUint(prefix + "cellsize", 10);
    itsCellCount = parset.getUint(prefix + "cellcount", 1);

    // Open ParmDB, SourceDB, and solver log. The latter should be made
    // optional.
    try
    {
        itsSourceDB = new SourceDB(ParmDBMeta("casa", "sky"));
        ParmManager::instance().initCategory(SKY, itsSourceDB->getParmDB());
    }
    catch(Exception &e)
    {
        THROW(Exception, "Failed to open sky model parameter database: "
            << "sky");
    }

    try
    {
        ParmManager::instance().initCategory(INSTRUMENT,
            ParmDB(ParmDBMeta("casa", "instrument")));
    }
    catch(Exception &e)
    {
        THROW(Exception, "Failed to open instrument model parameter database: "
            << "instrument");
    }

    casa::Path path("solver.log");
    itsLog = new ParmDBLog(path.absoluteName(), ParmDBLoglevel::NONE);

    // Create Instrument instance using information present in DPInput.
    size_t nStations = input->antennaNames().size();

    vector<Station::Ptr> stations;
    stations.reserve(nStations);
    for(size_t i = 0; i < nStations; ++i)
    {
        // Get station name and ITRF position.
        casa::MPosition position = MPosition::Convert(input->antennaPos()[i],
            MPosition::ITRF)();

        // Store station information.
        stations.push_back(Station::Ptr(new Station(input->antennaNames()(i),
            position)));
    }

    MPosition position = MPosition::Convert(input->arrayPos(),
        MPosition::ITRF)();
    itsInstrument = Instrument::Ptr(new Instrument("LOFAR", position,
        stations.begin(), stations.end()));

//    MeasurementAIPS __tmp("TEST.MS");
//    itsInstrument = __tmp.instrument();

    // Get phase reference from DPInput.
    itsPhaseReference = MDirection::Convert(input->phaseCenter(),
        MDirection::J2000)();

    // Construct frequency axis (need channel width information).
    itsFreqAxis = Axis::ShPtr(new RegularAxis(input->chanFreqs()(0) - 760.0
        / 2.0, 760.0, input->chanFreqs().size()));

    // Construct measurement expression (assume solving for direction
    // independent gain).
    ASSERT(input->getAnt1().size() == input->getAnt2().size());
    for(size_t i = 0; i < input->getAnt1().size(); ++i)
    {
        itsBaselines.append(baseline_t(input->getAnt1()(i),
            input->getAnt2()(i)));
    }

    ASSERT(input->ncorr() == 4);
    itsCorrelations.append(Correlation::XX);
    itsCorrelations.append(Correlation::XY);
    itsCorrelations.append(Correlation::YX);
    itsCorrelations.append(Correlation::YY);

    ModelConfig config;
    config.setGain();
    config.setCache();

    try
    {
        itsModel =
            MeasurementExprLOFAR::Ptr(new MeasurementExprLOFAR(*itsSourceDB,
            config, itsInstrument, itsBaselines, itsPhaseReference,
            itsReferenceFreq));
    }
    catch(Exception &e)
    {
        THROW(Exception, "Unable to construct the model expression.");
    }

    SolverOptions lsqOptions;
//    lsqOptions.maxIter = 200;
//    lsqOptions.epsValue = 1e-8;
//    lsqOptions.epsDerivative = 1e-8;
//    lsqOptions.colFactor = 1e-6;
//    lsqOptions.lmFactor = 1e-3;
//    lsqOptions.balancedEq = false;
//    lsqOptions.useSVD = true;

    lsqOptions.maxIter = 200;
    lsqOptions.epsValue = 1e-9;
    lsqOptions.epsDerivative = 1e-9;
    lsqOptions.colFactor = 1e-9;
    lsqOptions.lmFactor = 1.0;
    lsqOptions.balancedEq = false;
    lsqOptions.useSVD = true;

    itsOptions = EstimateOptions(EstimateOptions::COMPLEX, EstimateOptions::L2,
        false, itsCellCount, false, ~flag_t(0), flag_t(4), lsqOptions);
//        options.setThreshold(command.rmsThreshold().begin(),
//          command.rmsThreshold().end());
//        options.setEpsilon(command.epsilon().begin(), command.epsilon().end());

    vector<string> incl, excl;
    incl.push_back("Gain:{0:0,1:1}:*");
    itsParms = ParmManager::instance().makeSubset(incl, excl,
        itsModel->parms());

    // Deselect auto-correlations.
    itsBaselineMask = BaselineMask(true);
    for(size_t i = 0; i < nStations; ++i)
    {
        itsBaselineMask.clear(i, i);
    }

    itsCorrelationMask = CorrelationMask(true);

    itsTimes.resize(itsWindowSize);
    itsDPBuffers.resize(itsWindowSize);
}

Estimator::~Estimator()
{
    delete itsSourceDB;
    delete itsLog;
}

bool Estimator::process(const DPBuffer &buffer)
{
    itsTimer.start();

    itsTimes[itsTimeCount] = buffer.getTime();
    itsDPBuffers[itsTimeCount] = buffer;

    // Make sure the weights are available (why is there no
    // DPInfo::setNeedWeights() that takes care of this?).
    DPBuffer &myBuffer = itsDPBuffers[itsTimeCount];
    myBuffer.getWeights() = itsInput->fetchWeights(myBuffer,
        myBuffer.getRowNrs());

    ++itsTimeCount;

    if(itsTimeCount == itsWindowSize)
    {
        itsTimeCount = 0;

        // Copy and convert input data.
        VisBuffer::Ptr buffer = convert(itsWindowSize);

        // Construct solution grid.
        Grid grid = makeSolGrid(buffer->grid()[TIME]);

        // Set parameter domain.
        ParmManager::instance().setDomain(grid.getBoundingBox());

        // Estimate gains.
        itsEvalTimer.start();
        estimate(*itsLog, buffer, itsBaselineMask, itsCorrelationMask, itsModel,
            grid, itsParms, itsOptions);
        itsEvalTimer.stop();

        // Flush solutions to disk.
        itsFlushTimer.start();
        ParmManager::instance().flush();
        itsFlushTimer.stop();

        // Flush all time slots down the pipeline.
        itsTimer.stop();
        for(size_t i = 0; i < itsWindowSize; ++i)
        {
            getNextStep()->process(itsDPBuffers[i]);
        }
        itsTimer.start();
    }

    itsTimer.stop();

    return true;
}

void Estimator::finish()
{
    if(itsTimeCount > 0)
    {
        itsTimer.start();

        // Copy and convert input data.
        VisBuffer::Ptr buffer = convert(itsTimeCount);

        // Construct solution grid.
        Grid grid = makeSolGrid(buffer->grid()[TIME]);

        // Set parameter domain.
        ParmManager::instance().setDomain(grid.getBoundingBox());

        // Estimate gains.
        itsEvalTimer.start();
        estimate(*itsLog, buffer, itsBaselineMask, itsCorrelationMask, itsModel,
            grid, itsParms, itsOptions);
        itsEvalTimer.stop();

        // Flush solutions to disk.
        itsFlushTimer.start();
        ParmManager::instance().flush();
        itsFlushTimer.stop();

        itsTimer.stop();
    }

    // Flush all remaining time slots down the pipeline.
    for(size_t i = 0; i < itsTimeCount; ++i)
    {
        getNextStep()->process(itsDPBuffers[i]);
    }

    getNextStep()->finish();
}

void Estimator::updateInfo(DPInfo &info)
{
    info.setNeedVisData();
    itsTimeInterval = info.timeInterval();
    itsTimeWidths.resize(itsWindowSize, itsTimeInterval);

    // Test with multiple timeslots...
//    info.update(1, itsWindowSize);
}

void Estimator::show(std::ostream &os) const
{
    os << "Estimator" << endl;
    os << "    windowsize: " << itsWindowSize << endl;
    os << "    cellsize:   " << itsCellSize << endl;
    os << "    cellcount:  " << itsCellCount << endl;
}

void Estimator::showTimings(std::ostream &os, double duration) const
{
    os << "  ";
    FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
    os << " Estimator" << endl;

    os << "           ";
    FlagCounter::showPerc1(os, itsEvalTimer.getElapsed(),
        itsTimer.getElapsed());
    os << " Expression evaluation" << endl;

    os << "           ";
    FlagCounter::showPerc1(os, itsFlushTimer.getElapsed(),
        itsTimer.getElapsed());
    os << " Flush" << endl;

    os << "           ";
    FlagCounter::showPerc1(os, itsTimer.getElapsed() - itsEvalTimer.getElapsed()
        - itsFlushTimer.getElapsed(), itsTimer.getElapsed());
    os << " Other" << endl;
}

Grid Estimator::makeSolGrid(const Axis::ShPtr &timeAxis) const
{
    Axis::ShPtr freqAxis(new RegularAxis(itsFreqAxis->start(),
        itsFreqAxis->end(), 1, true));
    return Grid(freqAxis, timeAxis->compress(itsCellSize));
}

VisBuffer::Ptr Estimator::convert(size_t nTimeSlots) const
{
    ASSERT(nTimeSlots <= itsWindowSize);

    // Create time axis.
    Axis::ShPtr timeAxis;
    if(nTimeSlots == itsWindowSize)
    {
        timeAxis = Axis::ShPtr(new OrderedAxis(itsTimes, itsTimeWidths));
    }
    else
    {
        vector<double> times(itsTimes.begin(), itsTimes.begin() + nTimeSlots);
        vector<double> widths(itsTimeWidths.begin(), itsTimeWidths.begin()
            + nTimeSlots);
        timeAxis = Axis::ShPtr(new OrderedAxis(times, widths));
    }

    // Allocate buffer.
    VisDimensions dims;
    dims.setBaselines(itsBaselines);
    dims.setCorrelations(itsCorrelations);
    dims.setGrid(Grid(itsFreqAxis, timeAxis));

    VisBuffer::Ptr buffer = VisBuffer::Ptr(new VisBuffer(dims));
    buffer->setInstrument(itsInstrument);
    buffer->setPhaseReference(itsPhaseReference);
    buffer->setReferenceFreq(itsReferenceFreq);

    // Fill buffer from the accumulated DPBuffers.
    vector<DPBuffer>::const_iterator it = itsDPBuffers.begin();
    vector<DPBuffer>::const_iterator end = itsDPBuffers.begin() + nTimeSlots;
    for(size_t t = 0; it != end; ++it, ++t)
    {
        casa::Cube<casa::Complex> samples(it->getData());
        ASSERT(samples.shape().isEqual(IPosition(3, 4, itsFreqAxis->size(),
            itsBaselines.size())));

        casa::Cube<casa::Float> weights(it->getWeights());
        ASSERT(weights.shape().isEqual(IPosition(3, 4, itsFreqAxis->size(),
            itsBaselines.size())));

        casa::Cube<casa::Bool> flags(it->getFlags());
        ASSERT(flags.shape().isEqual(IPosition(3, 4, itsFreqAxis->size(),
            itsBaselines.size())));

        for(size_t i = 0; i < itsBaselines.size(); ++i)
        for(size_t j = 0; j < itsFreqAxis->size(); ++j)
        for(size_t k = 0; k < 4; ++k)
        {
            buffer->samples[i][t][j][k] = samples(k, j, i);
            buffer->flags[i][t][j][k] = flags(k, j, i);
            buffer->covariance[i][t][j][k][k] = 1.0 / weights(k, j, i);
//            buffer->covariance[i][t][j][k][k] = (float) (1.0 / weights(k, j, i));
        }
    }

    return buffer;
}

//void Estimator::copy(DPBuffer &out, size_t t)
//{
//    // Copy simulated visibilities to DPBuffer.
//    casa::Cube<casa::Complex> tmp(4, itsFreqAxis->size(), itsBaselines.size());
//    ASSERT(tmp.shape() == out.getData().shape());

//    for(size_t i = 0; i < itsBaselines.size(); ++i)
//    for(size_t j = 0; j < itsFreqAxis->size(); ++j)
//    for(size_t k = 0; k < 4; ++k)
//    {
//        tmp(k, j, i) = itsBuffer->samples[i][t][j][k];
//    }
//    out.setData(tmp);
//}

} //# namespace DPPP
} //# namespace LOFAR
