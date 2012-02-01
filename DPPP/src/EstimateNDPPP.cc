//# EstimateNDPPP.cc: NDPPP specific variant of BBS estimation routines.
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
#include <DPPP/EstimateNDPPP.h>
#include <BBSKernel/EstimateUtil.h>

namespace LOFAR
{
namespace DPPP
{
using namespace BBS;

namespace
{
    // Processing statistics and timers.
    class Statistics
    {
    public:
        enum Counter
        {
            C_ALL,
            C_FLAGGED,
            C_ZERO_WEIGHT,
            C_INVALID_RESIDUAL,
            C_INVALID_DERIVATIVE,
            C_INVALID_WEIGHT,
            C_OUTLIER,
            N_Counter
        };

        enum Timer
        {
            T_EVALUATE,
            T_MIX,
            T_MAKE_COEFF_MAP,
            T_PROCESS_CELL,
            T_MODIFIER,
            T_MAKE_NORM,
            N_Timer
        };

        Statistics();

        void inc(Counter counter);
        void inc(Counter counter, size_t count);
        void reset(Counter counter);

        void reset(Timer timer);
        void start(Timer timer);
        void stop(Timer timer);

        void reset();

        string counters() const;
        string timers() const;

    private:
        size_t          itsCounters[N_Counter];
        static string   theirCounterNames[N_Counter];

        NSTimer         itsTimers[N_Timer];
        static string   theirTimerNames[N_Timer];
    };

    // State kept for a single cell in the solution grid.
    struct Cell
    {
        // Flag that indicates if processing has completed.
        bool            done;

        // LSQ solver and current estimates for the coefficients.
        casa::LSQFit    solver;
        vector<double>  coeff;
    };

    // Cell counts separated by state. Used to indicate the status of all the
    // cells after performing an iteration.
    struct IterationStatus
    {
        unsigned int nActive, nConverged, nStopped, nNoReduction, nSingular;
    };

//    struct EstimateContext
//    {
//        vector<pair<size_t, size_t> >   baselineMap;
//        vector<pair<size_t, size_t> >   correlationMap;
//        vector<Interval<size_t> >       cellMap[2];
//        map<PValueKey, unsigned int>    coeffMap;

//        Statistics                      stats;
//        vector<Cell>                    cells;
//    };

    template <typename T_SAMPLE_MODIFIER>
    struct SampleProcessorComplex;

    JonesMatrix mix(const vector<JonesMatrix> &in,
        const vector<casa::Array<casa::DComplex> > &coeff,
        unsigned int target,
        const vector<unsigned int> &directions,
        unsigned int baseline);

    void mix(const vector<JonesMatrix> &in, vector<JonesMatrix> &out,
        const vector<casa::Array<casa::DComplex> > &coeff,
        const Location &start, const Location &end, unsigned int baseline);

    FlagArray mergeFlags(const vector<JonesMatrix> &in);
    void mixAndMerge(const Matrix &factor, const Element &in, Element &out);

    Matrix makeMixingFactor(const vector<casa::Array<casa::DComplex> > &coeff,
        const Location &start, const Location &end, unsigned int baseline,
        unsigned int correlation, unsigned int row, unsigned int column);

    void subtract2(vector<DPPP::DPBuffer> &buffer,
        const vector<MeasurementExpr::Ptr> &models,
        const vector<casa::Array<casa::DComplex> > &coeff,
        unsigned int target,
        const vector<unsigned int> &directions,
        const vector<pair<size_t, size_t> > &baselineMap,
        const vector<pair<size_t, size_t> > &correlationMap);

    template <typename T_SAMPLE_PROCESSOR>
    void equate(const Location &start, const Location &end,
        const vector<vector<DPPP::DPBuffer> > &buffers,
        const vector<MeasurementExpr::Ptr> &models,
        const vector<casa::Array<casa::DComplex> > &coeff,
        const vector<pair<size_t, size_t> > &baselineMap,
        const vector<pair<size_t, size_t> > &correlationMap,
        const vector<Interval<size_t> > (&cellMap)[2],
        const map<PValueKey, unsigned int> &coeffMap,
        vector<Cell> &cells,
        Statistics &stats);

    template <typename T_SAMPLE_PROCESSOR>
    void equate2(const Location &start, const Location &end, size_t blIndex,
        const vector<DPPP::DPBuffer> &obs, const JonesMatrix &sim,
        const vector<pair<size_t, size_t> > &correlationMap,
        const vector<Interval<size_t> > (&cellMap)[2],
        const map<PValueKey, unsigned int> &coeffMap,
        vector<Cell> &cells,
        Statistics &stats);

    void initCells(const Location &start, const Location &end,
        const ParmGroup &solvables, size_t nCoeff,
        const EstimateOptions &options,
        vector<Cell> &cells);

    // Perform a single iteration for all cells in the range [\p start, \p end)
    // that have not yet converged or failed, updating the coefficient values
    // to the new estimates found.
    // \pre The range starting at \p cell should contain exactly one Cell
    // instance for each cell in the range [\p start, \p end].
    IterationStatus iterate(const Location &start, const Location &end,
        const ParmGroup &solvables, const EstimateOptions &options,
        vector<Cell> &cells);

    void updateIterationStatus(const Cell &cell, IterationStatus &status);
} // unnamed namespace

void estimate(const vector<vector<DPPP::DPBuffer> > &buffers,
    const vector<MeasurementExpr::Ptr> &models,
    const vector<casa::Array<dcomplex> > &coeff,
    const BaselineSeq &baselines,
    const CorrelationSeq &correlations,
    const BaselineMask &baselineMask,
    const CorrelationMask &correlationMask,
    const Grid &visGrid,
    const Grid &solGrid,
    const EstimateOptions &options)
{
    // =========================================================================
    // CHECK PRECONDITIONS
    // =========================================================================
    {
        ASSERT(buffers.size() == models.size() && buffers.size() > 0);

        CorrelationSeq tmp;
        tmp.append(Correlation::XX);
        tmp.append(Correlation::XY);
        tmp.append(Correlation::YX);
        tmp.append(Correlation::YY);

        for(size_t i = 0; i < models.size(); ++i)
        {
            ASSERT(models[i]->baselines() == models.front()->baselines());
            ASSERT(models[i]->correlations() == tmp);
            ASSERT(models[i]->domain().contains(visGrid.getBoundingBox()));
        }
    }
    // =========================================================================

    const size_t nDirections = buffers.size();

    // Construct a sequence of pairs of indices of matching baselines (i.e.
    // baselines common to both buffer and model).
    vector<pair<size_t, size_t> > blMap;
    makeIndexMap(baselines, models.front()->baselines(), baselineMask,
        back_inserter(blMap));

    // Construct a sequence of pairs of indices of matching correlations
    // (i.e. correlations known by both buffer and model).
    vector<pair<size_t, size_t> > crMap;
    makeIndexMap(correlations, models.front()->correlations(), correlationMask,
        back_inserter(crMap));

    // Compute a mapping from cells of the solution grid to cell intervals
    // in the evaluation grid.
    vector<Interval<size_t> > cellMap[2];
    Interval<size_t> domain[2];
    domain[FREQ] = makeAxisMap(solGrid[FREQ], visGrid[FREQ],
        back_inserter(cellMap[FREQ]));
    domain[TIME] = makeAxisMap(solGrid[TIME], visGrid[TIME],
        back_inserter(cellMap[TIME]));

    ParmGroup solvables;
    for(size_t i = 0; i < nDirections; ++i)
    {
        ParmGroup tmp = models[i]->solvables();
        solvables.insert(tmp.begin(), tmp.end());
    }

//    LOG_INFO_STR("Selection: Baselines: " << blMap.size() << "/"
//        << baselines.size() << " Correlations: " << crMap.size() << "/"
//        << correlations.size() << " Parameters: " << solvables.size());
//        << "/" << models->nParm());

    // Make coefficient map.
    map<PValueKey, unsigned int> coeffMap;
    makeCoeffMap(solvables, inserter(coeffMap, coeffMap.begin()));
    LOG_INFO_STR("No. of coefficients to estimate: " << coeffMap.size());

    // Assign solution grid to solvables.
    ParmManager::instance().setGrid(solGrid, solvables);

//    // Log information that is valid for all chunks.
//    logCoefficientIndex(log, solvables);
//    logLSQOptions(log, options.lsqOptions());

    // ---------------------------------------------------------------------
    // Process each chunk of cells in a loop.
    // ---------------------------------------------------------------------

    // Clip chunk size to the size of the solution grid.
    size_t chunkSize = options.chunkSize() == 0 ? solGrid[TIME]->size()
        : std::min(options.chunkSize(), solGrid[TIME]->size());

    // Allocate cells.
//    boost::multi_array<Cell, 2> cells(boost::extents[chunkSize]
//        [solGrid[FREQ]->size()]);
    vector<Cell> cells(solGrid[FREQ]->size() * chunkSize);

    // Compute the number of cell chunks to process.
    size_t nChunks = (solGrid[TIME]->size() + chunkSize - 1) / chunkSize;

//    Timer::instance().reset();
    Statistics stats;

    NSTimer timer;
    timer.start();

    // Process the solution grid in chunks.
    for(size_t chunk = 0; chunk < nChunks; ++chunk)
    {
        NSTimer timerChunk, timerEquate, timerIterate;
        timerChunk.start();

        // Compute cell chunk boundaries in solution grid coordinates.
        Location chunkStart(0, chunk * chunkSize);
        Location chunkEnd(solGrid[FREQ]->size() - 1,
            std::min(chunkStart.second + chunkSize - 1,
            solGrid[TIME]->size() - 1));

        // Adjust cell chunk boundaries to exclude those cells for which no
        // visibility data is available.
        chunkStart =
            Location(std::max(chunkStart.first, domain[FREQ].start),
            std::max(chunkStart.second, domain[TIME].start));
        chunkEnd =
            Location(std::min(chunkEnd.first, domain[FREQ].end),
            std::min(chunkEnd.second, domain[TIME].end));

        // If there are no cells for which visibility data is available,
        // skip the chunk.
        if(chunkStart.first > chunkEnd.first
            || chunkStart.second > chunkEnd.second)
        {
//            LOG_DEBUG_STR("chunk: " << (chunk + 1) << "/" << nChunks
//                << " status: **skipped**");
            timerChunk.stop();
            continue;
        }

        // Ensure a model value is computed for all the visibility samples
        // within the chunk.
        Location reqStart(cellMap[FREQ][chunkStart.first].start,
            cellMap[TIME][chunkStart.second].start);
        Location reqEnd(cellMap[FREQ][chunkEnd.first].end,
            cellMap[TIME][chunkEnd.second].end);
        for(size_t i = 0; i < nDirections; ++i)
        {
            models[i]->setEvalGrid(visGrid.subset(reqStart, reqEnd));
        }

        // Initialize a cell instance for each cell in [chunkEnd,
        // chunkStart].
        initCells(chunkStart, chunkEnd, solvables, coeffMap.size(), options,
            cells);

        typedef SampleProcessorComplex<SampleModifierComplex> SampleProcessor;

//        Statistics stats;
        IterationStatus status = {0, 0, 0, 0, 0};
        unsigned int nIterations = 0;
        while(true)
        {
            // Construct normal equations from the data and an evaluation of
            // the model based on the current coefficient values.
            timerEquate.start();
            equate<SampleProcessor>(chunkStart, chunkEnd, buffers, models,
                coeff, blMap, crMap, cellMap, coeffMap, cells, stats);
            timerEquate.stop();

            // Perform a single iteration.
            timerIterate.start();
            status = iterate(chunkStart, chunkEnd, solvables, options, cells);
            timerIterate.stop();

            // Notify model that solvables have changed.
            for(size_t i = 0; i < nDirections; ++i)
            {
                models[i]->solvablesChanged();
            }

            // Update iteration count.
            ++nIterations;

            // If no active cells remain in this chunk (i.e. all cells have
            // converged or have been stopped), then move to the next chunk
            // of cells.
            if(status.nActive == 0)
            {
                break;
            }
        }
        timerChunk.stop();

        // Output statistics and timers.
        const size_t nCells = (chunkEnd.second - chunkStart.second + 1)
            * (chunkEnd.first - chunkStart.first + 1);
        LOG_DEBUG_STR("chunk: " << (chunk + 1) << "/" << nChunks
            << " cells: " << nCells << " iterations: " << nIterations
            << " status: " << status.nConverged << "/" << status.nStopped
            << "/" << status.nNoReduction << "/" << status.nSingular
            << " converged/stopped/noreduction/singular");
//        LOG_DEBUG_STR("\t" << stats.counters());
//        LOG_DEBUG_STR("\t" << stats.timers());
        LOG_DEBUG_STR("\ttimers: all: " << toString(timerChunk)
            << " equate: " << toString(timerEquate) << " iterate: "
            << toString(timerIterate) << " total/count/average");

        // Propagate solutions to the next chunk if required.
        if(options.propagate() && (chunk + 1) < nChunks)
        {
            Location srcStart(0, chunk * chunkSize);
            Location srcEnd(solGrid[FREQ]->size() - 1,
                srcStart.second + chunkSize - 1);

            Location destStart(0, (chunk + 1) * chunkSize);
            Location destEnd(solGrid[FREQ]->size() - 1,
                std::min(destStart.second + chunkSize - 1,
                solGrid[TIME]->size() - 1));

            passCoeff(solvables, srcStart, srcEnd, destStart, destEnd);
        }
    }

    timer.stop();

    ostringstream oss;
    oss << endl << "Estimate statistics:" << endl;
    {
        const double elapsed = timer.getElapsed();
        const unsigned long long count = timer.getCount();
        double average = count > 0 ? elapsed / count : 0.0;

        oss << "TIMER s ESTIMATE ALL" << " total " << elapsed << " count "
            << count << " avg " << average << endl;
    }

    LOG_DEBUG_STR("Estimate statistics:");
    LOG_DEBUG_STR("\t" << stats.counters());
    LOG_DEBUG_STR("\t" << stats.timers());

//    Timer::instance().dump(oss);
//    LOG_DEBUG(oss.str());

//    Timer::instance().reset();
}

void subtract(vector<DPPP::DPBuffer> &buffer,
    const vector<BBS::MeasurementExpr::Ptr> &models,
    const vector<casa::Array<casa::DComplex> > &coeff,
    const BBS::BaselineSeq &baselines,
    const BBS::CorrelationSeq &correlations,
    const BBS::BaselineMask &baselineMask,
    const BBS::CorrelationMask &correlationMask,
    const BBS::Grid &visGrid,
    unsigned int target,
    const vector<unsigned int> &directions)
{
    // =========================================================================
    // CHECK PRECONDITIONS
    // =========================================================================
    {
        ASSERT(target < models.size());

        ASSERT(directions.size() > 0);
        for(size_t i = 0; i < directions.size(); ++i)
        {
            ASSERT(directions[i] < models.size());
        }

        ASSERT(buffer.size() > 0);

        CorrelationSeq tmp;
        tmp.append(Correlation::XX);
        tmp.append(Correlation::XY);
        tmp.append(Correlation::YX);
        tmp.append(Correlation::YY);

        for(size_t i = 0; i < models.size(); ++i)
        {
            ASSERT(models[i]->baselines() == models.front()->baselines());
            ASSERT(models[i]->correlations() == tmp);
            ASSERT(models[i]->domain().contains(visGrid.getBoundingBox()));
        }
    }
    // =========================================================================

    // Construct a sequence of pairs of indices of matching baselines (i.e.
    // baselines common to both buffer and model).
    vector<pair<size_t, size_t> > blMap;
    makeIndexMap(baselines, models.front()->baselines(), baselineMask,
        back_inserter(blMap));

    // Construct a sequence of pairs of indices of matching correlations
    // (i.e. correlations known by both buffer and model).
    vector<pair<size_t, size_t> > crMap;
    makeIndexMap(correlations, models.front()->correlations(), correlationMask,
        back_inserter(crMap));

    for(size_t i = 0; i < models.size(); ++i)
    {
        models[i]->setEvalGrid(visGrid);
    }

    subtract2(buffer, models, coeff, target, directions, blMap, crMap);
}


namespace
{

template <typename T_SAMPLE_PROCESSOR>
void equate(const Location &start, const Location &end,
    const vector<vector<DPPP::DPBuffer> > &buffers,
    const vector<MeasurementExpr::Ptr> &models,
    const vector<casa::Array<casa::DComplex> > &coeff,
    const vector<pair<size_t, size_t> > &baselineMap,
    const vector<pair<size_t, size_t> > &correlationMap,
    const vector<Interval<size_t> > (&cellMap)[2],
    const map<PValueKey, unsigned int> &coeffMap,
    vector<Cell> &cells, Statistics &stats)

{
    ASSERT(buffers.size() == models.size());

    const size_t nDirections = buffers.size();

    vector<JonesMatrix> sim(nDirections);
    vector<JonesMatrix> mixed(nDirections);

    typedef vector<pair<size_t, size_t> >::const_iterator bl_iterator;
    for(bl_iterator it = baselineMap.begin(), it_end = baselineMap.end();
        it != it_end; ++it)
    {
        // Evaluate models.
        stats.start(Statistics::T_EVALUATE);
        for(size_t i = 0; i < nDirections; ++i)
        {
            sim[i] = models[i]->evaluate(it->second);
        }
        stats.stop(Statistics::T_EVALUATE);

        // Mix
        Location visStart(cellMap[FREQ][start.first].start,
            cellMap[TIME][start.second].start);
        Location visEnd(cellMap[FREQ][end.first].end,
            cellMap[TIME][end.second].end);
        stats.start(Statistics::T_MIX);
        mix(sim, mixed, coeff, visStart, visEnd, it->first);
        stats.stop(Statistics::T_MIX);

        // Flags will be equal for all mixed simulations. Skip baseline if all
        // grid points are flagged.
        if(sim.front().hasFlags())
        {
            FlagArray flags(sim.front().flags());
            if(flags.rank() == 0 && (*flags.begin() != 0))
            {
                continue;
            }
        }

        // Equate.
        for(size_t i = 0; i < nDirections; ++i)
        {
            equate2<T_SAMPLE_PROCESSOR>(start, end, it->first, buffers[i], mixed[i], correlationMap,
                cellMap, coeffMap, cells, stats);
        }
    } // baselines
}

template <typename T_SAMPLE_PROCESSOR>
void equate2(const Location &start, const Location &end, size_t blIndex,
    const vector<DPPP::DPBuffer> &obs, const JonesMatrix &sim,
    const vector<pair<size_t, size_t> > &correlationMap,
    const vector<Interval<size_t> > (&cellMap)[2],
    const map<PValueKey, unsigned int> &coeffMap,
    vector<Cell> &cells, Statistics &stats)
{
//    // can we somehow properly clear() if nUnkowns does not change?
//    // instead of using set() which reallocates?

    const unsigned int nFreq = cellMap[FREQ][end.first].end
        - cellMap[FREQ][start.first].start + 1;
    const unsigned int nTime = cellMap[TIME][end.second].end
        - cellMap[TIME][start.second].start + 1;

    double *reSim = 0, *imSim = 0;
    vector<unsigned int> coeffIndex(coeffMap.size());
    vector<double> reDerivative(coeffMap.size());
    vector<double> imDerivative(coeffMap.size());
    vector<double*> reSimDerivative(coeffMap.size(), 0);
    vector<double*> imSimDerivative(coeffMap.size(), 0);

    typedef vector<pair<size_t, size_t> >::const_iterator cr_iterator;
    for(cr_iterator cr_it = correlationMap.begin(), cr_end = correlationMap.end();
        cr_it != cr_end; ++cr_it)
    {
        const Element element = sim.element(cr_it->second);
        if(element.size() <= 1)
        {
            continue;
        }

        const size_t nCoeff = element.size() - 1;

        // -----------------------------------------------------------------
        // Setup pointers and strides to access the model value and
        // derivatives.
        // -----------------------------------------------------------------
        stats.start(Statistics::T_MAKE_COEFF_MAP);
        Matrix sim(element.value());
        ASSERT(sim.isComplex() && sim.isArray()
            && static_cast<unsigned int>(sim.nx()) == nFreq
            && static_cast<unsigned int>(sim.ny()) == nTime);
        sim.dcomplexStorage(reSim, imSim);

        size_t i = 0;
        for(Element::const_iterator el_it = element.begin(),
            el_end = element.end(); el_it != el_end; ++el_it, ++i)
        {
            // Look-up coefficient index for this coefficient.
            map<PValueKey, unsigned int>::const_iterator coeff_it =
                coeffMap.find(el_it->first);
            ASSERT(coeff_it != coeffMap.end());
            coeffIndex[i] = coeff_it->second;

            // Get pointers to the real and imaginary part of the partial
            // derivarive of the model with respect to this coefficient.
            Matrix derivative(el_it->second);
            ASSERT(derivative.isComplex() && derivative.isArray()
                && static_cast<unsigned int>(derivative.nx()) == nFreq
                && static_cast<unsigned int>(derivative.ny()) == nTime);
            derivative.dcomplexStorage(reSimDerivative[i], imSimDerivative[i]);
        }
        stats.stop(Statistics::T_MAKE_COEFF_MAP);

        size_t offset[2];
        offset[FREQ] = cellMap[FREQ][start.first].start;
        offset[TIME] = cellMap[TIME][start.second].start;

        vector<Cell>::iterator cell = cells.begin();
        for(CellIterator it(start, end); !it.atEnd(); ++it, ++cell)
        {
            if(cell->solver.isReady())
            {
                // Skip cell if it is inactive (converged or failed).
                continue;
            }

            stats.start(Statistics::T_PROCESS_CELL);

            const Interval<size_t> &freqInterval = cellMap[FREQ][it->first];
            const Interval<size_t> &timeInterval = cellMap[TIME][it->second];

            size_t index = (timeInterval.start - offset[TIME]) * nFreq
                + (freqInterval.start - offset[FREQ]);

            for(size_t t = timeInterval.start; t <= timeInterval.end; ++t)
            {
                const DPPP::DPBuffer &buffer = obs[t];

                for(size_t f = freqInterval.start; f <= freqInterval.end; ++f)
                {
                    const bool &flagged = buffer.getFlags()(cr_it->first, f, blIndex);
                    if(!flagged)
                    {
                        for(size_t i = 0; i < nCoeff; ++i)
                        {
                            reDerivative[i] = reSimDerivative[i][index];
                            imDerivative[i] = imSimDerivative[i][index];
                        }

                        const fcomplex &vis = buffer.getData()(cr_it->first, f, blIndex);
                        const float &weight = buffer.getWeights()(cr_it->first, f, blIndex);

                        T_SAMPLE_PROCESSOR::process(*cell, weight, real(vis),
                            imag(vis), reSim[index], imSim[index], nCoeff,
                            &(reDerivative[0]), &(imDerivative[0]),
                            &(coeffIndex[0]), stats);
                    }

                    ++index;
                }

                index -= (freqInterval.end - freqInterval.start + 1);
                index += nFreq;
            }

            // merge LSQFit object.
            stats.stop(Statistics::T_PROCESS_CELL);
        }
    }
}

    void subtract2(vector<DPPP::DPBuffer> &buffer,
        const vector<MeasurementExpr::Ptr> &models,
        const vector<casa::Array<casa::DComplex> > &coeff,
        unsigned int target,
        const vector<unsigned int> &directions,
        const vector<pair<size_t, size_t> > &baselineMap,
        const vector<pair<size_t, size_t> > &correlationMap)
    {
        vector<JonesMatrix> sim(directions.size());

        typedef vector<pair<size_t, size_t> >::const_iterator
            index_map_iterator;

        for(index_map_iterator bl_it = baselineMap.begin(),
            bl_end = baselineMap.end(); bl_it != bl_end; ++bl_it)
        {
            // Evaluate models.
            for(size_t i = 0; i < directions.size(); ++i)
            {
                sim[i] = models[directions[i]]->evaluate(bl_it->second);
            }

            // Mix
            JonesMatrix mixed = mix(sim, coeff, target, directions,
                bl_it->first);

            // Subtract.
            for(index_map_iterator cr_it = correlationMap.begin(),
                cr_end = correlationMap.end(); cr_it != cr_end; ++cr_it)
            {
                Matrix crMixed = mixed.element(cr_it->second).value();
                ASSERT(!crMixed.isNull());

                const unsigned int nFreq = crMixed.nx();
                const unsigned int nTime = crMixed.ny();
                ASSERT(crMixed.isComplex() && crMixed.isArray()
                    && nTime == buffer.size());

                const double *mixed_re = 0, *mixed_im = 0;
                crMixed.dcomplexStorage(mixed_re, mixed_im);

                for(vector<DPBuffer>::iterator buffer_it = buffer.begin(),
                    buffer_end = buffer.end(); buffer_it != buffer_end;
                    ++buffer_it)
                {
                    casa::Cube<casa::Complex> &data = buffer_it->getData();
                    ASSERT(data.shape()(1) == nFreq);

                    for(size_t i = 0; i < nFreq; ++i)
                    {
                        data(cr_it->first, i, bl_it->first) =
                            makedcomplex(*mixed_re++, *mixed_im++);
                    } // frequency
                } // time
            } // correlations
        } // baselines
    }

    JonesMatrix mix(const vector<JonesMatrix> &in,
        const vector<casa::Array<casa::DComplex> > &coeff,
        unsigned int target,
        const vector<unsigned int> &directions,
        unsigned int baseline)
    {
        const unsigned int nFreq = coeff.front().shape()(3);
        const unsigned int nTime = coeff.size();

        ASSERT(nFreq >= 1 && nTime >= 1);
        const Location start(0, 0);
        const Location end(nFreq - 1, nTime - 1);

        Matrix out[4];
        for(unsigned int i = 0; i < directions.size(); ++i)
        {
            for(unsigned int correlation = 0; correlation < 4; ++correlation)
            {
                Matrix weight = makeMixingFactor(coeff, start, end, baseline,
                    correlation, target, directions[i]);

                const Matrix sim = in[i].element(correlation).value();

                if(out[correlation].isNull())
                {
                    out[correlation] = weight * sim;
                }
                else
                {
                    out[correlation] += weight * sim;
                }
            } // correlations
        } // directions

        JonesMatrix result(out[0], out[1], out[2], out[3]);
        result.setFlags(mergeFlags(in));
        return result;
    }

    void mix(const vector<JonesMatrix> &in, vector<JonesMatrix> &out,
        const vector<casa::Array<casa::DComplex> > &coeff,
        const Location &start, const Location &end, unsigned int baseline)
    {
        // dims array: ndir x ndir x ncorr x nchan x nbl (minor -> major).
        // beter dims: nbl x ncorr x ndir x ndir x nchan ???
        // TODO: Diagonal of mixing matrix == 1, so could optimize for this.

        ASSERT(in.size() == out.size());

        FlagArray flags = mergeFlags(in);

        const unsigned int nDirections = in.size();
        for(unsigned int i = 0; i < nDirections; ++i)
        {
            Element element[4];
            for(unsigned int j = 0; j < nDirections; ++j)
            {
                for(unsigned int k = 0; k < 4; ++k)
                {
                    Matrix factor = makeMixingFactor(coeff, start, end,
                        baseline, k, i, j);

                    mixAndMerge(factor, in[j].element(k), element[k]);
                }
            }

            out[i] = JonesMatrix(element[0], element[1], element[2],
                element[3]);
            out[i].setFlags(flags);
        }
    }

    FlagArray mergeFlags(const vector<JonesMatrix> &in)
    {
        vector<JonesMatrix>::const_iterator first = in.begin();
        vector<JonesMatrix>::const_iterator last = in.end();

        for(; first != last && !first->hasFlags(); ++first)
        {
        }

        if(first == last)
        {
            return FlagArray();
        }

        FlagArray flags = first->flags().clone();
        ++first;

        for(; first != last; ++first)
        {
            if(first->hasFlags())
            {
                flags |= first->flags();
            }
        }

        return flags;
    }

    void mixAndMerge(const Matrix &factor, const Element &in, Element &out)
    {
        // Update value.
        Matrix value = out.value();
        if(value.isNull())
        {
            out.assign(factor * in.value());
        }
        else
        {
            value += factor * in.value();
        }

        // Update partial derivatives.
        Element::const_iterator inIter = in.begin();
        Element::const_iterator inEnd = in.end();

        Element::iterator outIter = out.begin();
        Element::iterator outEnd = out.end();

        while(inIter != inEnd && outIter != outEnd)
        {
            if(outIter->first == inIter->first)
            {
                outIter->second += factor * inIter->second;
                ++inIter;
                ++outIter;
            }
            else if(outIter->first < inIter->first)
            {
                ++outIter;
            }
            else
            {
                out.assign(inIter->first, factor * inIter->second);
                ++inIter;
            }
        }

        while(inIter != inEnd)
        {
            out.assign(inIter->first, factor * inIter->second);
            ++inIter;
        }
    }

    Matrix makeMixingFactor(const vector<casa::Array<casa::DComplex> > &coeff,
        const Location &start, const Location &end, unsigned int baseline,
        unsigned int correlation, unsigned int row, unsigned int column)
    {
        const unsigned int nFreq = end.first - start.first + 1;
        const unsigned int nTime = end.second - start.second + 1;

        Matrix factor(makedcomplex(0.0, 0.0), nFreq, nTime, false);
        double *re = 0, *im = 0;
        factor.dcomplexStorage(re, im);

        // dims array: ndir x ndir x ncorr x nchan x nbl (minor -> major).
        casa::IPosition index(5, column, row, correlation, 0, baseline);
        for(unsigned int t = start.second; t <= end.second; ++t)
        {
            const casa::Array<casa::DComplex> &tmp = coeff[t];
            ASSERT(tmp.shape()(3) == nFreq);
            for(index(3) = start.first; index(3) <= static_cast<int>(end.first);
                ++index(3))
            {
                const casa::DComplex &weight = tmp(index);
                *re++ = real(weight);
                *im++ = imag(weight);
            }
        }

        return factor;
    }

    template <typename T_SAMPLE_MODIFIER>
    struct SampleProcessorComplex
    {
        static inline void process(Cell &cell, double weight, double reObs,
            double imObs, double reSim, double imSim, unsigned int nDerivative,
            double *reDerivative, double *imDerivative,
            const unsigned int *index, Statistics &stats)
        {
            // Modify the observed and simulated data depending on the solving
            // mode (complex, phase only, amplitude only).
            stats.start(Statistics::T_MODIFIER);
            T_SAMPLE_MODIFIER::process(weight, reObs, imObs, reSim, imSim,
                reDerivative, imDerivative, nDerivative);
            stats.stop(Statistics::T_MODIFIER);

            if(weight == 0.0)
            {
                return;
            }

            // Compute the residual.
            double reResidual = reObs - reSim;
            double imResidual = imObs - imSim;

            // Update the normal equations.
            stats.start(Statistics::T_MAKE_NORM);
            cell.solver.makeNorm(nDerivative, index, reDerivative, weight,
                reResidual);
            cell.solver.makeNorm(nDerivative, index, imDerivative, weight,
                imResidual);
            stats.stop(Statistics::T_MAKE_NORM);
        }
    };

    void initCells(const Location &start, const Location &end,
        const ParmGroup &solvables, size_t nCoeff,
        const EstimateOptions &options,
        vector<Cell> &cells)
    {
        vector<Cell>::iterator cell = cells.begin();
        for(CellIterator it(start, end); !it.atEnd(); ++it, ++cell)
        {
            // Processing has not completed yet.
            cell->done = false;

            // Initalize LSQ solver.
            cell->solver = casa::LSQFit(static_cast<casa::uInt>(nCoeff));
            configLSQSolver(cell->solver, options.lsqOptions());

            // Initialize coefficients.
            cell->coeff.resize(nCoeff);
            loadCoeff(*it, solvables, cell->coeff.begin());

//            // Clear RMS and sample counts.
//            cell->rms = 0;
//            cell->count = 0;

//            // Initialize L1 epsilon value.
//            cell->epsilonIdx = 0;
//            cell->epsilon = options.algorithm() == EstimateOptions::L1
//                ? options.epsilon(0) : 0.0;

//            // Initialize RMS threshold and deactivate outlier detection.
//            cell->flag = false;
//            cell->thresholdIdx = 0;
//            cell->threshold = numeric_limits<double>::infinity();
//            cell->outliers = 0;
        }
    }


    IterationStatus iterate(const Location &start, const Location &end,
        const ParmGroup &solvables, const EstimateOptions &options,
        vector<Cell> &cells)
    {
        IterationStatus status = {0, 0, 0, 0, 0};

        vector<Cell>::iterator cell = cells.begin();
        for(CellIterator it(start, end); !it.atEnd(); ++it, ++cell)
        {
            // If processing on the cell is already done, only update the status
            // counts and continue to the next cell.
            if(cell->done)
            {
                updateIterationStatus(*cell, status);
                continue;
            }

//            // Compute RMS.
//            if(cell->count > 0)
//            {
//                cell->rms = sqrt(cell->rms / cell->count);
//            }

//            // Turn outlier detection off by default. May be enabled later on.
//            cell->flag = false;

            // Perform a single iteration if the cell has not yet converged or
            // failed.
            if(!cell->solver.isReady())
            {
                // LSQFit::solveLoop() only returns false if the normal
                // equations are singular. This can also be seen from the result
                // of LSQFit::isReady(), so we don't update the iteration status
                // here but do skip the update of the solvables.
                casa::uInt rank;
                if(cell->solver.solveLoop(rank, &(cell->coeff[0]),
                    options.lsqOptions().useSVD))
                {
                    // Store the updated coefficient values.
                    storeCoeff(*it, solvables, cell->coeff.begin());
                }
            }

//            // Handle L1 restart with a different epsilon value.
//            if(cell->solver.isReady()
//                && options.algorithm() == EstimateOptions::L1
//                && cell->epsilonIdx < options.nEpsilon())
//            {
//                // Move to the next epsilon value.
//                ++cell->epsilonIdx;

//                if(cell->epsilonIdx < options.nEpsilon())
//                {
//                    // Re-initialize LSQ solver.
//                    size_t nCoeff = cell->coeff.size();
//                    cell->solver =
//                        casa::LSQFit(static_cast<casa::uInt>(nCoeff));
//                    configLSQSolver(cell->solver, options.lsqOptions());

//                    // Update epsilon value.
//                    cell->epsilon = options.epsilon(cell->epsilonIdx);
//                }
//            }

//            // Handle restart with a new RMS threshold value.
//            if(cell->solver.isReady()
//                && options.robust()
//                && cell->thresholdIdx < options.nThreshold())
//            {
//                // Re-initialize LSQ solver.
//                size_t nCoeff = cell->coeff.size();
//                cell->solver = casa::LSQFit(static_cast<casa::uInt>(nCoeff));
//                configLSQSolver(cell->solver, options.lsqOptions());

//                // Reset L1 state.
//                cell->epsilonIdx = 0;
//                cell->epsilon = options.algorithm() == EstimateOptions::L1
//                    ? options.epsilon(0) : 0.0;

//                // Compute new RMS threshold and activate outlier detection.
//                cell->threshold = options.threshold(cell->thresholdIdx)
//                    * cell->rms;
//                cell->flag = true;

//                // Move to the next threshold.
//                 ++cell->thresholdIdx;
//            }

//            // Log solution statistics.
//            if(!options.robust()
//                && (options.algorithm() == EstimateOptions::L2))
//            {
//                logCellStats(log, grid.getCell(*it), *cell);
//            }

            if(cell->solver.isReady())
            {
                cell->done = true;
            }
//            else
//            {
//                // If not yet converged or failed, reset state for the next
//                // iteration.
//                cell->rms = 0.0;
//                cell->count = 0;
//                cell->outliers = 0;
//            }

            updateIterationStatus(*cell, status);
        }

        return status;
    }

    void updateIterationStatus(const Cell &cell, IterationStatus &status)
    {
        // casa::LSQFit::isReady() is incorrectly labelled non-const.
        casa::LSQFit &solver = const_cast<casa::LSQFit&>(cell.solver);

        // Decode and record the solver status.
        switch(solver.isReady())
        {
            case casa::LSQFit::NONREADY:
                ++status.nActive;
                break;

            case casa::LSQFit::SOLINCREMENT:
            case casa::LSQFit::DERIVLEVEL:
                ++status.nConverged;
                break;

            case casa::LSQFit::MAXITER:
                ++status.nStopped;
                break;

            case casa::LSQFit::NOREDUCTION:
                ++status.nNoReduction;
                break;

            case casa::LSQFit::SINGULAR:
                ++status.nSingular;
                break;

            default:
                // This assert triggers if an casa::LSQFit ready code is
                // encountered that is not covered above. The most likely cause
                // is that the casa::LSQFit::ReadyCode enumeration has changed
                // in which case the code above needs to be changed accordingly.
                ASSERT(false);
                break;
        }
    }

    Statistics::Statistics()
    {
        fill(itsCounters, itsCounters + N_Counter, 0);
    }

    inline void Statistics::inc(Statistics::Counter counter)
    {
        ++itsCounters[counter];
    }

    inline void Statistics::inc(Statistics::Counter counter, size_t count)
    {
        itsCounters[counter] += count;
    }

    inline void Statistics::reset(Statistics::Counter counter)
    {
        itsCounters[counter] = 0;
    }

    inline void Statistics::reset(Statistics::Timer timer)
    {
        itsTimers[timer].reset();
    }

    inline void Statistics::start(Statistics::Timer timer)
    {
        itsTimers[timer].start();
    }

    inline void Statistics::stop(Statistics::Timer timer)
    {
        itsTimers[timer].stop();
    }

    void Statistics::reset()
    {
        fill(itsCounters, itsCounters + N_Counter, 0);

        for(size_t i = 0; i < N_Timer; ++i)
        {
            itsTimers[i].reset();
        }
    }

    string Statistics::counters() const
    {
        ostringstream oss;
        oss << "counters:";
        for(size_t i = 0; i < N_Counter; ++i)
        {
            oss << " " << theirCounterNames[i] << ": " << itsCounters[i];
        }
        return oss.str();
    }

    string Statistics::timers() const
    {
        ostringstream oss;
        oss << "timers:";
        for(size_t i = 0; i < N_Timer; ++i)
        {
            oss << " " << theirTimerNames[i] << ": " << toString(itsTimers[i]);
        }
        oss << " total/count/average";
        return oss.str();
    }

    string Statistics::theirCounterNames[Statistics::N_Counter] =
        {"all",
         "flagged",
         "zero weight",
         "invalid residual",
         "invalid derivative",
         "invalid weight",
         "outlier"};

    string Statistics::theirTimerNames[Statistics::N_Timer] =
        {"evaluate",
        "mix",
        "coeff map",
        "process cell",
        "modify sample",
        "condition eq"};

} // unnamed namespace

////        // intervals in VisBuffer coordinates? (NB. need map of cells -> samples)
////        // need cell offset, because start cell may not be (0,0).
////        // need map sample -> cells, index OK (sample coordinates)
////        // IDEA: sample -> cell map as 1's and 0's, so position (offset) independent.

////        for(t)
////        {
////            for(f)
////            {
////                for(d)
////                {
////                    derivatives[d] = *derive_p[d]++;
////                }

////                if(!flags[t][f])
////                {
////                    LSQFit &tmp = solvers[sampleMap[TIME][t] - offset[TIME]][sampleMap[FREQ][t] - offset[FREQ]];
////                    process_sample(tmp, 1.0 / cov[t][f], obs[t][f], *value_p, derivatives);
////                }

////                value_p++;
////            }
////        }

////        // merge all LSQFit objects...

} //# namespace BBS
} //# namespace LOFAR
