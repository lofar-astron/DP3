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
    // State kept for a single cell in the solution grid.
    struct Cell
    {
        // LSQ solver and current estimates for the coefficients.
        casa::LSQFit    solver;
        vector<double>  coeff;
    };

    template <typename T_SAMPLE_MODIFIER>
    struct SampleProcessorComplex;

    JonesMatrix mix(const vector<JonesMatrix> &in,
        const vector<casa::Array<casa::DComplex> > &coeff,
        unsigned int target,
        unsigned int nsources,
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
        unsigned int nsources,
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
        vector<Cell> &cells);

    template <typename T_SAMPLE_PROCESSOR>
    void equate2(const Location &start, const Location &end, size_t blIndex,
        const vector<DPPP::DPBuffer> &obs, const JonesMatrix &sim,
        const vector<pair<size_t, size_t> > &correlationMap,
        const vector<Interval<size_t> > (&cellMap)[2],
        const map<PValueKey, unsigned int> &coeffMap,
        vector<Cell> &cells);

    void initCells(const Location &start, const Location &end,
        const ParmGroup &solvables, size_t nCoeff,
        const EstimateOptions &options,
        vector<Cell> &cells);

    // Perform a single iteration for all cells in the range [\p start, \p end)
    // that have not yet converged or failed, updating the coefficient values
    // to the new estimates found.
    // \pre The range starting at \p cell should contain exactly one Cell
    // instance for each cell in the range [\p start, \p end].
    bool iterate(const Location &start, const Location &end,
        const ParmGroup &solvables, const EstimateOptions &options,
        vector<Cell> &cells);
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
    const size_t nDirections = buffers.size();
    const size_t nModels = models.size();
    {
        ASSERT(nDirections >= nModels && nModels > 0);
        ASSERT(int(nDirections) == coeff[0].shape()[0]);
        ASSERT(int(nModels) == coeff[0].shape()[1]);

        CorrelationSeq tmp;
        tmp.append(Correlation::XX);
        tmp.append(Correlation::XY);
        tmp.append(Correlation::YX);
        tmp.append(Correlation::YY);

        for(size_t i = 0; i < nModels; ++i)
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

    // Compute a mapping from cells of the solution grid to cell intervals
    // in the evaluation grid.
    vector<Interval<size_t> > cellMap[2];
    Interval<size_t> domain[2];
    domain[FREQ] = makeAxisMap(solGrid[FREQ], visGrid[FREQ],
        back_inserter(cellMap[FREQ]));
    domain[TIME] = makeAxisMap(solGrid[TIME], visGrid[TIME],
        back_inserter(cellMap[TIME]));

    ParmGroup solvables;
    for(size_t i = 0; i < nModels; ++i)
    {
        ParmGroup tmp = models[i]->solvables();
        solvables.insert(tmp.begin(), tmp.end());
    }

    // Make coefficient map.
    map<PValueKey, unsigned int> coeffMap;
    makeCoeffMap(solvables, inserter(coeffMap, coeffMap.begin()));
    LOG_DEBUG_STR("No. of coefficients to estimate: " << coeffMap.size());

    // Assign solution grid to solvables.
    ParmManager::instance().setGrid(solGrid, solvables);

    // ---------------------------------------------------------------------
    // Process each chunk of cells in a loop.
    // ---------------------------------------------------------------------

    // Clip chunk size to the size of the solution grid.
    size_t chunkSize = options.chunkSize() == 0 ? solGrid[TIME]->size()
        : std::min(options.chunkSize(), solGrid[TIME]->size());

    // Allocate cells.
    vector<Cell> cells(solGrid[FREQ]->size() * chunkSize);

    // Compute the number of cell chunks to process.
    size_t nChunks = (solGrid[TIME]->size() + chunkSize - 1) / chunkSize;

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
            timerChunk.stop();
            continue;
        }

        // Ensure a model value is computed for all the visibility samples
        // within the chunk.
        Location reqStart(cellMap[FREQ][chunkStart.first].start,
            cellMap[TIME][chunkStart.second].start);
        Location reqEnd(cellMap[FREQ][chunkEnd.first].end,
            cellMap[TIME][chunkEnd.second].end);
        for(size_t i = 0; i < nModels; ++i)
        {
            models[i]->setEvalGrid(visGrid.subset(reqStart, reqEnd));
        }

        // Initialize a cell instance for each cell in [chunkEnd,
        // chunkStart].
        initCells(chunkStart, chunkEnd, solvables, coeffMap.size(), options,
            cells);

        typedef SampleProcessorComplex<SampleModifierComplex> SampleProcessor;

        bool done = false;
        unsigned int nIterations = 0;
        while(!done)
        {
            // Construct normal equations from the data and an evaluation of
            // the model based on the current coefficient values.
            timerEquate.start();
            equate<SampleProcessor>(chunkStart, chunkEnd, buffers, models,
                coeff, blMap, crMap, cellMap, coeffMap, cells);
            timerEquate.stop();

            // Perform a single iteration.
            timerIterate.start();
            done = iterate(chunkStart, chunkEnd, solvables, options, cells);
            timerIterate.stop();

            // Notify model that solvables have changed.
            for(size_t i = 0; i < nModels; ++i)
            {
                models[i]->solvablesChanged();
            }

            // Update iteration count.
            ++nIterations;
        }
        timerChunk.stop();

        // Output statistics and timers.
        const size_t nCells = (chunkEnd.second - chunkStart.second + 1)
            * (chunkEnd.first - chunkStart.first + 1);
        LOG_DEBUG_STR("chunk: " << (chunk + 1) << "/" << nChunks
            << " cells: " << nCells << " iterations: " << nIterations);
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
    unsigned int nsources)
{
    // =========================================================================
    // CHECK PRECONDITIONS
    // =========================================================================
    {
        ASSERT(nsources <= models.size());
        ASSERT(buffer.size() > 0);
        ASSERT(int(target) < coeff[0].shape()[1]);

        CorrelationSeq tmp;
        tmp.append(Correlation::XX);
        tmp.append(Correlation::XY);
        tmp.append(Correlation::YX);
        tmp.append(Correlation::YY);

        for(size_t i = 0; i < nsources; ++i)
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

    for(size_t i = 0; i < nsources; ++i)
    {
        models[i]->setEvalGrid(visGrid);
    }

    subtract2(buffer, models, coeff, target, nsources, blMap, crMap);
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
    vector<Cell> &cells)

{
    ASSERT(buffers.size() >= models.size());

    const size_t nDirections = buffers.size();
    const size_t nModels = models.size();

    vector<JonesMatrix> sim(nModels);
    vector<JonesMatrix> mixed(nDirections);

    typedef vector<pair<size_t, size_t> >::const_iterator bl_iterator;
    for(bl_iterator it = baselineMap.begin(), it_end = baselineMap.end();
        it != it_end; ++it)
    {
        // Evaluate models.
        for(size_t i = 0; i < nModels; ++i)
        {
            sim[i] = models[i]->evaluate(it->second);
        }

        // Mix
        Location visStart(cellMap[FREQ][start.first].start,
            cellMap[TIME][start.second].start);
        Location visEnd(cellMap[FREQ][end.first].end,
            cellMap[TIME][end.second].end);
        mix(sim, mixed, coeff, visStart, visEnd, it->first);

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
                cellMap, coeffMap, cells);
        }
    } // baselines
}

template <typename T_SAMPLE_PROCESSOR>
void equate2(const Location &start, const Location &end, size_t blIndex,
    const vector<DPPP::DPBuffer> &obs, const JonesMatrix &sim,
    const vector<pair<size_t, size_t> > &correlationMap,
    const vector<Interval<size_t> > (&cellMap)[2],
    const map<PValueKey, unsigned int> &coeffMap,
    vector<Cell> &cells)
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
                            &(coeffIndex[0]));
                    }

                    ++index;
                }

                index -= (freqInterval.end - freqInterval.start + 1);
                index += nFreq;
            }
        }
    }
}

void subtract2(vector<DPPP::DPBuffer> &buffer,
    const vector<MeasurementExpr::Ptr> &models,
    const vector<casa::Array<casa::DComplex> > &coeff,
    unsigned int target,
    unsigned int nsources,
    const vector<pair<size_t, size_t> > &baselineMap,
    const vector<pair<size_t, size_t> > &correlationMap)
{
    vector<JonesMatrix> sim(nsources);

    typedef vector<pair<size_t, size_t> >::const_iterator
        index_map_iterator;

    for(index_map_iterator bl_it = baselineMap.begin(),
        bl_end = baselineMap.end(); bl_it != bl_end; ++bl_it)
    {
        // Evaluate models.
        for(unsigned int i = 0; i < nsources; ++i)
        {
            sim[i] = models[i]->evaluate(bl_it->second);
        }

        // Mix
        JonesMatrix mixed = mix(sim, coeff, target, nsources, bl_it->first);

        // Subtract.
        for(index_map_iterator cr_it = correlationMap.begin(),
            cr_end = correlationMap.end(); cr_it != cr_end; ++cr_it)
        {
            Matrix crMixed = mixed.element(cr_it->second).value();
            ASSERT(!crMixed.isNull());

            const unsigned int nFreq = crMixed.nx();
            const unsigned int nTime = crMixed.ny();
            ASSERT(crMixed.isComplex() && crMixed.isArray());
            ASSERTSTR(nTime == buffer.size(), "nTime: " << nTime << " buffer size: " << buffer.size());

            const double *mixed_re = 0, *mixed_im = 0;
            crMixed.dcomplexStorage(mixed_re, mixed_im);

            for(vector<DPBuffer>::iterator buffer_it = buffer.begin(),
                buffer_end = buffer.end(); buffer_it != buffer_end;
                ++buffer_it)
            {
                casa::Cube<casa::Complex> &data = buffer_it->getData();
                ASSERT(data.shape()(1) == int(nFreq));

                for(size_t i = 0; i < nFreq; ++i)
                {
                    data(cr_it->first, i, bl_it->first) -=
                        makedcomplex(*mixed_re++, *mixed_im++);
                } // frequency
            } // time
        } // correlations
    } // baselines
}

JonesMatrix mix(const vector<JonesMatrix> &in,
    const vector<casa::Array<casa::DComplex> > &coeff,
    unsigned int target,
    unsigned int nsources,
    unsigned int baseline)
{
    const unsigned int nFreq = coeff.front().shape()(3);
    const unsigned int nTime = coeff.size();

    ASSERT(nFreq >= 1 && nTime >= 1);
    const Location start(0, 0);
    const Location end(nFreq - 1, nTime - 1);

    Matrix out[4];
    for(unsigned int i = 0; i < nsources; ++i)
    {
        for(unsigned int correlation = 0; correlation < 4; ++correlation)
        {
            // Exchanged target and i, because we want the effect of
            // direction i on the target direction.
            Matrix weight = makeMixingFactor(coeff, start, end, baseline,
                correlation, target, i);

            const Matrix sim = in[i].element(correlation).value();

            ASSERTSTR(sim.nx() == weight.nx(), "sim: " << sim.nx() << " weight: " << weight.nx());
            ASSERTSTR(sim.ny() == weight.ny(), "sim: " << sim.ny() << " weight: " << weight.ny());
            if(out[correlation].isNull())
            {
                out[correlation] = weight * sim;
            }
            else
            {
                out[correlation] += weight * sim;
            }
        } // correlations
    } // nsources

    JonesMatrix result(out[0], out[1], out[2], out[3]);
    result.setFlags(mergeFlags(in));
    return result;
}

void mix(const vector<JonesMatrix> &in, vector<JonesMatrix> &out,
    const vector<casa::Array<casa::DComplex> > &coeff,
    const Location &start, const Location &end, unsigned int baseline)
{
    // dims array: ndir x nmodel x ncorr x nchan x nbl (minor -> major).
    // better dims: nbl x ncorr x ndir x nmodel x nchan ???

    const unsigned int nModels = in.size();
    const unsigned int nDirections = coeff[0].shape()[0];
    ASSERT(nDirections == out.size());
    ASSERT(int(nModels) == coeff[0].shape()[1]);

    FlagArray flags = mergeFlags(in);

    for(unsigned int i = 0; i < nDirections; ++i)
    {
        Element element[4];
        for(unsigned int j = 0; j < nModels; ++j)
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

    // dims coeff: ndir x nmodel x ncorr x nchan x nbl (minor -> major).
    //  nmodel = nr of directions with source model (thus excl. target)
    casa::IPosition index(5, row, column, correlation, 0, baseline);
    for(unsigned int t = start.second; t <= end.second; ++t)
    {
        const casa::Array<casa::DComplex> &tmp = coeff[t];
        ASSERTSTR(tmp.shape()(3) == int(nFreq), "nFreq: "
		  << nFreq << ' ' << tmp.shape());
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
        const unsigned int *index)
    {
        // Modify the observed and simulated data depending on the solving
        // mode (complex, phase only, amplitude only).
        T_SAMPLE_MODIFIER::process(weight, reObs, imObs, reSim, imSim,
            reDerivative, imDerivative, nDerivative);

        if(weight == 0.0)
        {
            return;
        }

        // Compute the residual.
        double reResidual = reObs - reSim;
        double imResidual = imObs - imSim;

        // Update the normal equations.
        cell.solver.makeNorm(nDerivative, index, reDerivative, weight,
            reResidual);
        cell.solver.makeNorm(nDerivative, index, imDerivative, weight,
            imResidual);
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
        // Initalize LSQ solver.
        cell->solver = casa::LSQFit(static_cast<casa::uInt>(nCoeff));
        configLSQSolver(cell->solver, options.lsqOptions());

        // Initialize coefficients.
        cell->coeff.resize(nCoeff);
        loadCoeff(*it, solvables, cell->coeff.begin());
    }
}

bool iterate(const Location &start, const Location &end,
    const ParmGroup &solvables, const EstimateOptions &options,
    vector<Cell> &cells)
{
    const size_t nCellFreq = end.first - start.first + 1;
    const size_t nCellTime = end.second - start.second + 1;
    const size_t nCell = nCellFreq * nCellTime;

    bool done = true;
#pragma omp parallel for
    for(size_t i = 0; i < nCell; ++i)
    {
        Cell &cell = cells[i];

        // If processing on the cell is already done, only update the status
        // counts and continue to the next cell.
        if(cell.solver.isReady())
        {
            continue;
        }

        // Perform a single iteration if the cell has not yet converged or
        // failed.
        //
        // LSQFit::solveLoop() only returns false if the normal
        // equations are singular. This can also be seen from the result
        // of LSQFit::isReady(), so we don't update the iteration status
        // here but do skip the update of the solvables.
        casa::uInt rank;
        if(cell.solver.solveLoop(rank, &(cell.coeff[0]),
            options.lsqOptions().useSVD))
        {
            // Store the updated coefficient values.
            storeCoeff(Location(start.first + i % nCellFreq, start.second + i
                / nCellFreq), solvables, cell.coeff.begin());
        }

        if(!cell.solver.isReady())
        {
            done = false;
        }
    }

    return done;
}

} // unnamed namespace

} //# namespace BBS
} //# namespace LOFAR
