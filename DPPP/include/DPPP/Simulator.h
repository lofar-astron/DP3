//# Simulator.h: Compute visibilities for different model components types
//# (implementation of ModelComponentVisitor).
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

#ifndef DPPP_SIMULATOR_H
#define DPPP_SIMULATOR_H

// \file
// Compute visibilities for different model components types (implementation of
// ModelComponentVisitor).

#include <DPPP/Baseline.h>
#include <DPPP/ModelComponent.h>
#include <DPPP/ModelComponentVisitor.h>
#include <DPPP/Position.h>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>

#include <Common/lofar_complex.h>
#include <Common/lofar_vector.h>

namespace LOFAR
{
namespace DPPP
{

// \addtogroup NDPPP
// @{

class Simulator: public ModelComponentVisitor
{
public:
    Simulator(const Position &reference, size_t nStation, size_t nBaseline,
        size_t nChannel, const casacore::Vector<Baseline>& baselines,
        const casacore::Vector<double>& freq, const casacore::Matrix<double>& uvw,
        casacore::Cube<dcomplex>& buffer, bool stokesIOnly=false);

    void simulate(const ModelComponent::ConstPtr &component);

private:
    virtual void visit(const PointSource &component);
    virtual void visit(const GaussianSource &component);

private:
    Position                     itsReference;
    size_t                       itsNStation, itsNBaseline, itsNChannel;
    bool                         itsStokesIOnly;
    const casacore::Vector<Baseline> itsBaselines;
    const casacore::Vector<double>   itsFreq;
    const casacore::Matrix<double>   itsUVW;
    casacore::Cube<dcomplex>         itsBuffer;
    casacore::Matrix<dcomplex>       itsShiftBuffer;
    casacore::Matrix<dcomplex>       itsSpectrumBuffer;
};

// @}

} //# namespace DPPP
} //# namespace LOFAR

#endif
