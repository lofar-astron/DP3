//# Patch.h: A set of sources for which direction dependent effects are assumed
//# to be equal.
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

#ifndef DPPP_PATCH_H
#define DPPP_PATCH_H

// \file
// A set of sources for which direction dependent effects are assumed to be
// equal.

#include <DPPP/ModelComponent.h>
#include <DPPP/Position.h>
#include <Common/lofar_vector.h>
#include <Common/lofar_string.h>

namespace LOFAR
{
namespace DPPP
{

// \addtogroup NDPPP
// @{

class Patch
{
public:
    typedef shared_ptr<Patch>       Ptr;
    typedef shared_ptr<const Patch> ConstPtr;

    template <typename T>
    Patch(const string &name, T first, T last);

    const string &name() const;
    const Position &position() const;

    size_t nComponents() const;
    ModelComponent::ConstPtr component(size_t i) const;

private:
    typedef vector<ModelComponent::ConstPtr>::const_iterator const_iterator;
    const_iterator begin() const;
    const_iterator end() const;

    void computePosition();

    string                              itsName;
    Position                            itsPosition;
    vector<ModelComponent::ConstPtr>    itsComponents;
};


// @}

// -------------------------------------------------------------------------- //
// - Implementation: Patch                                                  - //
// -------------------------------------------------------------------------- //

template <typename T>
Patch::Patch(const string &name, T first, T last)
    :   itsName(name),
        itsComponents(first, last)
{
    computePosition();
}

inline const string &Patch::name() const
{
    return itsName;
}

inline const Position &Patch::position() const
{
    return itsPosition;
}

inline size_t Patch::nComponents() const
{
    return itsComponents.size();
}

inline ModelComponent::ConstPtr Patch::component(size_t i) const
{
    return itsComponents[i];
}


} //# namespace DPPP
} //# namespace LOFAR

#endif
