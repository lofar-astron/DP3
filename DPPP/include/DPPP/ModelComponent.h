//# ModelComponent.h: Base class for model components.
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

#ifndef DPPP_MODELCOMPONENT_H
#define DPPP_MODELCOMPONENT_H

// \file
// Base class for model components.

#include <Common/lofar_smartptr.h>

namespace LOFAR
{
namespace DPPP
{

class ModelComponentVisitor;
class Position;

// \addtogroup NDPPP
// @{

class ModelComponent
{
public:
    typedef shared_ptr<ModelComponent>          Ptr;
    typedef shared_ptr<const ModelComponent>    ConstPtr;

    virtual ~ModelComponent();
    virtual const Position &position() const = 0;
    virtual void accept(ModelComponentVisitor&) const = 0;
};

// @}

} //# namespace DPPP
} //# namespace LOFAR

#endif
