//# MWError.h: Basic exception for master/worker related errors
//#
//# Copyright (C) 2005
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
//# $Id: MWError.h 16886 2010-12-08 10:43:17Z diepen $

#ifndef LOFAR_LMWCOMMON_MWERROR_H
#define LOFAR_LMWCOMMON_MWERROR_H

// @file
// @brief Basic exception for master/worker related errors.
// @author Ger van Diepen (diepen AT astron nl)

namespace DP3 { namespace CEP {

  // @ingroup LMWCommon
  // @brief Basic exception for master/worker related errors.

  // This class defines the basic MW exception.
  // Only this basic exception is defined so far. In the future, some more 
  // fine-grained exceptions might be derived from it.
  typedef std::runtime_error MWError;

}} //# end namespaces

#endif
