//# OpenMP.h: Interface to some common OpenMP functions.
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
//# $Id: OpenMP.h 27642 2013-12-04 08:05:36Z diepen $
//#
//# @author Ger van Diepen <diepen at astron dot nl>

#ifndef LOFAR_COMMON_OPENMP_H
#define LOFAR_COMMON_OPENMP_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace DP3 {
  namespace OpenMP {

    // Get the maximum number of threads.
    // OpenMP sets it to the env.var. OMP_NUM_THREADS. If undefined, it is
    // the number of cores.
    // If OpenMP is not used, 1 is returned.
    inline uint maxThreads()
    {
#ifdef _OPENMP
      return omp_get_max_threads();
#else
      return 1;
#endif
    }

    // Set the number of threads to use. Note it can be overridden
    // for a parallel section by 'omp parallel num_threads(n)'.
    // Nothing is done if OpenMP is not used.
#ifdef _OPENMP
    inline void setNumThreads (uint n)
      { omp_set_num_threads (n); }
#else
    inline void setNumThreads (uint)
    {}
#endif

    // Get the number of threads used in a parallel piece of code.
    // If OpenMP is not used, 1 is returned.
    inline uint numThreads()
    {
#ifdef _OPENMP
      return omp_get_num_threads();
#else
      return 1;
#endif
    }

    // Get the thread number (0 till numThreads).
    // If OpenMP is not used, 0 is returned.
    inline uint threadNum()
    {
#ifdef _OPENMP
      return omp_get_thread_num();
#else
      return 0;
#endif
    }

    // Set if nested parallel sections are possible or not.
    // Nothing is done if OpenMP is not used.
#ifdef _OPENMP
    inline void setNested (bool nest)
      { omp_set_nested (nest); }
#else
    inline void setNested (bool)
    {}
#endif

    // Test if nested parallel sections are possible.
    // If OpenMP is not used, false is returned.
    inline bool nested()
    {
#ifdef _OPENMP
      return omp_get_nested();
#else
      return false;
#endif
    }

  } // end namespace
} // end namespace

#endif
