//# EstimateNew.h: Estimate Jones matrices for several directions and stations
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

#ifndef DPPP_ESTIMATENEW_H
#define DPPP_ESTIMATENEW_H

// \file
// Estimate Jones matrices for several directions simultaneously. A separate
// data stream is used for each direction. The mixing coefficients quantify the
// influence of each direction on each of the other directions (including time
// and frequency smearing).

#include <DPPP/Baseline.h>
#include <DPPP/Cursor.h>
#include <Common/lofar_complex.h>
#include <Common/lofar_vector.h>

//# Use Block<bool> instead of vector<bool> (because testing bits is slower).
#include <casa/Containers/Block.h>

namespace LOFAR {
  namespace DPPP {

    // \addtogroup NDPPP
    // @{

    // Estimate Jones matrices for several directions simultaneously. A separate
    // data stream is used for each direction. The mixing coefficients quantify the
    // influence of each direction on each of the other directions (including time
    // and frequency smearing).
    class EstimateNew
    {
    public:

      EstimateNew();

      // Update the object and size its internal buffers.
      void update (size_t maxndir, size_t nBaseline, size_t nStation,
                   size_t nChannel, size_t maxIter, bool propagateSolution);

      // \param[in]   data
      // Vector of length \p nDirection of cursors for 3-D buffers of observed
      // visiblity data of shape (\p nBaseline, \p nChannel, 4).
      // \param[in]   model
      // Vector of length \p nDirection of cursors for 3-D buffers of simulated
      // visiblity data of shape (\p nBaseline, \p nChannel, 4).
      // \param[in]   baselines
      // A cursor for a 1-D buffer of baselines of shape (\p nBaseline).
      // \param[in]   flag
      // A cursor for a 3-D buffer of observed visibility flags of shape
      // (\p nBaseline, \p nChannel, 4).
      // \param[in]   weight
      // A cursor for a 3-D buffer of observed visibility weights of shape
      // (\p nBaseline, \p nChannel, 4).
      // \param[in]   mix
      // A cursor for a 5-D buffer of mixing weights of shape
      // (\p nBaseline, \p nChannel, 4, \p nDirection, \p nDirection).
      // \param[in]   solveBoth
      // True = only use baseline if both stations are solvable
      //
      // <br>Note that the cursors are passed by value, so a copy is made.
      // In this way no reset of the cursor is needed.
      bool estimate (const vector<vector<int> >& unknownsIndex,
                     const vector<uint>& srcSet,
                     const_cursor<Baseline> baselines,
                     vector<const_cursor<fcomplex> > data,
                     vector<const_cursor<dcomplex> > model,
                     const_cursor<bool> flag,
                     const_cursor<float> weight,
                     const_cursor<dcomplex> mix,
		     double defaultGain,
                     bool solveBoth,
                     uint verbose);

      // Get the last solution.
      // It contains zeroes for the direction-stations not solved for.
      const vector<double>& getSolution() const
        { return itsSolution; }

      // Get the nr of iterations used.
      size_t nIterations() const
        { return itsNrIter; }

    private:
      // Initialize the solution. Nr must be a multiple of 8.
      // The diagonal is set to (diag,0) or (1e-8,0), off-diagonal to (0,0).
      void initSolution (const vector<vector<int> >& unknownsIndex,
                         const vector<uint>& srcSet,
			 double defaultGain);

      // Clear the solution for unsolvable stations
      // (essentially changing 1e-8 to 0).
      void clearNonSolvable (const vector<vector<int> >& unknownsIndex,
                             const vector<uint>& srcSet);

      // Update itsSolution from itsUnknowns for the unknowns to be used.
      void fillSolution (const vector<vector<int> >& unknownsIndex,
                         const vector<uint>& srcSet);

      // Fill itsDerivIndex for the unknowns of the given baseline
      // to be able to pass the equations to LSQFit::makeNorm.
      // It returns the number of unknowns.
      uint fillDerivIndex (size_t ndir,
                           const vector<vector<int> >& unknownsIndex,
                           const Baseline& baseline);

      //# Data members
      size_t itsNrBaselines;
      size_t itsNrStations;
      size_t itsNrChannels;
      size_t itsMaxIter;
      size_t itsNrIter;
      size_t itsNrDir;
      bool   itsPropagateSolution;
      casa::Block<bool>  itsSolveStation;  //# solve station i?
      vector<casa::uInt> itsDerivIndex;    //# index for LSQFit::makeIndex
      vector<double>     itsUnknowns;
      vector<double>     itsSolution;
      vector<dcomplex>   itsM;
      vector<dcomplex>   itsdM;
      vector<double>     itsdR;
      vector<double>     itsdI;
    };

    // @}

  } //# namespace DPPP
} //# namespace LOFAR

#endif
