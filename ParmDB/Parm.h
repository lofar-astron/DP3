//# Parm.h: Class giving access to a parameter
//#
//# Copyright (C) 2008
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
//# $Id: Parm.h 14038 2009-09-17 13:59:12Z diepen $

// @file
// @brief Class giving access to a parameter
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARM_H
#define LOFAR_PARMDB_PARM_H

#include "Grid.h"
#include "ParmSet.h"

#include <casacore/casa/Arrays/Array.h>

namespace DP3 {
namespace BBS {

  //# Forward declarations
  class ParmCache;
  class ParmValue;
  class ParmValueSet;
  class AxisMappingCache;

  // @ingroup ParmDB
  // @{

  // @brief Class giving access to a parameter
  // Parm makes it possible to get the parameter values for a given predict
  // domain. It uses a ParmCache to cache parameters for a work domain.
  // If a parameter is solvable, the solve grid must have been set before
  // getting the values and perturbed values. For solvable parameters it
  // is also possible to get or set the coefficients for each cell in the
  // solve grid.
  //
  // Parameter values are read from the ParmDB tables. If there are no values
  // for a parameter/domain, the appropriate value from the subtable with
  // default values is used. However, for a solvable parameter, the value
  // of the previous solve domain is used as long as that solve domain is
  // part of the work domain. Thus for a given work domain, the first solve
  // domain uses the default values, while the next domains use the values
  // of the previous domains.

  class Parm
  {
  public:
    // Construct the object for the given parmid.
    Parm (ParmCache&, ParmId parmid);

    // Construct the object for the given parm name.
    // It must have been added to the ParmSet used by the ParmCache.
    Parm (ParmCache&, const string& name);

    // Set the solve domains for the parm.
    // For an existing parm it must match the parm domains.
    // For a new parm it creates the domains.
    void setSolveGrid (const Grid& solveGrid);

    // Get the nr of coefficients.
    uint getCoeffSize (bool useMask=true);

    // Get the coefficients for the given location in the solve grid.
    // The solve grid must have been set before.
    std::vector<double> getCoeff (const Location&, bool useMask=true);

    // Get the errors for the given location in the solve grid.
    // The solve grid must have been set before.
    std::vector<double> getErrors (const Location&, bool useMask=true);

    // Set the coefficients for the given location in the solve grid.
    // If given, the errors are set too.
    // The solve grid must have been set before.
    // It sets the dirty flag, so the data are written when the ParmCache
    // is flushed.
    void setCoeff (const Location&, const double* values, uint nvalues,
                   const double* errors=0, bool useMask=true);

    // Revert to the original coefficients (as on disk).
    // (not implemented yet).
    // It clears the dirty flag.
    void revertCoeff();

    // Get the perturbations for the coefficients.
    // The possible mask is applied.
    const std::vector<double>& getPerturbations() const
      { return itsPerturbations; }

    // Get a particular perturbation.
    double getPerturbation (uint index)
      { return itsPerturbations.at (index); }

    // Get the result for the given grid. No perturbed values are calculated.
    // Normally the result has the same shape as the predict grid.
    // However, if there is a single constant value, it has shape [1,1].
    // In the future it might also get shape[nx,1] or [1,ny] in case the
    // values are constant in y or x.
    // <br>
    // The array will be resized if needed.
    // Note that it is possible to first create a correctly sized MeqMatrix
    // and create the Array from its raw storage. In that way the Array
    // data does not need to be copied to the MeqMatrix.
    // <br>
    // The argument <src>emptyResult</src> tells if an empty result can be
    // returned. Normally this is not the case (otherwise a default would
    // not be picked up), but in case of ParmFacade it is used.
    void getResult (casacore::Array<double>& result, const Grid& predictGrid,
                    bool emptyResult=false);

    // Get the values for the given predict grid.
    // The parm value is taken that contains the middle of a
    // predict interval.
    // If <src>perturb=true</src>, all perturbed values are calculated as well.
    // Otherwise only the first array in the vector is filled in.
    // As above, the shape of the array is normally [nx,ny],
    // but can be [1,1] if constant.
    void getResult (std::vector<casacore::Array<double> >& result,
                    const Grid& predictGrid, bool perturb);

    // Form the vector from values and mask.
    static std::vector<double> copyValues (const casacore::Array<double>& values,
                                      const casacore::Array<bool>& mask,
                                      bool useMask);

    // Evaluate the result for funklet coefficients.
    static void getResultCoeff (casacore::Array<double>* resultVec,
                                const Grid& predictGrid,
                                const ParmValueSet& pvset,
                                const std::vector<double>& perturbations,
                                AxisMappingCache& axisMappingCache);

    // Get the result for a single ParmValue with an array of scalars.
    static void getResultScalar (casacore::Array<double>& result,
                                 const Grid& predictGrid,
                                 const ParmValue& pval,
                                 AxisMappingCache& axisMappingCache);

    // Get the result for multiple ParmValues containing scalars.
    // If the <src>errors</src> argument is non-zero, the errors are filled too.
    static void getResultScalar (casacore::Array<double>& result,
                                 casacore::Array<double>* errors,
                                 const Grid& predictGrid,
                                 const ParmValueSet& pvset,
                                 AxisMappingCache& axisMappingCache);

    // Fill the result array partially for a single ParmValue.
    static void fillArrayPV (double* resData, int nrx, int stx, int sty,
                             int endx, int endy, const double* data,
                             const ParmValue& pval, const Grid& predictGrid);

    // Calculate the perturbations.
    void calcPerturbations();

  private:
    //# Data members
    ParmCache*     itsCache;
    ParmId         itsParmId;
    Grid           itsSolveGrid;
    std::vector<double> itsPerturbations;
  };

  // @}

} //# end namespace BBS
} //# end namspace LOFAR

#endif
