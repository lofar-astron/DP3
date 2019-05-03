//# ParmValue.h: A class containing the values of a parameter
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
//# $Id: ParmValue.h 16977 2010-12-20 08:40:36Z diepen $

// @file
// @brief A class containing the values of a parameter
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMVALUE_H
#define LOFAR_PARMDB_PARMVALUE_H

#include "Grid.h"
#include "Axis.h"

#include <casacore/casa/Arrays/Array.h>

#include <memory>
#include <string>
#include <vector>

namespace DP3 {
namespace BBS {

  // @ingroup ParmDB
  // @{

  // @brief A class containing the values of a parameter,
  // ParmValue holds the values of a given parameter and domain.
  // The object does not hold the name and domain info itself. Instead its
  // parent object ParmValueSet holds this information.
  //
  // The value is a 2-dim array holding scalar values or the coefficients
  // of a 2-dim funklet.
  // Thus if the parm is a scalar, a ParmValue object can hold the values of
  // multiple domains. If it's a funklet, only one domain is held.

  class ParmValue
  {
  public:
    // Define a shared pointer for this type.
    typedef std::shared_ptr<ParmValue> ShPtr;

    // Define the possible funklet types.
    enum FunkletType {
      // A constant scalar
      Scalar  = 0,
      // A polynomial
      Polc    = 1,
      // A polynomial of logs
      PolcLog = 2
    };

    // Construct with the given scalar value.
    explicit ParmValue (double value=0.);

    // Copy constructor makes a deep copy.
    ParmValue (const ParmValue&);

    ~ParmValue();

    // Assignment makes a deep copy.
    ParmValue& operator= (const ParmValue&);

    // Set as a single scalar value.
    void setScalar (double value);

    // Set as an array of coefficients.
    void setCoeff (const casacore::Array<double>&);

    // Set as an array of scalar values with the grid.
    // The shape of grid and values must match.
    void setScalars (const Grid&, const casacore::Array<double>&);

    // Set the errors.
    // They must have the same shape as the values, so the values must have
    // been set before.
    void setErrors (const casacore::Array<double>&);

    // Get the value shape.
    // <group>
    unsigned int nx() const
      { return static_cast<unsigned int>(itsValues.shape()[0]); }
    unsigned int ny() const
      { return static_cast<unsigned int>(itsValues.shape()[1]); }
    // </group>

    // Get the values.
    // <group>
    const casacore::Array<double>& getValues() const
      { return itsValues; }
    casacore::Array<double>& getValues()
      { return itsValues; }
    // </group>

    // Get the grid.
    const Grid& getGrid() const
      { return itsGrid; }

    // Are there errors? If false, the result of getErrors is undefined.
    bool hasErrors() const
      { return itsErrors != 0; }

    // Get the arrays with errors. Undefined if <src>getErrors()==false</src>.
    // <group>
    const casacore::Array<double>& getErrors() const
      { return *itsErrors; }
    casacore::Array<double>& getErrors()
      { return *itsErrors; }
    // </group>

    // Get/set the rowid to remember where the value is stored in the ParmDB.
    // <group>
    int getRowId() const
      { return itsRowId; }
    void setRowId (int rowId)
      { itsRowId = rowId; }
    void clearRowId()
      { itsRowId = -1; }
    // </group>
    
    // Return the scaled coefficients of a 2D polynomial using the
    // given offset and scale factor.
    static casacore::Matrix<double> scale2 (const casacore::Matrix<double>& coeff,
                                        double offx, double offy,
                                        double scalex, double scaley);

    // If needed rescale polynomial coefficients from the old domain to the
    // new domain given by the start and end values.
    // It returns true if rescaling was actually done.
    bool rescale (double sx, double ex, double sy, double ey,
                  const Box& oldDomain);

  private:
    // Make a deep copy of that.
    void copyOther (const ParmValue& that);

    // Fill Pascal's triangle till the given order.
    // The matrix will be resized as needed.
    static void fillPascal (casacore::Matrix<double>& pascal, int order);

    /// Data members.
    Grid                 itsGrid;          //# grid of the values
    casacore::Array<double>  itsValues;        //# scalar values or funklet coeff
    casacore::Array<double>* itsErrors;
    int                  itsRowId;         //# rowid in ParmDB
  };



  // @brief A class holding information of multiple domains of a parameter.
  // ParmValueSet holds the information of multiple domains of a parameter.
  // It has a grid defining the domains held.

  class ParmValueSet
  {
  public:

    // Create a parameterset with the given default
    // parm value (which is by default a scalar).
    // If the funklet type is a scalar, the value in the default must contain
    // one value only.
    // It is possible to specify the domain on which a funklet is scaled.
    explicit ParmValueSet (const ParmValue& defaultValue = ParmValue(),
                           ParmValue::FunkletType = ParmValue::Scalar,
                           double perturbation = 1e-6,
                           bool pertRel = true,
                           const Box& scaleDomain = Box());

    // Create the parameterset for the given domain grid and ParmValue objects.
    // If the funklet type is a scalar, the values in the ParmValues must
    // contain one value only.
    ParmValueSet (const Grid& domainGrid,
                  const std::vector<ParmValue::ShPtr>& values,
                  const ParmValue& defaultValue = ParmValue(),
                  ParmValue::FunkletType type = ParmValue::Scalar,
                  double perturbation = 1e-6,
                  bool pertRel = true);

    // Copy constructor.
    ParmValueSet (const ParmValueSet&);

    // Assignment.
    ParmValueSet& operator= (const ParmValueSet&);

    // Create the parameterset for the given grid from the given ParmValueSet
    // which should have only one ParmValue.
    // This is meant for solvable parameters using a default value.
    void setSolveGrid (const Grid& solveGrid);

    // Get the funklet type.
    ParmValue::FunkletType getType() const
      { return itsType; }

    // Get/set the mask telling which coefficients are solvable.
    // The array can be empty meaning that all coefficients are solvable.
    // <group>
    const casacore::Array<bool>& getSolvableMask() const
      { return itsSolvableMask; }
    void setSolvableMask (const casacore::Array<bool>& mask)
      { itsSolvableMask.assign (mask); }
    // <group>

    // Get the perturbation value.
    double getPerturbation() const
      { return itsPerturbation; }

    // Is the perturbation relative or absolute?
    bool getPertRel() const
      { return itsPertRel; }

    // Get access to the grid info, so a domain can be looked up.
    const Grid& getGrid() const
      { return itsDomainGrid; }

    // Get the nr of ParmValues.
    unsigned int size() const
      { return itsValues.size(); }

    // No ParmValues?
    bool empty() const
      { return itsValues.size() == 0; }

    // Get the default ParmValue.
    const ParmValue& getDefParmValue() const
      { return itsDefaultValue; }

    // Get access to the scale domain.
    // <group>
    const Box& getScaleDomain() const
      { return itsScaleDomain; }
    void setScaleDomain (const Box& domain)
      { itsScaleDomain = domain; }
    // </group>

    // Get the first ParmValue. If there are no ParmValues, the default
    // ParmValue is returned.
    const ParmValue& getFirstParmValue() const;

    // Get the i-th ParmValue.
    // <group>
    const ParmValue& getParmValue (int i) const
      { return *(itsValues[i]); }
    ParmValue& getParmValue (int i)
      { return *(itsValues[i]); }
    // </group>

    // Get/set the dirty flag.
    // The dirty flag has to be set when a new value is given to a ParmValue. 
    // It indicates that the value has to be written later on.
    // When written, the flag will be cleared.
    // <group>
    bool isDirty() const
      { return itsDirty; }
    void setDirty (bool dirty=true)
      { itsDirty = dirty; }
    // </group>

  private:
    // Helper functions for setSolveGrid.
    // <group>
    void createValues (const Grid& solveGrid);
    void checkGrid (const Grid& solveGrid);
    void addValues (const Grid& solveGrid);
    void addCoeffValues (const Grid& solveGrid);
    ParmValue::ShPtr copyParmCoeff (const ParmValue::ShPtr& pval);
    // </group>

    /// Data members.
    ParmValue::FunkletType itsType;
    double                 itsPerturbation;
    bool                   itsPertRel;
    casacore::Array<bool>      itsSolvableMask;
    Grid                   itsDomainGrid;
    std::vector<ParmValue::ShPtr> itsValues;
    ParmValue              itsDefaultValue;
    Box                    itsScaleDomain;
    bool                   itsDirty;
  };

  // @}

} //# end namespace BBS
} //# end namspace LOFAR

#endif
