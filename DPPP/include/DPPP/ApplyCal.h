//# ApplyCal.h: DPPP step class to apply a calibration correction to the data
//# Copyright (C) 2013
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
//# $Id: ApplyCal.h 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_APPLYCAL_H
#define DPPP_APPLYCAL_H

// @file
// @brief DPPP step class to apply a calibration correction to the data

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/H5Parm.h>
#include <ParmDB/ParmFacade.h>
#include <ParmDB/ParmSet.h>
#include <ParmDB/Parm.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <DPPP/FlagCounter.h>

namespace LOFAR {
  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class applying calibration parameters to the data.
    class ApplyCal: public DPStep
    {
    public:

      enum CorrectType {GAIN, FULLJONES, TEC, CLOCK, ROTATIONANGLE, SCALARPHASE, PHASE,
                        ROTATIONMEASURE, SCALARAMPLITUDE, AMPLITUDE};

      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      ApplyCal (DPInput*, const ParameterSet&, const string& prefix,
                bool substep=false, std::string predictDirection="");

      ApplyCal();

      virtual ~ApplyCal();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

      bool invert() {
        return itsInvert;
      }

      // Invert a 2x2 matrix in place
      static void invert (casacore::DComplex* v, double sigmaMMSE=0);

      // Apply a diagonal Jones matrix to the 2x2 visibilities matrix: A.V.B^H
      static void applyDiag (const casacore::DComplex* gainA,
                             const casacore::DComplex* gainB,
                             casacore::Complex* vis, float* weight, bool* flag,
                             uint bl, uint chan, bool updateWeights,
                             FlagCounter& flagCounter);

      // Apply a diagonal Jones matrix to the 2x2 visibilities matrix: A.V.B^H,
      // where the solution is equal for both polarizations
      static void applyScalar(const casacore::DComplex* gainA,
                              const casacore::DComplex* gainB,
                              casacore::Complex* vis, float* weight, bool* flag,
                              uint bl, uint chan, bool updateWeights,
                              FlagCounter& flagCounter);

      // Apply a full Jones matrix to the 2x2 visibilities matrix: A.V.B^H
      static void applyFull (const casacore::DComplex* gainA,
                             const casacore::DComplex* gainB,
                             casacore::Complex* vis, float* weight, bool* flag,
                             uint bl, uint chan, bool updateWeights,
                             FlagCounter& flagCounter);

      // Do the same as the combination of BBS + python script
      // covariance2weight.py (cookbook), except it stores weights per freq.
      // The diagonal of covariance matrix is transferred to the weights.
      // Note that the real covariance (mixing of noise terms after which they
      // are not independent anymore) is not stored.
      // The input covariance matrix C is assumed to be diagonal with elements
      // w_i (the weights), the result the diagonal of
      // (gainA kronecker gainB^H).C.(gainA kronecker gainB^H)^H
      static void applyWeights (const casacore::DComplex* gainA,
                                const casacore::DComplex* gainB,
                                float* weight);

    private:
      // Read parameters from the associated parmdb and store them in itsParms
      void updateParms (const double bufStartTime);

      // If needed, show the flag counts.
      virtual void showCounts (std::ostream&) const;

      void initDataArrays();

      // Check the number of polarizations in the parmdb or h5parm
      uint nPol(const std::string& parmName);

      static std::string correctTypeToString(CorrectType);
      static CorrectType stringToCorrectType(const string&);

      //# Data members.
      DPInput*         itsInput;
      DPBuffer         itsBuffer;
      string      itsName;
      string      itsParmDBName;
      bool             itsUseH5Parm;
      boost::shared_ptr<BBS::ParmFacade> itsParmDB;
      H5Parm           itsH5Parm;
      H5Parm::SolTab   itsSolTab;
      CorrectType      itsCorrectType;
      bool             itsInvert;
      uint             itsTimeSlotsPerParmUpdate;
      double           itsSigmaMMSE;
      bool             itsUpdateWeights;

      uint             itsCount; // number of steps

      // Expressions to search for in itsParmDB
      vector<casacore::String>   itsParmExprs;

      // parameters, numparms, antennas, time x frequency
      casacore::Cube<casacore::DComplex> itsParms;
      uint            itsTimeStep; // time step within current chunk
      uint            itsNCorr;
      double          itsTimeInterval;
      double          itsLastTime; // last time of current chunk
      FlagCounter     itsFlagCounter;
      bool            itsUseAP;      //# use ampl/phase or real/imag
      hsize_t         itsDirection;
      NSTimer         itsTimer;
    };

  } //# end namespace
}

#endif
