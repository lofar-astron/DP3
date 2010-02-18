//# MSReader.h: DPPP step reading from an MS
//# Copyright (C) 2010
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
//#
//# @author Ger van Diepen

#ifndef DPPP_MSREADER_H
#define DPPP_MSREADER_H

// @file
// @brief DPPP step reading from an MS

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/UVWCalculator.h>
#include <DPPP/FlagCounter.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/RefRows.h>
#include <casa/Arrays/Slicer.h>
#include <Common/lofar_vector.h>

namespace LOFAR {
  class ParameterSet;

  namespace DPPP {

    // @ingroup DPPP

    class MSReader: public DPInput
    {
    public:
      // Construct the object for the given MS.
      // Parameters are obtained from the parset using the given prefix.
      MSReader (const std::string& msName,
                const ParameterSet&, const string& prefix);

      virtual ~MSReader();

      // Process the next data chunk.
      // It returns false when at the end.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the average info (by initializing it).
      virtual void updateAverageInfo (AverageInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&);

      // Read the UVW at the given row numbers.
      virtual casa::Matrix<double> getUVW (const casa::RefRows& rowNrs);

      // Read the weights at the given row numbers.
      virtual casa::Cube<float> getWeights (const casa::RefRows& rowNrs);

      // Read the preAvg flags (LOFAR_PREAVG_FLAG) at the given row numbers.
      // It returns a 3-dim array [norigchan, ntimeavg, nbaseline].
      // If undefined, an empty array is returned.
      virtual casa::Cube<bool> getPreAvgFlags (const casa::RefRows& rowNrs);

      // Read the given data column at the given row numbers.
      virtual casa::Cube<casa::Complex> getData (const casa::String& columnName,
                                                 const casa::RefRows& rowNrs);

      // Write the flags at the given row numbers.
      void putFlags (const casa::RefRows& rowNrs,
                     const casa::Cube<bool>& flags);

      // Get the main MS table.
      casa::Table& table()
        { return itsMS; }

      // Get the name of the MS.
      const std::string& msName() const
        { return itsMS.tableName(); }

    private:
      // Prepare the access to the MS.
      // Return the first and last time and the interval.
      void prepare (double& firstTime, double& lastTime,
                    double& interval);

      // Calculate the UVWs for a missing time slot.
      void calcUVW();

      // Setup the UVW calculator.
      void setupUVWCalc();

      casa::Table         itsMS;
      casa::TableIterator itsIter;
      casa::Vector<int>   itsAnt1;
      casa::Vector<int>   itsAnt2;
      bool                itsHasWeightSpectrum;
      bool                itsHasPreAvgFlags;
      bool                itsUseFlags;
      uint                itsStartChan;
      uint                itsNrInserted;
      double              itsInterval;
      double              itsFirstTime;
      double              itsLastTime;
      double              itsNextTime;
      casa::String        itsDataColName;
      vector<int>         itsFlagChannels;
      casa::Slicer        itsSlicer;        //# slice in corr,chan
      casa::Slicer        itsFullSlicer;    //# slice in corr,chan,bl
      DPBuffer            itsBuffer;
      UVWCalculator       itsUVWCalc;
      bool                itsCountFlags;
      FlagCounter         itsFlagCounter;
    };

  } //# end namespace
}

#endif
