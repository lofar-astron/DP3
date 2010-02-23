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

    // This class is an input step reading the data from a MeasurementSet.
    // At the beginning it finds out the shape of the data; i.e., the
    // number of correlations, channels, baselines, and time slots.
    // It requires the data to be regularly shaped.
    //
    // The object is constructed from the 'msin' keywords in the parset file.
    // Currently the following can be given:
    // <ul>
    //  <li> msin: name of the MS
    //  <li> msin.startchan: first channel to use [0]
    //  <li> msin.nchan: number of channels to use [all]
    //  <li> msin.useflag: use the existing flags? [yes]
    //  <li> msin.datacolumn: the data column to use [DATA]
    //  <li> msin.starttime: first time to use [first time in MS]
    //  <li> msin.endtime: last time to use [last time in MS]
    // </ul>
    //
    // If a time slot is missing, it is inserted with flagged data set to zero.
    // Missing time slots can also be detected at the beginning or end of the
    // MS by giving the correct starttime and endtime.
    // The correct UVW coordinates are calculated for inserted time slots.
    //
    // The process function only reads the data and flags to avoid that
    // too much data is kept in memory.
    // Other columns (like WEIGHT, UVW) can be read when needed by using the
    // appropriate DPInput::fetch function.

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
      virtual void show (std::ostream&) const;

      // If needed, show the flag counts.
      virtual void showCounts (std::ostream& os) const;

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

      // Get the rownrs for meta info of missing time slots.
      const casa::Vector<uint> getBaseRowNrs() const
        { return itsBaseRowNrs; }

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
      uint                itsNrInserted;    //# nr of inserted time slots
      double              itsInterval;
      double              itsFirstTime;
      double              itsLastTime;
      double              itsNextTime;
      casa::String        itsDataColName;
      casa::Slicer        itsSlicer;        //# slice in corr,chan
      casa::Slicer        itsPreAvgSlicer;  //# slice in chan,timeavg,bl
      DPBuffer            itsBuffer;
      UVWCalculator       itsUVWCalc;
      casa::Vector<uint>  itsBaseRowNrs;    //# rownrs for meta of missing times
      bool                itsCountFlags;
      FlagCounter         itsFlagCounter;
    };

  } //# end namespace
}

#endif
