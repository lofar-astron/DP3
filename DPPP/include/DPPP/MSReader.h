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

    // @ingroup NDPPP

    // This class is a DPInput step reading the data from a MeasurementSet.
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
    //
    // The data columns are handled in the following way:
    // <table>
    //  <tr>
    //   <td>TIME</td>
    //   <td>The time slot center of the current data (in MJD seconds).
    //       It is assumed that all data have the same interval, which is
    //       used to find missing time slots.
    //   </td>
    //  </tr>
    //  <tr>
    //   <td>DATA</td>
    //   <td>The visibility data as [ncorr,nchan,nbaseline]. Only the
    //       part given by startchan and nchan is read. If a time slot is
    //       inserted, all its data are zero.
    //   </td>
    //  </tr>
    //  <tr>
    //   <td>FLAG</td>
    //   <td>The data flags as [ncorr,nchan,nbaseline] (True is bad).
    //       They are read from the FLAG column. If a FLAG_ROW is set, all
    //       flags for that baseline will be set. Also the flag of data
    //       containing NaN or infinite numbers will be set.
    //       All flags of an inserted time slot are set.
    //   </td>
    //  </tr>
    //  <tr>
    //   <td>WEIGHT</td>
    //   <td>The data weights as [ncorr,nchan,nbaseline]. Column
    //       WEIGHT_SPECTRUM is used if present and containing valid data,
    //       otherwise column WEIGHT is used. The weights of an inserted
    //       time slot are set to 0.
    //   </td>
    //  </tr>
    //  <tr>
    //   <td>UVW</td>
    //   <td>The UVW coordinates in meters as [3,nbaseline].
    //       They are calculated for a missing time slot.
    //   </td>
    //  </tr>
    //  <tr>
    //   <td>FULLRESFLAG</td>
    //   <td>For each baseline the LOFAR_FULL_RES_FLAG column is stored as
    //       a uChar array with shape [orignchan/8, ntimeavg]. The bits
    //       represent the flags. They are converted to a Bool array with shape
    //       [orignchan, ntimeavg, nbaseline].
    //       Note that orignchan is given by the NCHANNELS keyword in the
    //       LOFAR_FULL_RES_FLAG column.
    //       If column LOFAR_FULL_RES_FLAG is not present, the flags are used
    //       and it is assumed that no averaging was done yet (thus ntimeavg=1
    //       and orignchan=nchan).
    //   </td>
    //  </tr>
    // </table>

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

      // Read the FullRes flags (LOFAR_FULL_RES_FLAG) at the given row numbers.
      // It returns a 3-dim array [norigchan, ntimeavg, nbaseline].
      // If undefined, an empty array is returned.
      virtual casa::Cube<bool> getFullResFlags (const casa::RefRows& rowNrs);

      // Read the given data column at the given row numbers.
      virtual casa::Cube<casa::Complex> getData (const casa::String& columnName,
                                                 const casa::RefRows& rowNrs);

      // Write the flags at the given row numbers.
      // It is used by MSUpdater.
      void putFlags (const casa::RefRows& rowNrs,
                     const casa::Cube<bool>& flags);

      // Get the main MS table.
      casa::Table& table()
        { return itsMS; }

      // Get the rownrs for meta info of missing time slots.
      // It uses the rows of the first time slot.
      const casa::Vector<uint> getBaseRowNrs() const
        { return itsBaseRowNrs; }

      // Get the name of the MS.
      const std::string& msName() const
        { return itsMS.tableName(); }

      // Get the time information.
      double startTime() const
        { return itsFirstTime; }
      double endTime() const
        { return itsLastTime; }
      double timeInterval() const
        { return itsInterval; }

      // Get the nr of averaged full resolution channels.
      uint nchanAvg() const
        { return itsFullResNChanAvg; }
      // Get the nr of averaged full resolution time slots.
      uint ntimeAvg() const
        { return itsFullResNTimeAvg; }

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
      casa::String        itsDataColName;
      bool                itsHasWeightSpectrum;
      bool                itsUseFlags;
      uint                itsStartChan;
      double              itsInterval;
      double              itsFirstTime;
      double              itsLastTime;
      double              itsNextTime;
      uint                itsNrInserted;    //# nr of inserted time slots
      casa::Slicer        itsSlicer;        //# slice in corr,chan
      bool                itsHasFullResFlags;
      uint                itsFullResNChanAvg;
      uint                itsFullResNTimeAvg;
      casa::Slicer        itsFullResSlicer; //# slice in chan,timeavg,bl
      DPBuffer            itsBuffer;
      UVWCalculator       itsUVWCalc;
      casa::Vector<uint>  itsBaseRowNrs;    //# rownrs for meta of missing times
      bool                itsCountFlags;
      FlagCounter         itsFlagCounter;
    };

  } //# end namespace
}

#endif
