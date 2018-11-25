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

#include "DPInput.h"
#include "DPBuffer.h"
#include "UVWCalculator.h"
#include "FlagCounter.h"

#include <casacore/tables/Tables/TableIter.h>
#include <casacore/tables/Tables/RefRows.h>
#include <casacore/casa/Arrays/Slicer.h>

namespace DP3 {

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
    //  <li> msin.autoweight: calculate weights from autocorrelations? [no]
    //  <li> msin.startchan: first channel to use [0]
    //  <li> msin.nchan: number of channels to use [all]
    //  <li> msin.useflag: use the existing flags? [yes]
    //  <li> msin.datacolumn: the data column to use [DATA]
    //  <li> msin.weightcolumn: the weights column to use [WEIGHT_SPECTRUM or
    //           WEIGHT]
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
    //       If autoweight is on, the autocorrelations are used to
    //       calculate proper weights.
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
    //       If column LOFAR_FULL_RES_FLAG is not present, the flags are used
    //       and it is assumed that no averaging was done yet (thus ntimeavg=1
    //       and orignchan=nchan).
    //   </td>
    //  </tr>
    // </table>

    class MSReader: public DPInput
    {
    public:
      // Default constructor.
      MSReader();

      // Construct the object for the given MS.
      // Parameters are obtained from the parset using the given prefix.
      // The missingData argument is for MultiMSReader.
      MSReader (const std::string& msName,
                const ParameterSet&, const string& prefix,
                bool missingData = false);

      virtual ~MSReader();

      // Process the next data chunk.
      // It returns false when at the end.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Add some data to the MeasurementSet written/updated.
      // Do nothing.
      virtual void addToMS (const string&) {};

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // If needed, show the flag counts.
      virtual void showCounts (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

      // Read the UVW at the given row numbers into the buffer.
      virtual void getUVW (const casacore::RefRows& rowNrs,
                           double time,
                           DPBuffer&);

      // Read the weights at the given row numbers into the buffer.
      // Note: the buffer must contain DATA if autoweighting is in effect.
      virtual void getWeights (const casacore::RefRows& rowNrs,
                               DPBuffer&);

      // Read the fullRes flags (LOFAR_FULL_RES_FLAG) at the given row numbers
      // into the buffer.
      // If there is no such column, the flags are set to false and false is
      // returned.
      virtual bool getFullResFlags (const casacore::RefRows& rowNrs,
                                    DPBuffer&);

      // Read the model data at the given row numbers into the array.
      virtual void getModelData (const casacore::RefRows& rowNrs,
                                 casacore::Cube<casacore::Complex>&);

#ifdef HAVE_LOFAR_BEAM
      // Fill the vector with station beam info from the input MS.
      // Only fill it for the given station names.
      virtual void fillBeamInfo (std::vector<LOFAR::StationResponse::Station::Ptr>&,
                                 const casacore::Vector<casacore::String>& antNames);
#endif
      
      // Tell if the visibility data are to be read.
      virtual void setReadVisData (bool readVisData);

      // Get the main MS table.
      casacore::Table& table()
        { return itsMS; }

      // Get the name of the data column to be used.
      const std::string& dataColumnName() const
        { return itsDataColName; }

      const std::string& weightColumnName() const
        { return itsWeightColName; }

      const std::string& modelColumnName() const
        { return itsModelColName; }

      // Get the slicer in the FLAG and DATA column.
      const casacore::Slicer& colSlicer() const
        { return itsColSlicer; }

      // Get the rownrs for meta info of missing time slots.
      // It uses the rows of the first time slot.
      const casacore::Vector<uint>& getBaseRowNrs() const
        { return itsBaseRowNrs; }

      // Get the name of the MS.
      virtual std::string msName() const;

      // Get the time information.
      double firstTime() const
        { return itsFirstTime; }
      double lastTime() const
        { return itsLastTime; }

      // Get the selected spectral window.
      uint spectralWindow() const
        { return itsSpw; }

      // Get the baseline selection.
      const string& baselineSelection() const
        { return itsSelBL; }

      // Is the data column missing?
      bool missingData() const
        { return itsMissingData; }

      // Get the start channel.
      uint startChan() const
        { return itsStartChan; }

      // Get the nr of averaged full resolution channels.
      uint nchanAvgFullRes() const
        { return itsFullResNChanAvg; }
      // Get the nr of averaged full resolution time slots.
      uint ntimeAvgFullRes() const
        { return itsFullResNTimeAvg; }

      // Tell if the input MS has LOFAR_FULL_RES_FLAG.
      bool hasFullResFlags() const
        { return itsHasFullResFlags; }

      // Get access to the buffer.
      const DPBuffer& getBuffer() const
        { return itsBuffer; }

      // Flags inf and NaN
      static void flagInfNaN(const casacore::Cube<casacore::Complex>& dataCube,
                       casacore::Cube<bool>& flagsCube, FlagCounter& flagCounter);

    private:
      // Prepare the access to the MS.
      // Return the first and last time and the interval.
      void prepare (double& firstTime, double& lastTime,
                    double& interval);

      // Do the rest of the preparation.
      void prepare2();

      // Skip the first times in the MS in case a start time was given.
      // If needed, it sets itsFirstTime properly.
      void skipFirstTimes();

      // Calculate the UVWs for a missing time slot.
      void calcUVW (double time, DPBuffer&);

      // Calculate the weights from the autocorrelations.
      void autoWeight (casacore::Cube<float>& weights, const DPBuffer& buf);

    protected:
      //# Data members.
      std::string        itsMSName;
      casacore::Table         itsMS;
      casacore::Table         itsSelMS;         //# possible selection of spw, baseline
      casacore::TableIterator itsIter;
      std::string        itsDataColName;
      std::string        itsWeightColName;
      std::string        itsModelColName;
      std::string        itsStartChanStr;  //# startchan expression
      std::string        itsNrChanStr;     //# nchan expression
      std::string              itsSelBL;         //# Baseline selection string
      bool                itsReadVisData;   //# read visibility data?
      bool                itsNeedSort;      //# sort needed on time,baseline?
      bool                itsAutoWeight;    //# calculate weights from autocorr?
      bool                itsAutoWeightForce; //# always calculate weights?
      bool                itsHasWeightSpectrum;
      bool                itsUseFlags;
      bool                itsUseAllChan;    //# all channels (i.e. no slicer)?
      bool                itsMissingData;   //# allow missing data column?
      int                 itsSpw;           //# spw (band) to use (<0 no select)
      uint                itsNrBl;
      uint                itsNrCorr;
      uint                itsNrChan;
      uint                itsStartChan;
      double              itsTimeTolerance; //# tolerance for time comparison
      double              itsTimeInterval;
      double              itsStartTime;
      double              itsFirstTime;
      double              itsLastTime;
      double              itsNextTime;
      double              itsLastMSTime;
      uint                itsNrRead;        //# nr of time slots read from MS
      uint                itsNrInserted;    //# nr of inserted time slots
      casacore::Slicer        itsColSlicer;     //# slice in corr,chan column
      casacore::Slicer        itsArrSlicer;     //# slice in corr,chan,bl array
      bool                itsHasFullResFlags;
      uint                itsFullResNChanAvg;
      uint                itsFullResNTimeAvg;
      DPBuffer            itsBuffer;
      UVWCalculator       itsUVWCalc;
      casacore::Vector<uint>  itsBaseRowNrs;    //# rownrs for meta of missing times
      FlagCounter         itsFlagCounter;
      NSTimer             itsTimer;
    };

  } //# end namespace
}

#endif
