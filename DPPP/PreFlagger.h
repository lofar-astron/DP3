//# PreFlagger.h: DPPP step class to flag data on channel, baseline, or time
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

#ifndef DPPP_PREFLAGGER_H
#define DPPP_PREFLAGGER_H

// @file
// @brief DPPP step class to flag data on channel, baseline, or time

#include "DPInput.h"
#include "DPBuffer.h"
#include "BaselineSelection.h"

#include <casacore/measures/Measures/MDirection.h>

namespace DP3 {
  class ParameterSet;
  class ParameterValue;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class flagging data points based on data
    // selections given in the parset file.
    // The following selections can be given:
    // <ul>
    //  <li> minimum and/or maximum UV distance (projected baseline length)
    //  <li> minimum and/or maximum baseline length (intrinsic)
    //  <li> autocorrelations or crosscorrelations
    //  <li> baselines using names for antenna 1 and 2
    //  <li> antennae using antenna names
    //  <li> minimum and/or maximum amplitude/phase/real/imag per correlation
    //  <li> channel numbers
    //  <li> frequency ranges
    //  <li> sequence nr or time ranges
    //  <li> LST
    //  <li> azimuth/elevation
    // </ul>
    // The antenna names can contain shell-style wildcards (* ? [] {}).
    //
    // All selections are ANDed, thus only the data points matching all
    // selections are flagged. It is however, possible to specify a logical
    // expression of selections by means of the internal PSet class.
    // A PSet objects contains a set of ANDed selections. The PSets can
    // be logically combined by the user using the normal logical operators.

    class PreFlagger: public DPStep
    {
    // Make this Test class a friend, so it can access private code.
    friend class TestPSet;

    public:
      enum Mode {SetFlag, ClearFlag, SetComp, ClearComp};

      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      PreFlagger (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~PreFlagger();

      // Process the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the average info.
      // It is used to adjust the parms if needed.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the flag counts.
      virtual void showCounts (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      // This internal class represents a single set of ANDed selections.
      // PSets can be logically combined by the PreFlagger class.
      class PSet
      {
      // Make this Test class a friend, so it can access private code.
      friend class TestPSet;

      public:
        // Define the operators in pset expressions.
        // They have to have negative values in order of precedence.
        enum Oper {
          OpParen = -1,
          OpOr    = -2,
          OpAnd   = -3,
          OpNot   = -4
        };

        // Define the shared pointer for this type.
        typedef std::shared_ptr<PSet> ShPtr;

        // Default constructor (for test purposes).
        PSet()
        {}

        // Construct from the parset parameters.
        PSet (DPInput*, const ParameterSet& parset, const string& prefix);

        // Set and return the flags.
        casacore::Cube<bool>* process (const DPBuffer&, DPBuffer&, uint timeSlot,
                                   const casacore::Block<bool>& matchBL,
                                   NSTimer& timer);

        // Update the general info.
        // It is used to adjust the parms if needed.
        void updateInfo (const DPInfo&);

        // Show the pset parameters.
        void show (std::ostream&, bool showName) const;

      private:
        // Test if the time matches the time ranges.
        bool matchTime (double time, uint timeSlot) const;

        // Test if the value matches one of the ranges in the vector.
        bool matchRange (double v, const std::vector<double>& ranges) const;

        // Clear itsMatchBL for mismatching baselines.
        // If returns false if no matches were found.
        bool flagBL();

        // Clear itsMatchBL for baselines with mismatching UV distances.
        // If returns false if no matches were found.
        bool flagUV (const casacore::Matrix<double>& uvw);

        // Clear itsMatchBL for baselines with mismatching AzEl.
        // If returns false if no matches were found.
        bool flagAzEl (double time);

        // Test if azimuth or elevation of given antenna mismatches.
        // If so, clear itsMatchBL for all baselines containing the antenna.
        void testAzEl (casacore::MDirection::Convert& converter,
                       uint blnr, int ant,
                       const int* ant1, const int* ant2);

        // Set the flags based on amplitude threshold per correlation.
        void flagAmpl (const casacore::Cube<float>& amplitudes);

        // Set the flags based on phase threshold per correlation.
        void flagPhase (const casacore::Cube<casacore::Complex>& data);

        // Set the flags based on real/imaginary threshold per correlation.
        void flagReal (const casacore::Cube<casacore::Complex>& data);
        void flagImag (const casacore::Cube<casacore::Complex>& data);

        // Flag the channels given in itsChannels.
        void flagChannels();

        // Convert a string of (date)time ranges to double. Each range
        // must be given with .. or +-.
        // <tt>asTime=true</tt> means that the strings should contain times,
        // otherwise date/times.
        std::vector<double> fillTimes (const std::vector<string>& str, bool asTime,
                                  bool canEndBeforeStart);

        // Read the string as time or date/time and convert to seconds.
        // usepm indicates if the value is a plusminus value. If so, the
        // value must be a positive time.
        double getSeconds (const string& str, bool asTime, bool usepm);

        // Fill the baseline matrix; set true for baselines to flag.
        void fillBLMatrix();

        // Fill itsChannels if channel/freq selection is done.
        void fillChannels (const DPInfo&);

        // Return a vector with a value per correlation.
        // If no parm value given use the default.
        // If the parm value is a vector, use the given values. Use the
        // default for non-given correlations.
        // If the parm value is a single value, use it for all correlations.
        // <br>doFlag is set if values are given.
        std::vector<float> fillValuePerCorr (const ParameterValue& value,
                                        float defVal, bool& doFlag);

        // Handle the frequency ranges given and determine which channels
        // have to be flagged.
        casacore::Vector<bool> handleFreqRanges
        (const casacore::Vector<double>& chanFreqs);

        // Get the value and possible unit.
        // If no unit is given, the argument is left untouched.
        void getValue (const string& str, double& value, casacore::String& unit);

        // Get the frequency in Hz using the value and unit.
        double getFreqHz (double value, const casacore::String& unit);

        // Convert a PSet expression to Reversed Polish Notation in itsRpn.
        // It returns the names of all PSets.
        std::vector<string> exprToRpn (const string& expr);

        //# Data members of PreFlagger::PSet.
        DPInput*           itsInput;
        const DPInfo*      itsInfo;
        string             itsName;
        string             itsStrExpr;
        bool               itsFlagOnTimeOnly; //# true = only flag on time info
        bool               itsFlagOnTime; //# true = do time based flagging
        bool               itsFlagOnUV; //# true = do uv distance based flagging
        bool               itsFlagOnBL; //# true = do ant/bl based flagging
        bool               itsFlagOnAmpl; //# true = do amplitude based flagging
        bool               itsFlagOnPhase;//# true = do phase based flagging
        bool               itsFlagOnReal; //# true = do real based flagging
        bool               itsFlagOnImag; //# true = do imag based flagging
        bool               itsFlagOnAzEl; //# true = do Az/El based flagging
        BaselineSelection  itsSelBL;
        double             itsMinUV;    //# minimum UV distance; <0 means ignore
        double             itsMaxUV;    //# maximum UV distance; <0 means ignore
        casacore::Matrix<bool> itsFlagBL;   //# true = flag baseline [i,j]
        std::vector<double>     itsAzimuth;  //# azimuth ranges to be flagged
        std::vector<double>     itsElevation;//# elevation ranges to be flagged
        std::vector<double>     itsTimes;    //# time of day ranges to be flagged
        std::vector<double>     itsLST;      //# sidereal time ranges to be flagged
        std::vector<double>     itsATimes;   //# absolute time ranges to be flagged
        std::vector<double>     itsRTimes;   //# relative time ranges to be flagged
        std::vector<uint>       itsTimeSlot; //# time slots to be flagged
        std::vector<float>      itsAmplMin;  //# minimum amplitude for each corr
        std::vector<float>      itsAmplMax;  //# maximum amplitude for each corr
        std::vector<float>      itsPhaseMin; //# minimum phase for each corr
        std::vector<float>      itsPhaseMax; //# maximum phase for each corr
        std::vector<float>      itsRealMin;  //# minimum real for each corr
        std::vector<float>      itsRealMax;  //# maximum real for each corr
        std::vector<float>      itsImagMin;  //# minimum imaginary for each corr
        std::vector<float>      itsImagMax;  //# maximum imaginary for each corr
        std::vector<uint>       itsChannels; //# channels to be flagged.
        std::vector<string>     itsStrChan;  //# channel ranges to be flagged.
        std::vector<string>     itsStrFreq;  //# frequency ranges to be flagged
        std::vector<string>     itsStrTime;  //# time ranges to be flagged
        std::vector<string>     itsStrLST;   //# LST ranges to be flagged
        std::vector<string>     itsStrATime; //# absolute time ranges to be flagged
        std::vector<string>     itsStrRTime; //# relative time ranges to be flagged
        std::vector<string>     itsStrAzim;  //# azimuth ranges to be flagged
        std::vector<string>     itsStrElev;  //# elevation ranges to be flagged
        std::vector<int>         itsRpn;     //# PSet expression in RPN form
        std::vector<PSet::ShPtr> itsPSets;   //# PSets used in itsRpn
        casacore::Matrix<bool>  itsChanFlags; //# flags for channels to be flagged
        casacore::Cube<bool>    itsFlags;
        casacore::Block<bool>   itsMatchBL; //# true = baseline in buffer matches 
      };

      // Set the flags in outPtr where inPtr matches mode.
      void setFlags (const bool* inPtr, bool* outPtr,
                     uint nrcorr, uint nrchan, uint nrbl, bool mode);

      // Clear the flags in outPtr where inPtr matches mode.
      // If the corresponding data point of a flag is invalid
      // (non-finite or zero), it is always flagged.
      void clearFlags (const bool* inPtr, bool* outPtr,
                       uint nrcorr, uint nrchan, uint nrbl, bool mode,
                       const DPBuffer& buf);

      //# Data members of PreFlagger.
      string      itsName;
      DPInput*    itsInput;
      DPBuffer    itsBuffer;
      Mode        itsMode;
      NSTimer     itsTimer;
      PSet        itsPSet;
      uint        itsCount;
      FlagCounter itsFlagCounter;
    };
      
  } //# end namespace
}

#endif
