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

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <Common/lofar_vector.h>

namespace LOFAR {
  class ParameterSet;
  class ParameterValue;

  namespace DPPP {

    // @ingroup NDPPP

    // This class is a DPStep class flagging data points based on data
    // selections given in the parset file.
    // The following selections can be given:
    // <ul>
    //  <li> minimum and/or maximum UV distance
    //  <li> autocorrelations
    //  <li> baselines using names for antenna 1 and 2
    //  <li> antennae using antenna names
    //  <li> minimum and/or maximum amplitude per correlation
    //  <li> channel numbers
    //  <li> frequency ranges
    //  <li> sequence nr or time ranges
    // </ul>
    // The antenna names can contain shell-style wildcards (* ? [] {}).

    class PreFlagger: public DPStep
    {
    public:
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
      virtual void updateAverageInfo (AverageInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      class PSet
      {
      public:
        // Define the shared pointer for this type.
        typedef shared_ptr<PSet> ShPtr;

        // Construct from the parset parameters.
        PSet (DPInput*, const ParameterSet& parset, const string& prefix);

        // Set and return the flags.
        const casa::Cube<bool>& process (DPBuffer&,
                                         const casa::Block<bool>& matchBL);

        // Update the average info.
        // It is used to adjust the parms if needed.
        void updateInfo (const AverageInfo&);

        // Show the pset parameters.
        void show (std::ostream&) const;

      private:
        // Read and handle the Time parset parameters.
        void readTimeParms (const ParameterSet& parset);

        // Clear matchBL for mismatching baselines.
        void flagBL (const casa::Vector<int>& ant1,
                     const casa::Vector<int>& ant2);

        // Clear matchBL for baselines with mismatching UV distances.
        void flagUV (const casa::Matrix<double>& uvw);

        // Clear matchBL for baselines with mismatching AzEl.
        void flagAzEl ();

        // Set the flags based on amplitude threshold per correlation.
        void flagAmpl (const casa::Cube<float>& amplitudes);

        // Set the flags based on phase threshold per correlation.
        void flagPhase (const casa::Cube<casa::Complex>& data);

        // Set the flags based on real/imaginary threshold per correlation.
        void flagComplex (const casa::Cube<casa::Complex>& data);

        // Flag the channels given in itsChannels.
        void flagChannels();

        // Fill the baseline matrix; set true for baselines to flag.
        void fillBLMatrix (const casa::Vector<casa::String>& antNames);

        // Return a vector with a value per correlation.
        // If no parm value given use the default.
        // If the parm value is a vector, use the given values. Use the
        // default for non-given correlations.
        // If the parm value is a single value, use it for all correlations.
        // <br>doFlag is set if values are given.
        vector<float> fillValuePerCorr (const ParameterValue& value,
                                        float defVal, bool& doFlag);

        // Handle the frequency ranges given and determine which channels
        // have to be flagged.
        void handleFreqRanges (const casa::Vector<double>& chanFreqs);

        // Get the value and possible unit.
        // If no unit is given, the argument is left untouched.
        void getValue (const string& str, double& value, casa::String& unit);

        // Get the frequency in Hz using the value and unit.
        double getFreqHz (double value, const casa::String& unit);

        //# Data members of PreFlagger::PSet.
        DPInput*           itsInput;
        string             itsName;
        bool               itsFlagOnUV; //# true = do uv distance based flagging
        bool               itsFlagOnBL; //# true = do ant/bl based flagging
        bool               itsFlagOnAmpl; //# true = do amplitude based flagging
        bool               itsFlagOnPhase;//# true = do phase based flagging
        bool               itsFlagOnRI;   //# true = do real/imag based flagging
        bool               itsFlagOnAzEl; //# true = do Az/El based flagging
        bool               itsFlagOnLST;  //# true = do LST based flagging
        double             itsMinUV;    //# minimum UV distance; <0 means ignore
        double             itsMaxUV;    //# maximum UV distance; <0 means ignore
        casa::Matrix<bool> itsFlagBL;   //# true = flag baseline [i,j]
        vector<double>     itsElevation;//# elevation ranges to be flagged
        vector<double>     itsAzimuth;  //# azimuth ranges to be flagged
        vector<double>     itsTimes;    //# time ranges to be flagged
        vector<double>     itsLST;      //# sidereal time ranges to be flagged
        vector<float>      itsAmplMin;  //# minimum amplitude for each corr
        vector<float>      itsAmplMax;  //# maximum amplitude for each corr
        vector<float>      itsPhaseMin; //# minimum phase for each corr
        vector<float>      itsPhaseMax; //# maximum phase for each corr
        vector<float>      itsRealMin;  //# minimum real for each corr
        vector<float>      itsRealMax;  //# maximum real for each corr
        vector<float>      itsImagMin;  //# minimum imaginary for each corr
        vector<float>      itsImagMax;  //# maximum imaginary for each corr
        vector<uint>       itsChannels; //# channels to be flagged.
        vector<uint>       itsFlagChan; //# channels given to be flagged.
        vector<string>     itsFlagFreq; //# frequency ranges to be flagged
        string             itsCorrType; //# auto, cross, or all
        string             itsStrBL;    //# the baseline string
        vector<PSet::ShPtr> itsPSets;
        casa::Matrix<bool>  itsChanFlags; //# flags for channels to be flagged
        casa::Cube<bool>    itsFlags;
        casa::Block<bool>   itsMatchBL; //# true = baseline in buffer matches 
      };
        
      //# Data members of PreFlagger.
      string  itsName;
      NSTimer itsTimer;
      PSet    itsPSet;
    };
      
  } //# end namespace
}

#endif
