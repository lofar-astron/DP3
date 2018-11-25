//# DPInfo.h: General info about DPPP data processing attributes like averaging
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

#ifndef DPPP_DPINFO_H
#define DPPP_DPINFO_H

// @file
// @brief General info about DPPP data processing attributes like averaging

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MeasureHolder.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Containers/Record.h>

namespace DP3 {
  namespace DPPP {

    //# Forward declarations.
    class DPInput;
    
    enum BeamCorrectionMode {
      NoBeamCorrection = 0,
      FullBeamCorrection = 1,
      ArrayFactorBeamCorrection = 2,
      ElementBeamCorrection = 3
    };

    // @ingroup NDPPP

    // This class contains the information about the number of correlations,
    // channels, baselines, and times.
    // It is initialized by the first step and updated by steps like
    // Averager that change the number of channels or times.
    // Steps can take information from it to know about shapes.

    class DPInfo
    {
    public:
      // Default constructor.
      DPInfo();

      // Set the initial info from the input.
      void init (uint ncorr, uint startChan, uint nchan,
                 uint ntime, double startTime, double timeInterval,
                 const string& msName, const string& antennaSet);

      // Set nr of channels.
      void setNChan (uint nchan)
        { itsNChan = nchan; }

      // Set time interval
      void setTimeInterval (double timeInterval)
        { itsTimeInterval = timeInterval; }

      // Set the frequency info.
      // An empty resolutions or effectiveBW is default to chanWidths.
      // If totalBW is 0, it is set to the sum of effectiveBW.
      // If refFreq is 0, it is set to the middle of chanFreqs (mean if even).
      void set (const casacore::Vector<double>& chanFreqs,
                const casacore::Vector<double>& chanWidths,
                const casacore::Vector<double>& resolutions= casacore::Vector<double>(),
                const casacore::Vector<double>& effectiveBW= casacore::Vector<double>(),
                double totalBW = 0,
                double refFreq = 0);

      // Set array info.
      void set (const casacore::MPosition& arrayPos,
                const casacore::MDirection& phaseCenter,
                const casacore::MDirection& delayCenter,
                const casacore::MDirection& tileBeamDir);

      // Set the info for the given antennae and baselines.
      void set (const casacore::Vector<casacore::String>& antNames,
                const casacore::Vector<casacore::Double>& antDiam,
                const std::vector<casacore::MPosition>& antPos,
                const casacore::Vector<casacore::Int>& ant1,
                const casacore::Vector<casacore::Int>& ant2);

      // Set the name of the data column
      void setDataColName(const std::string& dataColName) {
        itsDataColName=dataColName;
      }

      // Set the name of the weight column
      void setWeightColName(const std::string& weightColName) {
        itsWeightColName=weightColName;
      }

      // Update the info for the given average factors.
      // If chanAvg is higher than the actual nr of channels, it is reset.
      // The same is true for timeAvg.
      // It returns the possibly reset nr of channels to average.
      uint update (uint chanAvg, uint timeAvg);

      // Update the info from the given selection parameters.
      // Optionally unused stations are really removed from the antenna lists.
      void update (uint startChan, uint nchan,
                   const std::vector<uint>& baselines, bool remove);

      // Remove unused stations from the antenna lists.
      void removeUnusedAnt();

      // Set the phase center.
      // If original=true, it is set to the original phase center.
      void setPhaseCenter (const casacore::MDirection& phaseCenter, bool original)
        { itsPhaseCenter=phaseCenter; itsPhaseCenterIsOriginal = original; }


      // Get the info.
      const string& msName() const
        { return itsMSName; }
      const string& antennaSet() const
        { return itsAntennaSet; }
      uint ncorr() const
        { return itsNCorr; }
      uint nchan() const
        { return itsNChan; }
      uint startchan() const
        { return itsStartChan; }
      uint origNChan() const
        { return itsOrigNChan; }
      uint nchanAvg() const
        { return itsChanAvg; }
      uint nantenna() const
        { return itsAntNames.size(); }
      uint nbaselines() const
        { return itsAnt1.size(); }
      uint ntime() const
        { return itsNTime; }
      uint ntimeAvg() const
        { return itsTimeAvg; }
      double startTime() const
        { return itsStartTime; }
      double timeInterval() const
        { return itsTimeInterval; }
      const casacore::Vector<casacore::Int>& getAnt1() const
        { return itsAnt1; }
      const casacore::Vector<casacore::Int>& getAnt2() const
        { return itsAnt2; }
      const casacore::Vector<casacore::String>& antennaNames() const
        { return itsAntNames; }
      const casacore::Vector<casacore::Double>& antennaDiam() const
        { return itsAntDiam; }
      const std::vector<casacore::MPosition>& antennaPos() const
        { return itsAntPos; }
      const casacore::MPosition& arrayPos() const
        { return itsArrayPos; }
      const casacore::MPosition arrayPosCopy() const
        {  return copyMeasure(casacore::MeasureHolder(itsArrayPos)).asMPosition(); }
      const casacore::MDirection& phaseCenter() const
        { return itsPhaseCenter; }
      const casacore::MDirection phaseCenterCopy() const
      {  return copyMeasure(casacore::MeasureHolder(itsPhaseCenter)).asMDirection(); }
      bool phaseCenterIsOriginal() const
        { return itsPhaseCenterIsOriginal; }
      const casacore::MDirection& delayCenter() const
        { return itsDelayCenter; }
      const casacore::MDirection delayCenterCopy() const
        { return copyMeasure(casacore::MeasureHolder(itsDelayCenter)).asMDirection(); }
      const casacore::MDirection& tileBeamDir() const
        { return itsTileBeamDir; }
      const casacore::MDirection tileBeamDirCopy() const
        { return copyMeasure(casacore::MeasureHolder(itsTileBeamDir)).asMDirection(); }
      const casacore::Vector<double>& chanFreqs() const
        { return itsChanFreqs; }
      const casacore::Vector<double>& chanWidths() const
        { return itsChanWidths; }
      const casacore::Vector<double>& resolutions() const
        { return itsResolutions; }
      const casacore::Vector<double>& effectiveBW() const
        { return itsEffectiveBW; }
      const std::string& getDataColName() const
        { return itsDataColName; }
      const std::string& getWeightColName() const
        { return itsWeightColName; }
      double totalBW() const
        { return itsTotalBW; }
      double refFreq() const
        { return itsRefFreq; }

      // Get the antenna numbers actually used in the (selected) baselines.
      // E.g. [0,2,5,6]
      const std::vector<int>& antennaUsed() const
        { return itsAntUsed; }

      // Get the indices of all antennae in the used antenna vector above.
      // -1 means that the antenna is not used.
      // E.g. [0,-1,1,-1,-1,2,3] for the example above.
      const std::vector<int>& antennaMap() const
        { return itsAntMap; }

      // Are the visibility data needed?
      bool needVisData() const
        { return itsNeedVisData; }
      // Does the last step need to write data and/or flags?
      bool needWrite() const
        { return itsWriteData || itsWriteFlags || itsWriteWeights; }
      bool writeData() const
        { return itsWriteData; }
      bool writeFlags() const
        { return itsWriteFlags; }
      bool writeWeights() const
        { return itsWriteWeights; }
      // Has the meta data been changed in a step (precluding an update)?
      bool metaChanged() const
        { return itsMetaChanged; }

      // Set if visibility data needs to be read.
      void setNeedVisData()
        { itsNeedVisData = true; }
      // Set if data needs to be written.
      void setWriteData()
        { itsWriteData = true; }
      void setWriteFlags()
        { itsWriteFlags = true; }
      void setWriteWeights()
        { itsWriteWeights = true; }
      // Clear all write flags.
      void clearWrites()
        { itsWriteData = itsWriteFlags = itsWriteWeights = false; }
      // Set change of meta data.
      void setMetaChanged()
        { itsMetaChanged = true; }
      void clearMetaChanged()
        { itsMetaChanged = false; }

      // Get the baseline table index of the autocorrelations.
      // A negative value means there are no autocorrelations for that antenna.
      const std::vector<int>& getAutoCorrIndex() const;

      // Get the lengths of the baselines (in meters).
      const std::vector<double>& getBaselineLengths() const;
      
      void setNThreads(uint nThreads)
      { itsNThreads = nThreads; }
      
      uint nThreads() const
        { return itsNThreads; }

      // Convert to a Record.
      // The names of the fields in the record are the data names without 'its'.
      casacore::Record toRecord() const;

      // Update the DPInfo object from a Record.
      // It is possible that only a few fields are defined in the record.
      void fromRecord (const casacore::Record& rec);
      
      enum BeamCorrectionMode beamCorrectionMode() const
      { return itsBeamCorrectionMode; };
      void setBeamCorrectionMode(enum BeamCorrectionMode mode)
      { itsBeamCorrectionMode = mode; }
      
      const casacore::MDirection& beamCorrectionDir() const
      { return itsBeamCorrectionDir; }
      void setBeamCorrectionDir(const casacore::MDirection& dir)
      { itsBeamCorrectionDir = dir; }

    private:
      // Set which antennae are actually used.
      void setAntUsed();

      // Creates a real copy of a casacore::Measure by exporting to a Record
      static casacore::MeasureHolder copyMeasure(const casacore::MeasureHolder fromMeas);

      //# Data members.
      bool   itsNeedVisData;    //# Are the visibility data needed?
      bool   itsWriteData;      //# Must the data be written?
      bool   itsWriteFlags;     //# Must the flags be written?
      bool   itsWriteWeights;   //# Must the weights be written?
      bool   itsMetaChanged;    //# Are meta data changed? (e.g., by averaging)
      string itsMSName;
      std::string itsDataColName;
      std::string itsWeightColName;
      string itsAntennaSet;
      uint   itsNCorr;
      uint   itsStartChan;
      uint   itsOrigNChan;
      uint   itsNChan;
      uint   itsChanAvg;
      uint   itsNTime;
      uint   itsTimeAvg;
      double itsStartTime;
      double itsTimeInterval;
      casacore::MDirection itsPhaseCenter;
      bool             itsPhaseCenterIsOriginal;
      casacore::MDirection itsDelayCenter;
      casacore::MDirection itsTileBeamDir;
      enum BeamCorrectionMode itsBeamCorrectionMode;
      casacore::MDirection itsBeamCorrectionDir;
      casacore::MPosition  itsArrayPos;
      casacore::Vector<double>       itsChanFreqs;
      casacore::Vector<double>       itsChanWidths;
      casacore::Vector<double>       itsResolutions;
      casacore::Vector<double>       itsEffectiveBW;
      double                     itsTotalBW;
      double                     itsRefFreq;
      casacore::Vector<casacore::String> itsAntNames;
      casacore::Vector<casacore::Double> itsAntDiam;
      std::vector<casacore::MPosition>    itsAntPos;
      std::vector<int>                itsAntUsed;
      std::vector<int>                itsAntMap;
      casacore::Vector<casacore::Int>    itsAnt1;          //# ant1 of all baselines
      casacore::Vector<casacore::Int>    itsAnt2;          //# ant2 of all baselines
      mutable std::vector<double>     itsBLength;       //# baseline lengths
      mutable std::vector<int>        itsAutoCorrIndex; //# autocorr index per ant
      uint itsNThreads;
    };

  } //# end namespace
}

#endif
