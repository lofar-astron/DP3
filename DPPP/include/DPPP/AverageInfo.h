//# AverageInfo.h: Info how the data are averaged in time and frequency
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

#ifndef DPPP_AVERAGEINFO_H
#define DPPP_AVERAGEINFO_H

// @file
// @brief Info how the data are averaged in time and frequency.

#include <Common/LofarTypes.h>

namespace LOFAR {
  namespace DPPP {

    // @ingroup NDPPP

    // This class contains the information about the number of correlations,
    // channels, baselines, and times.
    // It is initialized by the first step and updated by steps like
    // Averager that change the number of channels or times.
    // Steps can take information from it to know about shapes.

    class AverageInfo
    {
    public:
      // Default constructor.
      AverageInfo();

      // Set the initial info from the input.
      void init (uint ncorr, uint startChan, uint nchan, uint nbaselines,
                 uint ntime, double timeInterval);

      // Update the info from the given average factors.
      // If chanAvg is higher than the actual nr of channels, it is reset.
      // It returns the possibly reset nr of channels to average.
      uint update (uint chanAvg, uint timeAvg);

      // Get the info.
      uint ncorr() const
        { return itsNCorr; }
      uint startChan() const
        { return itsStartChan; }
      uint origNChan() const
        { return itsOrigNChan; }
      uint nchan() const
        { return itsNChan; }
      uint nchanAvg() const
        { return itsChanAvg; }
      uint nbaselines() const
        { return itsNBl; }
      uint ntime() const
        { return itsNTime; }
      uint ntimeAvg() const
        { return itsTimeAvg; }
      double timeInterval() const
        { return itsTimeInterval; }

    private:
      uint   itsNCorr;
      uint   itsStartChan;
      uint   itsOrigNChan;
      uint   itsNChan;
      uint   itsChanAvg;
      uint   itsNBl;
      uint   itsNTime;
      uint   itsTimeAvg;
      double itsTimeInterval;
    };

  } //# end namespace
}

#endif
