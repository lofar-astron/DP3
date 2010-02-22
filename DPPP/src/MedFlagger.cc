//# MedFlagger.h: DPPP step class to average in time and/or freq
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

#include <lofar_config.h>
#include <DPPP/MedFlagger.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/AverageInfo.h>
#include <Common/ParameterSet.h>
#include <Common/StreamUtil.h>
#include <Common/LofarLogger.h>
#include <casa/Arrays/ArrayMath.h>
#include <iostream>
#include <algorithm>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MedFlagger::MedFlagger (const ParameterSet& parset, const string& prefix)
      : itsName       (prefix),
        itsThreshold  (parset.getUint (prefix+"threshold", 1)),
        itsFreqWindow (parset.getUint (prefix+"freqwindow", 1)),
        itsTimeWindow (parset.getUint (prefix+"timewindow", 1)),
        itsNTimes     (0),
        itsNTimesDone (0),
        itsFlagCounter("MedFlagger " + prefix)
    {
      ASSERT (itsFreqWindow > 0  &&  itsTimeWindow > 0);
      itsBuf.resize (itsTimeWindow);
      itsFlagCorr = parset.getUintVector (prefix+"correlations",
                                          vector<uint>());
    }

    MedFlagger::~MedFlagger()
    {}

    void MedFlagger::show (std::ostream& os) const
    {
      os << "MADFlagger " << itsName << std::endl;
      os << "  freqwindow:     " << itsFreqWindow << std::endl;
      os << "  timewindow:     " << itsTimeWindow << std::endl;
      os << "  threshold:      " << itsThreshold << std::endl;
      os << "  flagcorr:       " << itsFlagCorr << std::endl;
    }

    void MedFlagger::showCounts (std::ostream& os) const
    {
    }

    void MedFlagger::updateAverageInfo (AverageInfo& info)
    {
      // Check if frequency window is not too large.
      // Check if window sizes are odd.
      itsFreqWindow = std::min(itsFreqWindow, info.nchan());
      if (itsFreqWindow%2 == 0) {
        itsFreqWindow -= 1;
      }
      if (itsTimeWindow%2 == 0) {
        itsTimeWindow -= 1;
      }
      // Set or check the correlations to flag on.
      vector<uint> flagCorr;
      uint ncorr = info.ncorr();
      if (itsFlagCorr.empty()) {
        // No correlations given means use them all.
        for (uint i=0; i<ncorr; ++i) {
          flagCorr.push_back (i);
        }
      } else {
        for (uint i=0; i<itsFlagCorr.size(); ++i) {
          // Only take valid corrrelations.
          if (itsFlagCorr[i] < ncorr) {
            flagCorr.push_back (itsFlagCorr[i]);
          }
        }
        // If no valid left, use first one.
        if (flagCorr.empty()) {
          LOG_INFO_STR ("No valid correlations given in MedFlagger "
                        << itsName << "; first one will be used");
          flagCorr.push_back (0);
        }
      }
      itsFlagCorr = flagCorr;
      // Initialize the flag counters.
      itsFlagCounter.init (info.nbaselines(), info.nchan(), info.ncorr());
    }

    bool MedFlagger::process (const DPBuffer& buf)
    {
      // Accumulate in the time window.
      // The buffer is wrapped, thus oldest entries are overwritten.
      uint index = itsNTimes % itsTimeWindow;
      itsBuf[index] = buf;
      // Calculate amplitudes if needed.
      DPBuffer& dbuf = itsBuf[index];
      if (dbuf.getAmplitudes().empty()) {
        dbuf.setAmplitudes (amplitude(dbuf.getData()));
      }
      // Fill flags if needed.
      if (dbuf.getFlags().empty()) {
        dbuf.getFlags().resize (dbuf.getData().shape());
        dbuf.getFlags() = false;
      }
      itsNTimes++;
      ///cout << "medproc: " << itsNTimes << endl;
      // Flag if there are enough time entries in the buffer.
      if (itsNTimes > itsTimeWindow/2) {
        // Fill the vector telling which time entries to use for the medians.
        // Usually it is simple and all time entires are needed.
        vector<uint> timeEntries;
        timeEntries.reserve (itsTimeWindow);
        for (uint i=0; i<itsTimeWindow; ++i) {
          timeEntries.push_back (i);
        }
        // If window not entirely full, use copies as needed.
        // This is done as follows:
        // Suppose timewindow=7 and we only have entries 0,1,2,3,4.
        // The entries are mirrored, thus we get 4,3,2,1,0,1,2,3,4
        // to obtain sufficient time entries.
        for (uint i=itsNTimes; i<itsTimeWindow; ++i) {
          timeEntries[i] = i-itsNTimes+1;
        }
        flag (itsNTimesDone%itsTimeWindow, timeEntries);
        itsNTimesDone++;
      }
      return true;
    }

    void MedFlagger::finish()
    {
      // Adjust window size if there are fewer time entries.
      if (itsNTimes < itsTimeWindow) {
        itsTimeWindow = itsNTimes;
      }
      uint halfWindow = itsTimeWindow/2;
      vector<uint> timeEntries(itsTimeWindow);
      // Process possible leading entries.
      // This can happen if the window was larger than number of times.
      while (itsNTimesDone <= halfWindow) {
        // Process in the same way as in process where nt replaces itsNTimes.
        uint nt = itsNTimesDone + halfWindow/2 + 1;
        for (uint i=0; i<nt; ++i) {
          timeEntries[i] = i;
        }
        // Mirror as needed.
        for (uint i=nt; i<itsTimeWindow; ++i) {
          timeEntries[i] = i-nt+1;
        }
        flag (itsNTimesDone, timeEntries);
        itsNTimesDone++;
      }
      ASSERT (itsNTimes - itsNTimesDone == halfWindow);
      // Process the remaining time entries.
      // Time entries have to be mirrored.
      while (itsNTimesDone < itsNTimes) {
        for (uint i=0; i<itsTimeWindow; ++i) {
          uint t = itsNTimesDone + i - halfWindow;
          if (t >= itsNTimes) {
            t = itsNTimes + itsNTimes - t - 2;
          }
          ///cout << t << ' ';
          timeEntries[i] = t%itsTimeWindow;
        }
        ///cout << endl;
        flag (itsNTimesDone%itsTimeWindow, timeEntries);
        itsNTimesDone++;
      }
      // Let the next step finish its processing.
      getNextStep()->finish();
    }

    void MedFlagger::flag (uint index, const vector<uint>& timeEntries)
    {
      ///cout << "flag: " <<itsNTimes<<' '<<itsNTimesDone<<' ' <<index << timeEntries << endl;
      // Result is 'copy' of the entry at the given time index.
      DPBuffer buf (itsBuf[index]);
      IPosition shp = buf.getData().shape();
      uint ncorr = shp[0];
      uint nchan = shp[1];
      uint nrbl  = shp[2];
      uint ntime = timeEntries.size();
      // Create a temporary buffer to hold data for determining medians.
      Block<float> tempBuf(itsFreqWindow*ntime);
      // Get pointers to data and flags.
      const float* dataPtr = buf.getAmplitudes().data();
      bool* flagPtr = buf.getFlags().data();
      float Z1, Z2;
      float MAD = 1.4826;   //# constant determined by Pandey
      // Now flag each baseline, channel and correlation for this time window.
      for (uint ib=0; ib<nrbl; ++ib) {
        for (uint ic=0; ic<nchan; ++ic) {
          bool corrIsFlagged = false;
          // Iterate over given correlations.
          for (vector<uint>::const_iterator iter = itsFlagCorr.begin();
               iter != itsFlagCorr.end(); ++iter) {
            uint ip = *iter;
            // If one correlation is flagged, all of them will be flagged.
            // So no need to check others.
            if (flagPtr[ip]) {
              corrIsFlagged = true;
              break;
            }
            // Calculate values from the median.
            computeFactors (timeEntries, ib, ic, ip, nchan, ncorr,
                            Z1, Z2, tempBuf.storage());
            if (dataPtr[ip] > Z1 + itsThreshold * Z2 * MAD) {
              corrIsFlagged = true;
              itsFlagCounter.incrBaseline(ib);
              itsFlagCounter.incrChannel(ic);
              itsFlagCounter.incrCorrelation(ip);
              break;
            }
          }
          if (corrIsFlagged) {
            for (uint ip=0; ip<ncorr; ++ip) {
              flagPtr[ip] = true;
            }
          }
          dataPtr += ncorr;
          flagPtr += ncorr;
        }
      }
      // Process the result in the next step.
      getNextStep()->process (buf);
    }
            
    void MedFlagger::computeFactors (const vector<uint>& timeEntries,
                                     uint bl, int chan, int corr,
                                     int nchan, int ncorr,
                                     float& Z1, float& Z2,
                                     float* tempBuf)
    {
      // Collect all non-flagged data points for given baseline, channel,
      // and correlation in the window around the channel.
      uint np = 0;
      // At the beginning or end of the window the values are wrapped.
      // So we might need to move in two parts.
      // This little piece of code is tested in tMirror.cc.
      int hw = itsFreqWindow/2;
      int s1 = chan - hw;
      int e1 = chan + hw + 1;
      int s2 = 1;
      int e2 = 1;
      if (s1 < 0) {
        e2 = -s1 + 1;
        s1 = 0;
      } else if (e1 > nchan) {
        s2 = nchan + nchan - e1 - 1; // e1-nchan+1 too far, go back that amount
        e2 = nchan-1;
        e1 = nchan;
      }
      // Iterate over all time entries.
      for (vector<uint>::const_iterator iter=timeEntries.begin();
           iter != timeEntries.end(); ++iter) {
        const DPBuffer& inbuf = itsBuf[*iter];
        // Get pointers to given baseline and correlation.
        uint offset = bl*nchan*ncorr + corr;
        const float* dataPtr = inbuf.getAmplitudes().data() + offset;
        const bool*  flagPtr = inbuf.getFlags().data() + offset;
        // Now move data from the two channel parts.
        for (int i=s1*ncorr; i<e1*ncorr; i+=ncorr) {
          if (!flagPtr[i]) {
            tempBuf[np++] = dataPtr[i];
          }
        }
        for (int i=s2*ncorr; i<e2*ncorr; i+=ncorr) {
          if (!flagPtr[i]) {
            tempBuf[np++] = dataPtr[i];
          }
        }
      }
      // If only flagged data, don't do anything.
      if (np == 0) {
        Z1 = -1.0;
        Z2 = 0.0;
      } else {
        // Get median of data and get median of absolute difference.
        ///std::nth_element (tempBuf, tempBuf+np/2, tempBuf+np);
        ///Z1 = *(tempBuf+np/2);
        Z1 = GenSort<float>::kthLargest (tempBuf, np, np/2);
        for (uint i=0; i<np; ++i) {
          tempBuf[i] = std::abs(tempBuf[i] - Z1);
        }
        ///std::nth_element (tempBuf, tempBuf+np/2, tempBuf+np);
        ///Z2 = *(tempBuf+np/2);
        Z2 = GenSort<float>::kthLargest (tempBuf, np, np/2);
      }
    }


  } //# end namespace
}
