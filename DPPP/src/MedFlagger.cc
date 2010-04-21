//# MedFlagger.cc: DPPP step class to flag data based on median filtering
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
#include <casa/Containers/Record.h>
#include <casa/Containers/RecordField.h>
#include <tables/Tables/ExprNode.h>
#include <tables/Tables/RecordGram.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCPosition.h>
#include <iostream>
#include <algorithm>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MedFlagger::MedFlagger (DPInput* input,
                            const ParameterSet& parset, const string& prefix)

      : itsInput         (input),
        itsName          (prefix),
        itsThresholdStr  (parset.getString (prefix+"threshold", "1")),
        itsFreqWindowStr (parset.getString (prefix+"freqwindow", "1")),
        itsTimeWindowStr (parset.getString (prefix+"timewindow", "1")),
        itsNTimes        (0),
        itsNTimesDone    (0)
    {
      itsFlagCorr = parset.getUintVector (prefix+"correlations",
                                          vector<uint>());
      itsApplyAutoCorr = parset.getBool  (prefix+"applyautocorr", false);
      itsMinBLength    = parset.getDouble(prefix+"blmin", -1);
      itsMaxBLength    = parset.getDouble(prefix+"blmax", 1e30);
      // Determine baseline indices of autocorrelations.
      const Vector<int>& ant1 = itsInput->getAnt1();
      const Vector<int>& ant2 = itsInput->getAnt2();
      int nant = 1 + std::max (max(ant1), max(ant2));
      itsAutoCorrIndex.resize (nant);
      std::fill (itsAutoCorrIndex.begin(), itsAutoCorrIndex.end(), -1);
      itsNrAutoCorr = 0;
      for (uint i=0; i<ant1.size(); ++i) {
        if (ant1[i] == ant2[i]) {
          itsAutoCorrIndex[ant1[i]] = i;
          itsNrAutoCorr++;
        }
      }
      if (itsApplyAutoCorr) {
        ASSERTSTR (itsNrAutoCorr > 0, "applyautocorr=True cannot be used if "
                   "the data does not contain autocorrelations");
      }
      // Calculate the baseline lengths.
      // First get the antenna positions.
      const vector<MPosition>& antPos = itsInput->antennaPos();
      vector<Vector<double> > antVec;
      antVec.reserve (antPos.size());
      for (vector<MPosition>::const_iterator iter = antPos.begin();
           iter != antPos.end(); ++iter) {
        // Convert to ITRF and keep as x,y,z in m.
        antVec.push_back
          (MPosition::Convert(*iter, MPosition::ITRF)().getValue().getValue());
      }
      // Fill in the length of each baseline.
      itsBLength.reserve (ant1.size());
      for (uint i=0; i<ant1.size(); ++i) {
        Array<double> diff(antVec[ant2[i]] - antVec[ant1[i]]);
        itsBLength.push_back (sqrt(sum(diff*diff)));
      }
    }

    MedFlagger::~MedFlagger()
    {}

    void MedFlagger::show (std::ostream& os) const
    {
      os << "MADFlagger " << itsName << std::endl;
      os << "  freqwindow:     " << itsFreqWindowStr
         << "   (max = " << itsFreqWindow << ')' << std::endl;
      os << "  timewindow:     " << itsTimeWindowStr
         << "   (max = " << itsTimeWindow << ')' << std::endl;
      os << "  threshold:      " << itsThresholdStr
         << "   (max = " << itsThreshold << ')' << std::endl;
      os << "  correlations:   " << itsFlagCorr << std::endl;
      os << "  applyautocorr:  " << itsApplyAutoCorr
         << "   (nautocorr = " << itsNrAutoCorr << ')' << std::endl;
      os << "  blmin:          " << itsMinBLength << " m" << std::endl;
      os << "  blmax:          " << itsMaxBLength << " m" << std::endl;
    }

    void MedFlagger::showCounts (std::ostream& os) const
    {
      os << endl << "Flags set by MADFlagger " << itsName;
      os << endl << "=======================" << endl;
      itsFlagCounter.showBaseline (os, itsInput->getAnt1(),
                                   itsInput->getAnt2(), itsNTimes);
      itsFlagCounter.showChannel  (os, itsNTimes);
      itsFlagCounter.showCorrelation (os, itsNTimes);
    }

    void MedFlagger::showTimings (std::ostream& os, double duration) const
    {
      double flagDur = itsTimer.getElapsed();
      os << "  ";
      FlagCounter::showPerc1 (os, flagDur, duration);
      os << " MADFlagger " << itsName << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsMoveTimer.getElapsed(), flagDur);
      os << " of it spent in shuffling data" << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsMedianTimer.getElapsed(), flagDur);
      os << " of it spent in calculating medians" << endl;
    }

    void MedFlagger::updateAverageInfo (AverageInfo& info)
    {
      // Evaluate the window size expressions.
      getExprValues (info.nchan(), info.ntime());
      itsBuf.resize (itsTimeWindow);
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
      itsTimer.start();
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
        // Arrange indices such that any width can be used; thus first center,
        // then one left and right, second left and right, etc.
        // If window not entirely full, use copies as needed.
        // This is done as follows:
        // Suppose timewindow=9 and we only have entries 0,1,2,3,4.
        // The entries are mirrored, thus we get 4,3,2,1,0,1,2,3,4
        // to obtain sufficient time entries.
        vector<uint> timeEntries;
        timeEntries.reserve (itsTimeWindow);
        ///uint rinx = itsNTimesDone % itsTimeWindow;
        ///timeEntries.push_back (rinx);   // center
        ///uint linx = rinx;
        ///for (uint i=1; i<=itsTimeWindow/2; ++i) {
        ///if (linx == 0) linx = itsTimeWindow;
        ///linx--;
        ///rinx++;
        ///if (rinx == itsTimeWindow) rinx = 0;
        ///if (i >= itsNTimes) rinx = linx;
        ///timeEntries.push_back (linx);
        ///timeEntries.push_back (rinx);
        timeEntries.push_back (itsNTimesDone % itsTimeWindow);   // center
        for (uint i=1; i<=itsTimeWindow/2; ++i) {
          timeEntries.push_back
            (std::abs(int(itsNTimesDone) - int(i)) % itsTimeWindow);
          timeEntries.push_back ((itsNTimesDone + i) % itsTimeWindow);
        }
        flag (itsNTimesDone%itsTimeWindow, timeEntries);
        itsNTimesDone++;
      }
      itsTimer.stop();
      return true;
    }

    void MedFlagger::finish()
    {
      itsTimer.start();
      // Adjust window size if there are fewer time entries.
      if (itsNTimes < itsTimeWindow) {
        itsTimeWindow = itsNTimes;
      }
      uint halfWindow = itsTimeWindow/2;
      vector<uint> timeEntries(itsTimeWindow);
      // Process possible leading entries.
      // This can happen if the window was larger than number of times.
      while (itsNTimesDone <= halfWindow) {
        // Process in the same way as in process.
        uint inx = 0;
        timeEntries[inx++] = itsNTimesDone % itsTimeWindow;   // center
        for (uint i=1; i<=halfWindow; ++i) {
          timeEntries[inx++] =
            std::abs(int(itsNTimesDone) - int(i)) % itsTimeWindow;
          timeEntries[inx++] = (itsNTimesDone + i) % itsTimeWindow;
        }
        flag (itsNTimesDone, timeEntries);
        itsNTimesDone++;
      }
      ASSERT (itsNTimes - itsNTimesDone == halfWindow);
      // Process the remaining time entries.
      while (itsNTimesDone < itsNTimes) {
        uint inx = 0;
        timeEntries[inx++] = itsNTimesDone % itsTimeWindow;   // center
        for (uint i=1; i<=halfWindow; ++i) {
          timeEntries[inx++] =
            std::abs(int(itsNTimesDone) - int(i)) % itsTimeWindow;
          // Time entries might need to be mirrored at the end.
          uint ri = itsNTimesDone + i;
          if (ri >= itsNTimes) {
            ri = 2*(itsNTimes-1) - ri;
          }
          timeEntries[inx++] = ri % itsTimeWindow;
        }
        flag (itsNTimesDone%itsTimeWindow, timeEntries);
        itsNTimesDone++;
      }
      itsTimer.stop();
      // Let the next step finish its processing.
      getNextStep()->finish();
    }

    void MedFlagger::flag (uint index, const vector<uint>& timeEntries)
    {
      ///cout << "flag: " <<itsNTimes<<' '<<itsNTimesDone<<' ' <<index << timeEntries << endl;
      // Get antenna numbers in case applyautocorr is true.
      const Vector<int>& ant1 = itsInput->getAnt1();
      const Vector<int>& ant2 = itsInput->getAnt2();
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
        double threshold = itsThresholdArr[ib];
        // Do only autocorrelations if told so.
        // Otherwise do baseline only if length within min-max.
        if ((!itsApplyAutoCorr  &&  itsBLength[ib] >= itsMinBLength  &&
             itsBLength[ib] <= itsMaxBLength)  ||
            (itsApplyAutoCorr  &&  ant1[ib] == ant2[ib])) {
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
              if (dataPtr[ip] > Z1 + threshold * Z2 * MAD) {
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
        } else {
          dataPtr += nchan*ncorr;
          flagPtr += nchan*ncorr;
        }
      }
      // Apply autocorrelations flags if needed.
      // Only to baselines with length within min-max.
      if (itsApplyAutoCorr) {
        flagPtr = buf.getFlags().data();
        for (uint ib=0; ib<nrbl; ++ib) {
          // Flag crosscorr if at least one autocorr is present.
          // Only if baseline length within min-max.
          if (ant1[ib] != ant2[ib]  &&  itsBLength[ib] >= itsMinBLength  &&
              itsBLength[ib] <= itsMaxBLength) {
            int inx1 = itsAutoCorrIndex[ant1[ib]];
            int inx2 = itsAutoCorrIndex[ant2[ib]];
            if (inx1 >= 0  ||  inx2 >= 0) {
              // Find flags of the autocorrelations of both antennae.
              // Use other autocorr if one autocorr does not exist.
              // In this way inx does not need to be tested in the inner loop.
              if (inx1 < 0) {
                inx1 = inx2;
              } else if (inx2 < 0) {
                inx2 = inx1;
              }
              bool* flagAnt1 = buf.getFlags().data() + inx1*nchan*ncorr;
              bool* flagAnt2 = buf.getFlags().data() + inx2*nchan*ncorr;
              // Flag if not flagged yet and if one of autocorr is flagged.
              for (uint ic=0; ic<nchan; ++ic) {
                if (!*flagPtr  &&  (*flagAnt1 || *flagAnt2)) {
                  for (uint ip=0; ip<ncorr; ++ip) {
                    flagPtr[ip] = true;
                  }
                  itsFlagCounter.incrBaseline(ib);
                  itsFlagCounter.incrChannel(ic);
                }
                flagPtr  += ncorr;
                flagAnt1 += ncorr;
                flagAnt2 += ncorr;
              }
            }
          } else {
            flagPtr += nchan*ncorr;
          }
        }
      }
      // Process the result in the next step.
      itsTimer.stop();
      getNextStep()->process (buf);
      itsTimer.start();
    }
            
    void MedFlagger::computeFactors (const vector<uint>& timeEntries,
                                     uint bl, int chan, int corr,
                                     int nchan, int ncorr,
                                     float& Z1, float& Z2,
                                     float* tempBuf)
    {
      itsMoveTimer.start();
      // Collect all non-flagged data points for given baseline, channel,
      // and correlation in the window around the channel.
      uint np = 0;
      // At the beginning or end of the window the values are wrapped.
      // So we might need to move in two parts.
      // This little piece of code is tested in tMirror.cc.
      int hw = itsFreqWindowArr[bl]/2;
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
      const uint* iter = &(timeEntries[0]);
      const uint* endIter = iter + itsTimeWindowArr[bl];
      for (; iter!=endIter; ++iter) {
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
      itsMoveTimer.stop();
      // If only flagged data, don't do anything.
      if (np == 0) {
        Z1 = -1.0;
        Z2 = 0.0;
      } else {
        itsMedianTimer.start();
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
        itsMedianTimer.stop();
      }
    }

    void MedFlagger::getExprValues (int maxNChan, int maxNTime)
    {
      // Parse the expressions.
      // Baseline length can be used as 'bl' in the expressions.
      Record rec;
      rec.define ("bl", double(0));
      TableExprNode node1 (RecordGram::parse(rec, itsFreqWindowStr));
      TableExprNode node2 (RecordGram::parse(rec, itsTimeWindowStr));
      TableExprNode node3 (RecordGram::parse(rec, itsThresholdStr));
      // Size the arrays.
      uInt nrbl = itsBLength.size();
      itsThresholdArr.reserve  (nrbl);
      itsTimeWindowArr.reserve (nrbl);
      itsFreqWindowArr.reserve (nrbl);
      itsFreqWindow = 0;
      itsTimeWindow = 0;
      itsThreshold  = -1e30;
      // Evaluate the expression for each baseline.
      double result;
      RecordFieldPtr<double> blref(rec, "bl");
      for (uint i=0; i<nrbl; ++i) {
        // Put the length of each baseline in the record used to evaluate.
        *blref = itsBLength[i];
        // Evaluate freqwindow size and make it odd if needed.
        node1.get (rec, result);
        int freqWindow = std::min (std::max(1, int(result+0.5)), maxNChan);
        if (freqWindow%2 == 0) {
          freqWindow--;
        }
        itsFreqWindowArr.push_back (freqWindow);
        itsFreqWindow = std::max(itsFreqWindow, uint(freqWindow));
        // Evaluate timewindow size and make it odd if needed.
        node2.get (rec, result);
        int timeWindow = std::max(1, int(result+0.5));
        if (maxNTime > 0  &&  timeWindow > maxNTime) {
          timeWindow = maxNTime;
        }
        if (timeWindow%2 == 0) {
          timeWindow--;
        }
        itsTimeWindowArr.push_back (timeWindow);
        itsTimeWindow = std::max (itsTimeWindow, uint(timeWindow));
        // Evaluate threshold.
        node3.get (rec, result);
        itsThresholdArr.push_back (result);
        if (result > itsThreshold) {
          itsThreshold = result;
        }
      }
    }

  } //# end namespace
}
