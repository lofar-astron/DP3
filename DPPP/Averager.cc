//# Averager.cc: DPPP step class to average in time and/or freq
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

#include "Averager.h"
#include "DPBuffer.h"
#include "DPInfo.h"

#include "../Common/ParallelFor.h"
#include "../Common/ParameterSet.h"
#include "../Common/StringUtil.h"

#include <casacore/casa/Arrays/ArrayMath.h>

#include <cassert>
#include <iostream>
#include <iomanip>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    Averager::Averager (DPInput* input,
                        const ParameterSet& parset,
                        const string& prefix)
      : itsInput     (input),
        itsName      (prefix),
        itsMinNPoint (parset.getUint  (prefix+"minpoints", 1)),
        itsMinPerc   (parset.getFloat (prefix+"minperc", 0.) / 100.),
        itsNTimes    (0),
        itsTimeInterval (0),
        itsNoAvg (true)
    {
      string freqResolutionStr = parset.getString(prefix+"freqresolution","0");
      itsFreqResolution = getFreqHz(freqResolutionStr);

      if (itsFreqResolution > 0) {
        itsNChanAvg = 0; // Will be set later in updateinfo
      } else {
        itsNChanAvg = parset.getUint  (prefix+"freqstep", 1);
      }

      itsTimeResolution = parset.getFloat(prefix+"timeresolution", 0.);
      if (itsTimeResolution > 0) {
        itsNTimeAvg = 0; // Will be set later in updateInfo
      } else {
        itsNTimeAvg = parset.getUint(prefix+"timestep", 1);
      }
    }

    Averager::Averager (DPInput* input, const string& stepName,
                        uint nchanAvg, uint ntimeAvg)
      : itsInput     (input),
        itsName      (stepName),
        itsFreqResolution (0),
        itsTimeResolution (0),
        itsNChanAvg  (nchanAvg),
        itsNTimeAvg  (ntimeAvg),
        itsMinNPoint (1),
        itsMinPerc   (0),
        itsNTimes    (0),
        itsTimeInterval (0)
    {
      if (itsNChanAvg <= 0) itsNChanAvg = 1;
      if (itsNTimeAvg <= 0) itsNTimeAvg = 1;
      itsNoAvg = (itsNChanAvg == 1  &&  itsNTimeAvg == 1);
    }

    Averager::~Averager()
    {}

    void Averager::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();
      info().setWriteFlags();
      info().setMetaChanged();

      if (itsNChanAvg <= 0) {
        if (itsFreqResolution > 0) {
          double chanwidth = infoIn.chanWidths()[0];
          itsNChanAvg = std::max(1, (int)(itsFreqResolution / chanwidth + 0.5));
        } else {
          itsNChanAvg = 1;
        }
      }

      itsTimeInterval = infoIn.timeInterval();
      if (itsNTimeAvg <= 0) {
        if (itsTimeResolution > 0) {
          itsNTimeAvg = std::max(1, (int)(itsTimeResolution / itsTimeInterval + 0.5));
        } else {
          itsNTimeAvg = 1;
        }
      }

      itsNoAvg = (itsNChanAvg == 1  &&  itsNTimeAvg == 1);

      // Adapt averaging to available nr of channels and times.
      itsNTimeAvg = std::min (itsNTimeAvg, infoIn.ntime());
      itsNChanAvg = info().update (itsNChanAvg, itsNTimeAvg);
    }

    void Averager::show (std::ostream& os) const
    {
      os << "Averager " << itsName << std::endl;
      os << "  freqstep:       " << itsNChanAvg;
      if (itsFreqResolution>0) {
        os << " (set by freqresolution: " << itsFreqResolution << " Hz)" << std::endl;
      }
      os << "  timestep:       " << itsNTimeAvg;
      if (itsTimeResolution>0) {
        os << " (set by timeresolution: " << itsTimeResolution << ")";
      }
      os << std::endl;
      os << "  minpoints:      " << itsMinNPoint << std::endl;
      os << "  minperc:        " << 100*itsMinPerc << std::endl;
    }

    void Averager::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " Averager " << itsName << endl;
    }

    bool Averager::process (const DPBuffer& buf)
    {
      // Nothing needs to be done if no averaging.
      if (itsNoAvg) {
        getNextStep()->process (buf);
        return true;
      }
      itsTimer.start();
      // Sum the data in time applying the weights.
      // The summing in channel and the averaging is done in function average.
      if (itsNTimes == 0) {
        // The first time we assign because that is faster than first clearing
        // and adding thereafter.
        itsBuf.getData().assign (buf.getData());
        itsBuf.getFlags().assign (buf.getFlags());
        itsBuf.getUVW().assign (itsInput->fetchUVW(buf, itsBuf, itsTimer));
        itsBuf.getWeights().assign (itsInput->fetchWeights(buf, itsBuf, itsTimer));
        IPosition shapeIn = buf.getData().shape();
        itsNPoints.resize (shapeIn);
        itsAvgAll.reference (buf.getData() * itsBuf.getWeights());
        itsWeightAll.resize (shapeIn);
        itsWeightAll = itsBuf.getWeights();
        // Take care of the fullRes flags.
        // We have to shape the output array and copy to a part of it.
        const Cube<bool>& fullResFlags =
          itsInput->fetchFullResFlags (buf, itsBufTmp, itsTimer);
        IPosition ofShape = fullResFlags.shape();
        ofShape[1] *= itsNTimeAvg;      // more time entries, same chan and bl
        itsBuf.getFullResFlags().resize (ofShape);
        itsBuf.getFullResFlags() = true; // initialize for times missing at end
        copyFullResFlags (fullResFlags, buf.getFlags(), 0);
        // Set middle of new interval.
        double time = buf.getTime() + 0.5*(itsNTimeAvg-1)*itsTimeInterval;
        itsBuf.setTime     (time);
        itsBuf.setExposure (itsNTimeAvg*itsTimeInterval);
        // Only set.
        itsNPoints = 1;
        // Set flagged points to zero.
        Array<bool>::const_contiter infIter = buf.getFlags().cbegin();
        Array<Complex>::contiter    dataIter = itsBuf.getData().cbegin();
        Array<float>::contiter      wghtIter = itsBuf.getWeights().cbegin();
        Array<int>::contiter        outnIter = itsNPoints.cbegin();
        Array<int>::contiter        outnIterEnd = itsNPoints.cend();
        while (outnIter != outnIterEnd) {
          if (*infIter) {
            // Flagged data point
            *outnIter = 0;
            *dataIter = Complex();
            *wghtIter = 0;
          } else {
            // Weigh the data point
            *dataIter *= *wghtIter;
          }
          ++infIter;
          ++dataIter;
          ++wghtIter;
          ++outnIter;
        }
      } else {
        // Not the first time.
        // For now we assume that all timeslots have the same nr of baselines,
        // so check if the buffer sizes are the same.
        assert (itsBuf.getData().shape() == buf.getData().shape());
        itsBufTmp.referenceFilled (buf);
        itsBuf.getUVW() += itsInput->fetchUVW (buf, itsBufTmp, itsTimer);
        copyFullResFlags (itsInput->fetchFullResFlags (buf, itsBufTmp, itsTimer),
                          buf.getFlags(), itsNTimes);
        const Cube<float>& weights =
          itsInput->fetchWeights (buf, itsBufTmp, itsTimer);
        // Ignore flagged points.
        Array<Complex>::const_contiter indIter = buf.getData().cbegin();
        Array<float>::const_contiter   inwIter = weights.cbegin();
        Array<bool>::const_contiter    infIter = buf.getFlags().cbegin();
        Array<Complex>::contiter outdIter = itsBuf.getData().cbegin();
        Array<Complex>::contiter alldIter = itsAvgAll.cbegin();
        Array<float>::contiter   outwIter = itsBuf.getWeights().cbegin();
        Array<float>::contiter   allwIter = itsWeightAll.cbegin();
        Array<int>::contiter outnIter    = itsNPoints.cbegin();
        Array<int>::contiter outnIterEnd = itsNPoints.cend();
        while (outnIter != outnIterEnd) {
          *alldIter += *indIter * *inwIter;
          *allwIter += *inwIter;
          if (!*infIter) {
            *outdIter += *indIter * *inwIter;
            *outwIter += *inwIter;
            (*outnIter)++;
          }
          ++indIter;
          ++inwIter;
          ++infIter;
          ++outdIter;
          ++alldIter;
          ++outwIter;
          ++allwIter;
          ++outnIter;
        }
      }
      // Do the averaging if enough time steps have been processed.
      itsNTimes += 1;
      if (itsNTimes >= itsNTimeAvg) {
        average();
        itsTimer.stop();
        getNextStep()->process (itsBufOut);
        itsNTimes = 0;
      } else {
        itsTimer.stop();
      }
      return true;
    }

    void Averager::finish()
    {
      // Average remaining entries.
      if (itsNTimes > 0) {
        itsTimer.start();
        average();
        itsTimer.stop();
        getNextStep()->process (itsBufOut);
        itsNTimes = 0;
      }
      // Let the next steps finish.
      getNextStep()->finish();
    }

    void Averager::average()
    {
      IPosition shp = itsBuf.getData().shape();
      uint nchanin = shp[1];
      uint npin = shp[0] * nchanin;
      shp[1] = (shp[1] + itsNChanAvg - 1) / itsNChanAvg;
      itsBufOut.getData().resize (shp);
      itsBufOut.getWeights().resize (shp);
      itsBufOut.getFlags().resize (shp);
      uint ncorr = shp[0];
      uint nchan = shp[1];
      uint  nbl   = shp[2];
      uint npout = ncorr * nchan;
      ParallelFor<uint> loop(getInfo().nThreads());
      loop.Run(0, nbl, [&](uint k, size_t /*thread*/) {
        const Complex* indata = itsBuf.getData().data() + k*npin;
        const Complex* inalld = itsAvgAll.data() + k*npin;
        const float* inwght = itsBuf.getWeights().data() + k*npin;
        const float* inallw = itsWeightAll.data() + k*npin;
        const int* innp = itsNPoints.data() + k*npin;
        Complex* outdata = itsBufOut.getData().data() + k*npout;
        float* outwght = itsBufOut.getWeights().data() + k*npout;
        bool* outflags = itsBufOut.getFlags().data() + k*npout;
        for (uint i=0; i<ncorr; ++i) {
          uint inxi = i;
          uint inxo = i;
          for (uint ch=0; ch<nchan; ++ch) {
            uint nch = std::min(itsNChanAvg, nchanin - ch*itsNChanAvg);
            uint navgAll = nch * itsNTimes;
            Complex sumd;
            Complex sumad;
            float   sumw  = 0;
            float   sumaw = 0;
            uint    np = 0;
            for (uint j=0; j<nch; ++j) {
              sumd  += indata[inxi]; // Note: weight is accounted for in process
              sumad += inalld[inxi];
              sumw  += inwght[inxi];
              sumaw += inallw[inxi];
              np    += innp[inxi];
              inxi  += ncorr;
            }
            // Flag the point if insufficient unflagged data.
            if (sumw == 0  ||  np < itsMinNPoint  || np < navgAll*itsMinPerc) {
              outdata[inxo]  = (sumaw==0  ?  Complex() : sumad/sumaw);
              outflags[inxo] = true;
              outwght[inxo]  = sumaw;
            } else {
              outdata[inxo]  = sumd / sumw;
              outflags[inxo] = false;
              outwght[inxo]  = sumw;
            }
            inxo += ncorr;
          }
        }
      });
      // Set the remaining values in the output buffer.
      itsBufOut.setTime     (itsBuf.getTime());
      itsBufOut.setExposure (itsBuf.getExposure());
      itsBufOut.setFullResFlags (itsBuf.getFullResFlags());
      // The result UVWs are the average of the input.
      // If ever needed, UVWCalculator can be used to calculate the UVWs.
      itsBufOut.setUVW (itsBuf.getUVW() / double(itsNTimes));
    }

    void Averager::copyFullResFlags (const Cube<bool>& fullResFlags,
                                     const Cube<bool>& flags,
                                     int timeIndex)
    {
      // Copy the fullRes flags to the given index.
      // Furthermore the appropriate FullRes flags are set for a
      // flagged data point. It can be the case that an input data point
      // has been averaged before, thus has fewer channels than FullResFlags.
      // nchan and nbl are the same for in and out.
      // ntimout is a multiple of ntimavg.
      IPosition shapeIn  = fullResFlags.shape();
      IPosition shapeOut = itsBuf.getFullResFlags().shape();
      IPosition shapeFlg = flags.shape();
      uint nchan    = shapeIn[0];    // original nr of channels
      uint ntimavg  = shapeIn[1];    // nr of averaged times in input data
      uint nchanavg = nchan / shapeFlg[1]; // nr of avg chan in input data
      uint ntimout  = shapeOut[1];   // nr of averaged times in output data
      uint nbl      = shapeIn[2];    // nr of baselines
      uint ncorr    = shapeFlg[0];   // nr of correlations (in FLAG)
      // in has to be copied to the correct time index in out.
      bool* outBase = itsBuf.getFullResFlags().data() + nchan*ntimavg*timeIndex;
      for (uint k=0; k<nbl; ++k) {
        const bool* inPtr   = fullResFlags.data() + k*nchan*ntimavg;
        const bool* flagPtr = flags.data() + k*ncorr*shapeFlg[1];
        bool* outPtr = outBase + k*nchan*ntimout;
        memcpy (outPtr, inPtr, nchan*ntimavg*sizeof(bool));
        // Applying the flags only needs to be done if the input data
        // was already averaged before.
        if (ntimavg > 1  ||  nchanavg > 1) {
          for (int j=0; j<shapeFlg[1]; ++j) {
            // If a data point is flagged, the flags in the corresponding
            // FullRes window have to be set.
            // Only look at the flags of the first correlation.
            if (*flagPtr) {
              bool* avgPtr = outPtr + j*nchanavg;
              for (uint i=0; i<ntimavg; ++i) {
                std::fill (avgPtr, avgPtr+nchanavg, true);
                avgPtr += nchan;
              }
            }
            flagPtr += ncorr;
          }
        }
      }
    }

    double Averager::getFreqHz(const string& freqstr) {
      String unit;
      // See if a unit is given at the end.
      String v(freqstr);
      // Remove possible trailing blanks.
      rtrim(v);
      Regex regex("[a-zA-Z]+$");
      string::size_type pos = v.index (regex);
      if (pos != String::npos) {
        unit = v.from   (pos);
        v    = v.before (pos);
      }
      // Set value and unit.

      double value = strToDouble(v);
      if (unit.empty()) {
        return value;
      } else {
        Quantity q(value, unit);
        return q.getValue("Hz", true);
      }
    }

  } //# end namespace
}
