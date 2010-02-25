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

#include <lofar_config.h>
#include <DPPP/Averager.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/AverageInfo.h>
#include <Common/ParameterSet.h>
#include <Common/LofarLogger.h>
#include <casa/Arrays/ArrayMath.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    Averager::Averager (DPInput* input,
                        const ParameterSet& parset, const string& prefix)
      : itsInput    (input),
        itsName     (prefix),
        itsNChanAvg (parset.getUint (prefix+"freqstep", 1)),
        itsNTimeAvg (parset.getUint (prefix+"timestep", 1)),
        itsNTimes   (0),
        itsTimeInterval (0)
    {
      ASSERTSTR (itsNChanAvg > 1  ||  itsNTimeAvg > 1,
                 "freqstep and/or timestep has to be specified when averaging");
    }

    Averager::~Averager()
    {}

    void Averager::updateAverageInfo (AverageInfo& info)
    {
      itsTimeInterval = info.timeInterval();
      itsNChanAvg = info.update (itsNChanAvg, itsNTimeAvg);
    }

    void Averager::show (std::ostream& os) const
    {
      os << "Averager " << itsName << std::endl;
      os << "  freqstep:       " << itsNChanAvg << std::endl;
      os << "  timestep:       " << itsNTimeAvg << std::endl;
    }

    bool Averager::process (const DPBuffer& buf)
    {
      RefRows rowNrs(buf.getRowNrs());
      if (itsNTimes == 0) {
        // The first time we assign because that is faster than first clearing
        // and adding thereafter.
        itsBuf.getUVW()     = itsInput->fetchUVW (buf, rowNrs);
        itsBuf.getWeights() = itsInput->fetchWeights (buf, rowNrs);
        itsBuf.getData()    = buf.getData();
        IPosition shapeIn   = buf.getData().shape();
        itsNPoints.resize (shapeIn);
        // Take care of the fullRes flags.
        // We have to shape the output array and copy to a part of it.
        Cube<bool> fullResFlags(itsInput->fetchFullResFlags (buf, rowNrs));
        IPosition ofShape = fullResFlags.shape();
        ofShape[1] *= itsNTimeAvg;      // more time entries, same chan and bl
        // Make it unique in case FullRes is referenced elsewhere.
        // (itsBuf.FullRes is referenced when set in buf by average())
        itsBuf.getFullResFlags().unique();
        itsBuf.getFullResFlags().resize (ofShape);
        itsBuf.getFullResFlags() = true; // initialize for times missing at end
        copyFullResFlags (fullResFlags, buf.getFlags(), 0);
        // Set middle of new interval.
        double time = buf.getTime() + 0.5*(itsNTimeAvg-1)*itsTimeInterval;
        itsBuf.setTime (time);
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
        ASSERT (itsBuf.getData().shape() == buf.getData().shape());
        itsBuf.getUVW() += itsInput->fetchUVW (buf, rowNrs);
        copyFullResFlags (itsInput->fetchFullResFlags(buf, rowNrs),
                          buf.getFlags(), itsNTimes);
        Cube<float> weights(itsInput->fetchWeights(buf, rowNrs));
        // Ignore flagged points.
        Array<Complex>::const_contiter indIter = buf.getData().cbegin();
        Array<float>::const_contiter   inwIter = weights.cbegin();
        Array<bool>::const_contiter    infIter = buf.getFlags().cbegin();
        Array<Complex>::contiter outdIter = itsBuf.getData().cbegin();
        Array<float>::contiter   outwIter = itsBuf.getWeights().cbegin();
        Array<int>::contiter outnIter = itsNPoints.cbegin();
        Array<int>::contiter outnIterEnd = itsNPoints.cend();
        while (outnIter != outnIterEnd) {
          if (!*infIter) {
            *outdIter += *indIter * *inwIter;
            *outwIter += *inwIter;
            (*outnIter)++;
          }
          ++indIter;
          ++inwIter;
          ++infIter;
          ++outdIter;
          ++outwIter;
          ++outnIter;
        }
      }
      // Do the averaging if enough time steps have been processed.
      itsNTimes += 1;
      if (itsNTimes >= itsNTimeAvg) {
        DPBuffer buf = average();
        getNextStep()->process (buf);
        itsNTimes = 0;
      }
      return true;
    }

    void Averager::finish()
    {
      // Average remaining entries.
      if (itsNTimes > 0) {
        DPBuffer buf = average();
        getNextStep()->process (buf);
        itsNTimes = 0;
      }
      // Let the next steps finish.
      getNextStep()->finish();
    }

    DPBuffer Averager::average() const
    {
      IPosition shp = itsBuf.getData().shape();
      uint nchanin = shp[1];
      uint npin = shp[0] * nchanin;
      shp[1] = (shp[1] + itsNChanAvg - 1) / itsNChanAvg;
      DPBuffer buf;
      buf.getData().resize (shp);
      buf.getWeights().resize (shp);
      buf.getFlags().resize (shp);
      uint ncorr = shp[0];
      uint nchan = shp[1];
      uint nbl   = shp[2];
      uint npout = ncorr * nchan;
      const Complex* indata = itsBuf.getData().data();
      const float* inwght = itsBuf.getWeights().data();
      const int* innp = itsNPoints.data();
      Complex* outdata = buf.getData().data();
      float* outwght = buf.getWeights().data();
      bool* outflags = buf.getFlags().data();
      for (uint k=0; k<nbl; ++k) {
        for (uint i=0; i<ncorr; ++i) {
          uint inxi = i;
          uint inxo = i;
          for (uint ch=0; ch<nchan; ++ch) {
            uint nch = std::min(itsNChanAvg, nchanin - ch*itsNChanAvg);
            Complex sumd;
            float   sumw=0;
            uint np = 0;
            for (uint j=0; j<nch; ++j) {
              sumd += indata[inxi];
              sumw += inwght[inxi];
              np   += innp[inxi];
              inxi += ncorr;
            }
            if (sumw == 0) {
              outdata[inxo] = Complex();
              outflags[inxo] = true;
            } else {
              outdata[inxo] = sumd / sumw;
              outflags[inxo] = false;
            }
            outwght[inxo] = sumw;
            inxo += ncorr;
          }
        }
        // Increment data pointers for the next baseline.
        indata   += npin;
        inwght   += npin;
        innp     += npin;
        outdata  += npout;
        outwght  += npout;
        outflags += npout;
      }
      // Make sure the loops ended correctly.
      DBGASSERT (indata == itsBuf.getData().data() + itsBuf.getData().size());
      DBGASSERT (outdata == buf.getData().data() + buf.getData().size());
      // Set the remaining values in the output buffer.
      buf.setTime (itsBuf.getTime());
      buf.setFullResFlags (itsBuf.getFullResFlags());
      // The result UVWs are the average of the input.
      // If ever needed, UVWCalculator can be used to calculate the UVWs.
      buf.setUVW (itsBuf.getUVW() / double(itsNTimes));
      return buf;
    }

    void Averager::copyFullResFlags (const Cube<bool>& fullResFlags,
                                     const Cube<bool>& flags,
                                     int timeIndex)
    {
      // Copy the fullRes flags to the given index.
      // Furthermore the appropriate FullRes flags are set for a
      // flagged data point. It can be the case that an input data point
      // has been averaged before, thus has fewer channels than FullResFlags.
      // nrchan and nrbl are the same for in and out.
      // nrtimout is a multiple of nrtimavg.
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
      bool* outPtr = itsBuf.getFullResFlags().data() + nchan*ntimavg*timeIndex;
      const bool* inPtr   = fullResFlags.data();
      const bool* flagPtr = flags.data();
      for (uint k=0; k<nbl; ++k) {
        memcpy (outPtr, inPtr, nchan*ntimavg*sizeof(bool));
        // Applying the flags only needs to be done if the input data
        // was already averaged before.
        if (ntimavg > 1  &&  nchanavg > 1) {
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
        outPtr += nchan*ntimout;
        inPtr  += nchan*ntimavg;
      }
    }

  } //# end namespace
}
