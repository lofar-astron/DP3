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
      info.update (itsNChanAvg, itsNTimeAvg);
    }

    void Averager::show (std::ostream& os)
    {
      os << "Averager " << itsName << std::endl;
      os << "  freqstep        " << itsNChanAvg << std::endl;
      os << "  timestep        " << itsNTimeAvg << std::endl;
    }

    bool Averager::process (const DPBuffer& buf)
    {
      RefRows rowNrs(buf.getRowNrs());
      // If first time, resize and clear the buffer.
      if (itsNTimes == 0) {
        // The first time we assign because that is faster than first clearing
        // and adding thereafter.
        itsBuf.getWeights() = itsInput->fetchWeights (buf, rowNrs);
        itsBuf.getData()    = buf.getData();
        IPosition shapeIn = buf.getData().shape();
        itsNPoints.resize (shapeIn);
        // Take care of the preAvg flags.
        // If not averaging in time, we can simply use them from the input.
        // Note that this is possible because we do not keep nChanAvg in the
        // preAvgFlags cube and because AverageInfo.cc checks if averaging
        // in channel fits integrally.
        if (itsNTimeAvg == 1) {
          itsBuf.setPreAvgFlags (itsInput->fetchPreAvgFlags (buf, rowNrs));
        } else {
          // We have to shape the output array and copy to a part of it.
          Array<bool> preAvgFlags(itsInput->fetchPreAvgFlags (buf, rowNrs));
          IPosition ofShape = preAvgFlags.shape();
          ofShape[1] *= itsNTimeAvg;      // more time entries, same chan and bl
          // Make it unique in case PreAvg is referenced elsewhere.
          itsBuf.getPreAvgFlags().unique();
          itsBuf.getPreAvgFlags().resize (ofShape);
          itsBuf.getPreAvgFlags() = true; // initialize for times missing at end
          copyPreAvgFlags (preAvgFlags, 0);
        }
        // Set middle of new interval.
        double time = buf.getTime() + 0.5*(itsNTimeAvg-1)*itsTimeInterval;
        itsBuf.setTime (time);
        // Only set.
        itsNPoints = 1;
        if (! buf.hasNoFlags()) {
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
              // Weight the data point
              *dataIter *= *wghtIter;
            }
            ++infIter;
            ++dataIter;
            ++wghtIter;
            ++outnIter;
          }
        }
      } else {
        // Not the first time.
        // For now we assume that all timeslots have the same nr of baselines,
        // so check if the buffer sizes are the same.
        ASSERT (itsBuf.getData().shape() == buf.getData().shape());
        copyPreAvgFlags (itsInput->fetchPreAvgFlags(buf, rowNrs), itsNTimes);
        Cube<float> weights(itsInput->fetchWeights(buf, rowNrs));
        if (buf.hasNoFlags()) {
          // No flags, so we can simply add.
          itsBuf.getData()    += buf.getData() * weights;
          itsBuf.getWeights() += weights;
          itsNPoints += 1;
        } else {
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
      if (itsNTimes > 0) {
        DPBuffer buf = average();
        getNextStep()->process (buf);
        itsNTimes = 0;
      }
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
      buf.setTime (itsBuf.getTime());
      buf.setPreAvgFlags (itsBuf.getPreAvgFlags());
      DPBuffer::mergePreAvgFlags (buf.getPreAvgFlags(), itsBuf.getFlags());
      return buf;
    }

    void Averager::copyPreAvgFlags (const Cube<bool>& preAvgFlags,
                                    int timeIndex)
    {
      // Copy the preAvg flags.
      // nrchan and nrbl are the same for in and out.
      // nrtimout is a multiple of nrtimin.
      IPosition shapeIn = preAvgFlags.shape();
      IPosition shapeOut = itsBuf.getPreAvgFlags().shape();
      uint nrchan = shapeIn[0];
      uint nrtimin= shapeIn[1];
      uint nrtimout = shapeOut[1];
      uint nrbl = shapeIn[2];
      // in has to be copied to the correct time index in out.
      bool* outPtr = itsBuf.getPreAvgFlags().data() + nrchan*nrtimin*timeIndex;
      const bool* inPtr = preAvgFlags.data();
      for (uint i=0; i<nrbl; ++i) {
        memcpy (outPtr, inPtr, nrchan*nrtimin*sizeof(bool));
        outPtr += nrchan*nrtimout;
        inPtr  += nrchan*nrtimin;
      }
    }

  } //# end namespace
}
