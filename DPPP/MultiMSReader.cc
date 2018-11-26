//# MultiMSReader.cc: DPPP step reading from multiple MSs
//# Copyright (C) 2011
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
//# $Id: MSReader.cc 19257 2011-11-14 14:36:08Z diepen $
//#
//# @author Ger van Diepen

#include "MultiMSReader.h"
#include "DPBuffer.h"
#include "DPLogger.h"
#include "DPInfo.h"
#include "Exceptions.h"

#include "../Common/ParameterSet.h"
#include "../Common/StreamUtil.h"

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/Utilities/GenSort.h>
#include <casacore/casa/OS/Conversion.h>

#include <iostream>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    MultiMSReader::MultiMSReader (const vector<string>& msNames,
                                  const ParameterSet& parset,
                                  const string& prefix)
      : itsFirst    (-1),
        itsNMissing (0),
        itsMSNames  (msNames),
        itsRegularChannels (true)
    {
      if(msNames.size() <= 0)
        throw Exception("No names of MeasurementSets given");
      itsMSName           = itsMSNames[0];
      itsStartChanStr     = parset.getString (prefix+"startchan", "0");
      itsNrChanStr        = parset.getString (prefix+"nchan", "0");
      itsUseFlags         = parset.getBool   (prefix+"useflag", true);
      itsDataColName      = parset.getString (prefix+"datacolumn", "DATA");
      itsWeightColName    = parset.getString (prefix+"weightcolumn",
                                                "WEIGHT_SPECTRUM"),
      itsMissingData      = parset.getBool   (prefix+"missingdata", false);
      itsAutoWeight       = parset.getBool   (prefix+"autoweight", false);
      itsNeedSort         = parset.getBool   (prefix+"sort", false);
      itsOrderMS          = parset.getBool   (prefix+"orderms", true);
      // Open all MSs.
      DPStep::ShPtr nullStep (new NullStep());
      itsReaders.reserve (msNames.size());
      itsSteps.reserve   (msNames.size());
      for (uint i=0; i<msNames.size(); ++i) {
        itsReaders.push_back (new MSReader (msNames[i], parset, prefix,
                                            itsMissingData));
        // itsSteps takes care of deletion of the MSReader object.
        itsSteps.push_back   (DPStep::ShPtr(itsReaders[i]));
        // Add a null step for the reader.
        itsSteps[i]->setNextStep (nullStep);
        // Ignore if the MS is missing.
        if (itsReaders[i]->table().isNull()) {
          itsReaders[i] = 0;
          itsNMissing++;
        } else if (itsFirst < 0) {
          itsFirst = i;
        }
      }

      // TODO: check if frequencies are regular, insert some empy readers
      // if necessary

      if(itsFirst<0)
        throw Exception("All input MeasurementSets do not exist");
      itsBuffers.resize (itsReaders.size());
    }

    MultiMSReader::~MultiMSReader()
    {}

    void MultiMSReader::setReadVisData (bool readVisData)
    {
      itsReadVisData = readVisData;
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          itsReaders[i]->setReadVisData (readVisData);
        }
      }
    }

    void MultiMSReader::handleBands()
    {
      if (itsNMissing > 0) {
        fillBands();
        return;
      }
      if (itsOrderMS) {
        sortBands();
      }

      // Collect the channel info of all MSs.
      Vector<double> chanFreqs  (itsNrChan);
      Vector<double> chanWidths (itsNrChan);
      Vector<double> resolutions(itsNrChan);
      Vector<double> effectiveBW(itsNrChan);
      uint inx = 0;
      for (uint i=0; i<itsReaders.size(); ++i) {
        uint nchan = itsReaders[i]->getInfo().nchan();
        objcopy (chanFreqs.data()  + inx,
                 itsReaders[i]->getInfo().chanFreqs().data(),  nchan);
        objcopy (chanWidths.data() + inx,
                 itsReaders[i]->getInfo().chanWidths().data(), nchan);
        objcopy (resolutions.data() + inx,
                 itsReaders[i]->getInfo().resolutions().data(), nchan);
        objcopy (effectiveBW.data() + inx,
                 itsReaders[i]->getInfo().effectiveBW().data(), nchan);
        inx += nchan;
      }
      info().set (chanFreqs, chanWidths, resolutions, effectiveBW, 0., 0.);
    }

    void MultiMSReader::sortBands()
    {
      // Order the bands (MSs) in order of frequency.
      int nband = itsReaders.size();
      Vector<double> freqs(nband);
      for (int i=0; i<nband; ++i) {
        freqs[i] = itsReaders[i]->getInfo().chanFreqs().data()[0];
      }
      Vector<uInt> index;
      GenSortIndirect<double>::sort (index, freqs);
      vector<MSReader*> oldReaders (itsReaders);
      for (int i=0; i<nband; ++i) {
        itsReaders[i] = oldReaders[index[i]];
      }      
    }

    void MultiMSReader::fillBands()
    {
      if (itsOrderMS)
        throw Exception("Cannot order the MSs if some are missing");
      // Get channel width (which should be the same for all bands).
      double chanw = itsReaders[itsFirst]->getInfo().chanWidths().data()[0];
      // Get frequency for first subband.
      double freq  = itsReaders[itsFirst]->getInfo().chanFreqs().data()[0];
      freq -= itsFirst*itsFillNChan*chanw;
      // Add missing channels to the total nr.
      itsNrChan += itsNMissing*itsFillNChan;
      // Collect the channel info of all MSs.
      Vector<double> chanFreqs (itsNrChan);
      Vector<double> chanWidths(itsNrChan);
      uint inx = 0;
      // Data for a missing MS can only be inserted if all other MSs have
      // the same nr of channels and are in increasing order of freq.
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          if (itsReaders[i]->getInfo().nchan() != itsFillNChan)
            throw Exception("An MS is missing; the others should have equal nchan");
          // Check if all channels have the same width and are consecutive.
          const Vector<double>& freqs = itsReaders[i]->getInfo().chanFreqs();
          const Vector<double>& width = itsReaders[i]->getInfo().chanWidths();
          if(freqs[0] < freq && !near(freqs[0], freq, 1e-5))
            throw Exception("Subbands should be in increasing order of frequency; found "
              + std::to_string(freqs[0]) + ", expected " + std::to_string(freq) + " (diff="
              + std::to_string(freqs[0]-freq) + ')');
          freq = freqs[itsFillNChan-1] + width[itsFillNChan-1];
          objcopy (chanFreqs.data()  + inx, freqs.data(), itsFillNChan);
          objcopy (chanWidths.data() + inx, width.data(), itsFillNChan);
          inx += itsFillNChan;
        } else {
          // Insert channel info for missing MSs.
          for (uint j=0; j<itsFillNChan; ++j) {
            chanFreqs[inx]  = freq;
            chanWidths[inx] = chanw;
            freq += chanw;
            inx++;
          }
        }
      }

      info().set (chanFreqs, chanWidths);
    }

    bool MultiMSReader::process (const DPBuffer& buf)
    {
      // Stop if at end.
      if (! itsReaders[itsFirst]->process (buf)) {
        return false;   // end of input
      }
      const DPBuffer& buf1 = itsReaders[itsFirst]->getBuffer();
      itsBuffer.setTime     (buf1.getTime());
      itsBuffer.setExposure (buf1.getExposure());
      itsBuffer.setRowNrs   (buf1.getRowNrs());
      // Size the buffers.
      if (itsBuffer.getFlags().empty()) {
        if (itsReadVisData) {
          itsBuffer.getData().resize (IPosition(3, itsNrCorr,
                                                itsNrChan, itsNrBl));
        }
        itsBuffer.getFlags().resize (IPosition(3, itsNrCorr,
                                               itsNrChan, itsNrBl));
      }
      // Loop through all readers and get data and flags.
      IPosition s(3, 0, 0, 0);
      IPosition e(3, itsNrCorr-1, 0, itsNrBl-1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          if (int(i) != itsFirst) {
            itsReaders[i]->process (buf);
          }
          const DPBuffer& msBuf = itsReaders[i]->getBuffer();
          if (msBuf.getRowNrs().empty())
            throw Exception(
                     "When using multiple MSs, the times in all MSs have to be "
                     "consecutive; this is not the case for MS " + std::to_string(i));
          // Copy data and flags.
          e[1] = s[1] + itsReaders[i]->getInfo().nchan() - 1;
          if (itsReadVisData) {
            itsBuffer.getData()(s,e) = msBuf.getData();
          }
          itsBuffer.getFlags()(s,e) = msBuf.getFlags();
        } else {
          e[1] = s[1] + itsFillNChan - 1;
          if (itsReadVisData) {
            itsBuffer.getData()(s,e) = Complex();
          }
          itsBuffer.getFlags()(s,e) = true;
        }          
        s[1] = e[1] + 1;
      }
      getNextStep()->process (itsBuffer);
      return true;
    }

    void MultiMSReader::finish()
    {
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          itsReaders[i]->finish();
        }
      }
      getNextStep()->finish();
    }

    void MultiMSReader::updateInfo (const DPInfo& infoIn)
    {
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          itsReaders[i]->updateInfo (infoIn);
        }
      }
      info() = itsReaders[itsFirst]->getInfo();
      // Use the first valid MS as the standard MS (for meta data)
      // Get meta data and check they are equal for all MSs.
      itsMS              = itsReaders[itsFirst]->table();
      itsStartTime       = getInfo().startTime();
      itsFirstTime       = itsReaders[itsFirst]->firstTime();
      itsLastTime        = itsReaders[itsFirst]->lastTime();
      itsTimeInterval    = getInfo().timeInterval();
      itsSelBL           = itsReaders[itsFirst]->baselineSelection();
      itsSpw             = itsReaders[itsFirst]->spectralWindow();
      itsNrCorr          = getInfo().ncorr();
      itsNrBl            = getInfo().nbaselines();
      itsNrChan          = 0;
      itsFillNChan       = getInfo().nchan();
      itsStartChan       = itsReaders[itsFirst]->startChan();
      itsFullResNChanAvg = itsReaders[itsFirst]->nchanAvgFullRes();
      itsFullResNTimeAvg = itsReaders[itsFirst]->ntimeAvgFullRes();
      itsHasFullResFlags = itsReaders[itsFirst]->hasFullResFlags();
      itsBaseRowNrs      = itsReaders[itsFirst]->getBaseRowNrs();
      for (uint i=0; i<itsMSNames.size(); ++i) {
        if (itsReaders[i]) {
          const DPInfo& rdinfo = itsReaders[i]->getInfo();
          if (!near(itsStartTime, rdinfo.startTime()))
            throw Exception("Start time of MS " + itsMSNames[i]
              + " differs from " + itsMSNames[itsFirst]);
          if (!near(itsLastTime, itsReaders[i]->lastTime()))
            throw Exception("Last time of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          if (!near(itsTimeInterval, rdinfo.timeInterval()))
            throw Exception("Time interval of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          if (itsNrCorr != rdinfo.ncorr())
            throw Exception("Number of correlations of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          if (itsNrBl != rdinfo.nbaselines())
            throw Exception("Number of baselines of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          if (itsFullResNChanAvg != itsReaders[i]->nchanAvgFullRes())
            throw Exception("FullResNChanAvg of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          if (itsFullResNTimeAvg != itsReaders[i]->ntimeAvgFullRes())
            throw Exception("FullResNTimeAvg of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          if (getInfo().antennaSet() != rdinfo.antennaSet())
            throw Exception("Antenna set of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          if (!allEQ (getInfo().getAnt1(), rdinfo.getAnt1()))
            throw Exception("Baseline order (ant1) of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          if (!allEQ (getInfo().getAnt2(), rdinfo.getAnt2()))
            throw Exception("Baseline order (ant2) of MS " + itsMSNames[i]
                     + " differs from " + itsMSNames[itsFirst]);
          itsNrChan += rdinfo.nchan();
          itsHasFullResFlags = (itsHasFullResFlags  &&
                                itsReaders[i]->hasFullResFlags());
        }
      }
      // Handle the bands and take care of missing MSs.
      // Sort them if needed.
      handleBands();

      // check that channels are regularly spaced, give warning otherwise
      if (itsNrChan>1) {
        Vector<Double> upFreq = info().chanFreqs()(
                                  Slicer(IPosition(1,1),
                                         IPosition(1,itsNrChan-1)));
        Vector<Double> lowFreq = info().chanFreqs()(
                                  Slicer(IPosition(1,0),
                                         IPosition(1,itsNrChan-1)));
        Double freqstep0=upFreq(0)-lowFreq(0);
        // Compare up to 1kHz accuracy
        itsRegularChannels=allNearAbs(upFreq-lowFreq, freqstep0, 1.e3) &&
                           allNearAbs(info().chanWidths(),
                                      info().chanWidths()(0), 1.e3);
      }

      // Set correct nr of channels.
      info().setNChan (itsNrChan);
      // Initialize the flag counters.
      itsFlagCounter.init (getInfo());
    }

    void MultiMSReader::show (std::ostream& os) const
    {
      os << "MultiMSReader" << std::endl;
      os << "  input MSs:      " << itsMSNames[0] << std::endl;
      for (uint i=1; i<itsMSNames.size(); ++i) {
        os << "                  " << itsMSNames[i] << std::endl;
      }
      if (! itsSelBL.empty()) {
        os << "  baseline:       " << itsSelBL << std::endl;
      }
      os << "  band            " << itsSpw << std::endl;
      os << "  startchan:      " << itsStartChan << "  (" << itsStartChanStr
         << ')' << std::endl;
      os << "  nchan:          " << itsNrChan << "  (" << itsNrChanStr
         << ')';
      if (itsRegularChannels) {
        os <<" (regularly spaced)" << std::endl;
      } else {
        os <<" (NOT regularly spaced)" << std::endl;
      }
      os << "  ncorrelations:  " << itsNrCorr << std::endl;
      os << "  nbaselines:     " << itsNrBl << std::endl;
      os << "  ntimes:         " << itsMS.nrow() / itsNrBl << std::endl;
      os << "  time interval:  " << itsTimeInterval << std::endl;
      os << "  DATA column:    " << itsDataColName << std::endl;
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          if (itsReaders[i]->missingData()) {
            os << "      column missing in  " << itsMSNames[i] << std::endl;
          }
        } else {
          os << "      MS missing         " << itsMSNames[i] << std::endl;
        }
      }
      os << "  WEIGHT column:  " << itsWeightColName << std::endl;
      os << "  autoweight:     " << boolalpha << itsAutoWeight << std::endl;
    }

    void MultiMSReader::showCounts (std::ostream& os) const
    {
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          itsReaders[i]->showCounts (os);
        }
      }
    }

    void MultiMSReader::showTimings (std::ostream& os, double duration) const
    {
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          itsReaders[i]->showTimings (os, duration);
        }
      }
    }

    void MultiMSReader::getUVW (const RefRows& rowNrs,
                                double time, DPBuffer& buf)
    {
      // All MSs have the same UVWs, so use first one.
      itsReaders[itsFirst]->getUVW (rowNrs, time, buf);
    }

    void MultiMSReader::getWeights (const RefRows& rowNrs, DPBuffer& buf)
    {
      Cube<float>& weights = buf.getWeights();
      // Resize if needed (probably when called for first time).
      if (weights.empty()) {
        weights.resize (itsNrCorr, itsNrChan, itsNrBl);
      }
      IPosition s(3, 0, 0, 0);
      IPosition e(3, itsNrCorr-1, 0, itsNrBl-1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          uint nchan = itsReaders[i]->getInfo().nchan();
          e[1] = s[1] + nchan-1;
          itsReaders[i]->getWeights (rowNrs, itsBuffers[i]);
          weights(s,e) = itsBuffers[i].getWeights();
        } else {
          e[1] = s[1] + itsFillNChan-1;
          weights(s,e) = float(0);
        }
        s[1] = e[1] + 1;
      }
    }

    bool MultiMSReader::getFullResFlags (const RefRows& rowNrs,
                                         DPBuffer& buf)
    {
      Cube<bool>& flags = buf.getFullResFlags();
      // Resize if needed (probably when called for first time).
      if (flags.empty()) {
        int norigchan = itsNrChan * itsFullResNChanAvg;
        flags.resize (norigchan, itsFullResNTimeAvg, itsNrBl);
      }
      // Return false if no fullRes flags available.
      if (!itsHasFullResFlags) {
        flags = false;
        return false;
      }
      // Flag everything if data rows are missing.
      if (rowNrs.rowVector().empty()) {
        flags = true;
        return true;
      }
      // Get the flags from all MSs and combine them.
      IPosition s(3, 0);
      IPosition e(flags.shape() - 1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          itsReaders[i]->getFullResFlags (rowNrs, itsBuffers[i]);
          e[0] = s[0] + itsBuffers[i].getFullResFlags().shape()[0] - 1;
          flags(s,e) = itsBuffers[i].getFullResFlags();
        } else {
          e[0] = s[0] + itsFillNChan - 1;
          flags(s,e) = true;
        }
        s[0] = e[0] + 1;
      }
      return true;
    }

  } //# end namespace
}
