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

#include <lofar_config.h>
#include <DPPP/MultiMSReader.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <Common/StreamUtil.h>
#include <Common/LofarLogger.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ExprNode.h>
#include <tables/Tables/RecordGram.h>
#include <measures/Measures/MeasTable.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>
#include <measures/TableMeasures/ArrayMeasColumn.h>
#include <casa/Containers/Record.h>
#include <casa/Quanta/MVTime.h>
#include <casa/Utilities/GenSort.h>
#include <casa/OS/Conversion.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MultiMSReader::MultiMSReader (const vector<string>& msNames,
                                  const ParSet& parset, const string& prefix)
      : itsFirst   (-1),
        itsMSNames (msNames)
    {
      ASSERTSTR (msNames.size() > 0, "No names of MeasurementSets given");
      itsStartChanStr     = parset.getString (prefix+"startchan", "0");
      itsNrChanStr        = parset.getString (prefix+"nchan", "0");
      itsUseFlags         = parset.getBool   (prefix+"useflag", true);
      itsDataColName      = parset.getString (prefix+"datacolumn", "DATA");
      itsMissingData      = parset.getBool   (prefix+"missingdata", false);
      itsAutoWeight       = parset.getBool   (prefix+"autoweight", false);
      itsNeedSort         = parset.getBool   (prefix+"sort", false);
      itsOrderMS          = parset.getBool   (prefix+"orderms", true);
      // Open all MSs.
      DPStep::ShPtr nullStep (new NullStep());
      itsReaders.reserve (msNames.size());
      itsSteps.reserve   (msNames.size());
      int nmissing = 0;
      for (uint i=0; i<msNames.size(); ++i) {
        itsReaders.push_back (new MSReader (msNames[i], parset, prefix,
                                            itsMissingData));
        itsSteps.push_back   (DPStep::ShPtr(itsReaders[i]));
        // Add a null step for the reader.
        itsSteps[i]->setNextStep (nullStep);
        // Ignore if the MS is missing.
        if (itsReaders[i]->table().isNull()) {
          itsReaders[i] = 0;
          nmissing++;
        } else if (itsFirst < 0) {
          itsFirst = i;
        }
      }
      ASSERTSTR (itsFirst>=0, "All input MeasurementSets do not exist");
      // Use the first valid MS as the standard MS (for meta data)
      // Get meta data and check they are equal for all MSs.
      itsMS              = itsReaders[itsFirst]->table();
      itsStartTime       = itsReaders[itsFirst]->startTime();
      itsFirstTime       = itsReaders[itsFirst]->firstTime();
      itsLastTime        = itsReaders[itsFirst]->lastTime();
      itsInterval        = itsReaders[itsFirst]->timeInterval();
      itsSelBL           = itsReaders[itsFirst]->baselineSelection();
      itsSpw             = itsReaders[itsFirst]->spectralWindow();
      itsNrCorr          = itsReaders[itsFirst]->ncorr();
      itsNrBl            = itsReaders[itsFirst]->nbaselines();
      itsNrChan          = 0;
      itsFillNChan       = itsReaders[itsFirst]->nchan();
      itsStartChan       = itsReaders[itsFirst]->startChan();
      itsFullResNChanAvg = itsReaders[itsFirst]->nchanAvg();
      itsFullResNTimeAvg = itsReaders[itsFirst]->ntimeAvg();
      itsAnt1            = itsReaders[itsFirst]->getAnt1();
      itsAnt2            = itsReaders[itsFirst]->getAnt2();
      itsAntNames        = itsReaders[itsFirst]->antennaNames();
      itsAntPos          = itsReaders[itsFirst]->antennaPos();
      itsArrayPos        = itsReaders[itsFirst]->arrayPos();
      itsPhaseCenter     = itsReaders[itsFirst]->phaseCenter();
      itsDelayCenter     = itsReaders[itsFirst]->delayCenter();
      itsTileBeamDir     = itsReaders[itsFirst]->tileBeamDir();
      itsHasFullResFlags = itsReaders[itsFirst]->hasFullResFlags();
      itsBaseRowNrs      = itsReaders[itsFirst]->getBaseRowNrs();
      for (uint i=0; i<msNames.size(); ++i) {
        if (itsReaders[i]) {
          ASSERTSTR (near(itsStartTime, itsReaders[i]->startTime())  &&
                     near(itsLastTime, itsReaders[i]->lastTime())  &&
                     near(itsInterval, itsReaders[i]->timeInterval())  &&
                     itsNrCorr == itsReaders[i]->ncorr()  &&
                     itsNrBl   == itsReaders[i]->nbaselines()  &&
                     itsFullResNChanAvg == itsReaders[i]->nchanAvg()  &&
                     itsFullResNTimeAvg == itsReaders[i]->ntimeAvg()  &&
                     allEQ (itsAnt1, itsReaders[i]->getAnt1())  &&
                     allEQ (itsAnt2, itsReaders[i]->getAnt2()),
                     "Meta data of MS " << msNames[i]
                     << " differs from " << msNames[itsFirst]);
          itsNrChan += itsReaders[i]->nchan();
          itsHasFullResFlags = (itsHasFullResFlags  &&
                                itsReaders[i]->hasFullResFlags());
        }
      }
      // Handle the bands and take care of missing MSs.
      // Sort them if needed.
      handleBands (nmissing);
      // Initialize the flag counters.
      itsFlagCounter.init (itsNrBl, itsNrChan, itsNrCorr);
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

    void MultiMSReader::handleBands (uint nmissing)
    {
      if (nmissing > 0) {
        fillBands (nmissing);
        return;
      }
      if (itsOrderMS) {
        sortBands();
      }
      // Collect the channel info of all MSs.
      itsChanFreqs.resize  (itsNrChan);
      itsChanWidths.resize (itsNrChan);
      uint inx = 0;
      for (uint i=0; i<itsReaders.size(); ++i) {
        uint nchan = itsReaders[i]->nchan();
        objcopy (itsChanFreqs.data()  + inx,
                 itsReaders[i]->chanFreqs().data(),  nchan);
        objcopy (itsChanWidths.data() + inx,
                 itsReaders[i]->chanWidths().data(), nchan);
        inx += nchan;
      }      
    }

    void MultiMSReader::sortBands()
    {
      // Order the bands (MSs) in order of frequency.
      int nband = itsReaders.size();
      Vector<double> freqs(nband);
      for (int i=0; i<nband; ++i) {
        freqs[i] = itsReaders[i]->chanFreqs().data()[0];
      }
      Vector<uInt> index;
      GenSortIndirect<double>::sort (index, freqs);
      vector<MSReader*> oldReaders (itsReaders);
      for (int i=0; i<nband; ++i) {
        itsReaders[i] = oldReaders[index[i]];
      }      
    }

    void MultiMSReader::fillBands (uint nmissing)
    {
      ASSERTSTR (!itsOrderMS, "Cannot order the MSs if some are missing");
      // Get channel width (which should be the same for all bands).
      double chanw = itsReaders[itsFirst]->chanWidths().data()[0];
      // Get frequency for first subband.
      double freq  = itsReaders[itsFirst]->chanFreqs().data()[0];
      freq -= itsFirst*itsFillNChan*chanw;
      // Add missing channels to the total nr.
      itsNrChan += nmissing*itsFillNChan;
      // Collect the channel info of all MSs.
      itsChanFreqs.resize  (itsNrChan);
      itsChanWidths.resize (itsNrChan);
      uint inx = 0;
      // Data for a missing MS can only be inserted if all other MSs have
      // the same nr of channels and are in increasing order of freq.
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          ASSERTSTR (itsReaders[i]->nchan() == itsFillNChan,
                     "An MS is missing; the others should have equal size");
          // Check if all channels have the same width and are consecutive.
          const Vector<double>& freqs = itsReaders[i]->chanFreqs();
          const Vector<double>& width = itsReaders[i]->chanWidths();
          ASSERTSTR (freqs[0] > freq  ||  near(freqs[0], freq),
                     "Subbands should be in increasing order of frequency");
          freq = freqs[itsFillNChan-1] + width[itsFillNChan-1];
          objcopy (itsChanFreqs.data()  + inx, freqs.data(), itsFillNChan);
          objcopy (itsChanWidths.data() + inx, width.data(), itsFillNChan);
          inx += itsFillNChan;
        } else {
          // Insert channel info for missing MSs.
          for (uint j=0; j<itsFillNChan; ++j) {
            itsChanFreqs[inx]  = freq;
            itsChanWidths[inx] = chanw;
            freq += chanw;
            inx++;
          }
        }
      }
    }

    void MultiMSReader::getFreqInfo (Vector<double>& freq,
                                     Vector<double>& width,
                                     Vector<double>& effBW,
                                     Vector<double>& resolution,
                                     double& refFreq) const
    {
      freq.resize (itsNrChan);
      width.resize (itsNrChan);
      effBW.resize (itsNrChan);
      resolution.resize (itsNrChan);
      IPosition s(1,0);
      IPosition e(1,0);
      double rf;
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          e[0] = s[0] + itsReaders[i]->nchan() - 1;
          Vector<double> subfreq (freq(s,e));
          Vector<double> subwidth (width(s,e));
          Vector<double> subeffBW (effBW(s,e));
          Vector<double> subresolution (resolution(s,e));
          itsReaders[i]->getFreqInfo (subfreq, subwidth,
                                      subeffBW, subresolution, rf);
        } else {
          e[0] = s[0] + itsFillNChan - 1;
          freq(s,e)       = itsChanFreqs(s,e);
          width(s,e)      = itsChanWidths(s,e);
          effBW(s,e)      = itsChanWidths(s,e);
          resolution(s,e) = itsChanWidths(s,e);
        }
        s[0] = e[0] + 1;
      }
      // Take the middle channel for the reference frequency.
      refFreq = freq[freq.size()/2];
    }

    bool MultiMSReader::process (const DPBuffer& buf)
    {
      // Stop if at end.
      if (! itsReaders[itsFirst]->process (buf)) {
        return false;   // end of input
      }
      const DPBuffer& buf1 = itsReaders[itsFirst]->getBuffer();
      itsBuffer.setTime   (buf1.getTime());
      itsBuffer.setRowNrs (buf1.getRowNrs());
      itsBuffer.setUVW    (buf1.getUVW());
      // Size the buffers.
      if (itsReadVisData) {
        itsBuffer.getData().resize (IPosition(3, itsNrCorr,
                                              itsNrChan, itsNrBl));
      }
      itsBuffer.getFlags().resize (IPosition(3, itsNrCorr,
                                             itsNrChan, itsNrBl));
      // Lopp through all readers and get data and flags.
      IPosition s(3, 0, 0, 0);
      IPosition e(3, itsNrCorr-1, 0, itsNrBl-1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          if (int(i) != itsFirst) {
            itsReaders[i]->process (buf);
          }
          const DPBuffer& msBuf = itsReaders[i]->getBuffer();
          ASSERTSTR (! msBuf.getRowNrs().empty(),
                     "When using multiple MSs, the times in all MSs have to be "
                     "consecutive; this is not the case for MS " << i);
          // Copy data and flags.
          e[1] = s[1] + itsReaders[i]->nchan() - 1;
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

    void MultiMSReader::updateInfo (DPInfo& info)
    {
      info.init (itsNrCorr, itsStartChan, itsNrChan, itsNrBl,
                 int((itsLastTime - itsFirstTime)/itsInterval + 1.5),
                 itsInterval);
      info.setPhaseCenter (itsPhaseCenter, true);
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
         << ')' << std::endl;
      os << "  ncorrelations:  " << itsNrCorr << std::endl;
      os << "  nbaselines:     " << itsNrBl << std::endl;
      os << "  ntimes:         " << itsMS.nrow() / itsNrBl << std::endl;
      os << "  time interval:  " << itsInterval << std::endl;
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
      os << "  autoweight:     " << itsAutoWeight << std::endl;
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

    Matrix<double> MultiMSReader::getUVW (const RefRows& rowNrs)
    {
      // All MSs have the same UVWs, so use first one.
      return itsReaders[itsFirst]->getUVW (rowNrs);
    }

    Cube<float> MultiMSReader::getWeights (const RefRows& rowNrs,
                                           const DPBuffer& buf)
    {
      Cube<float> weights(itsNrCorr, itsNrChan, itsNrBl);
      IPosition s(3, 0, 0, 0);
      IPosition e(3, itsNrCorr-1, 0, itsNrBl-1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          uint nchan = itsReaders[i]->nchan();
          e[1] = s[1] + nchan-1;
          weights(s,e) = itsReaders[i]->getWeights (rowNrs, buf);
        } else {
          e[1] = s[1] + itsFillNChan-1;
          weights(s,e) = float(0);
        }
        s[1] = e[1] + 1;
      }
      return weights;
    }

    Cube<bool> MultiMSReader::getFullResFlags (const RefRows& rowNrs)
    {
      // Return empty array if no fullRes flags.
      if (!itsHasFullResFlags  ||  rowNrs.rowVector().empty()) {
        return Cube<bool>();
      }
      Cube<bool> flags;
      vector<Cube<bool> > fullResFlags;
      fullResFlags.reserve (itsReaders.size());
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          fullResFlags.push_back (itsReaders[i]->getFullResFlags (rowNrs));
        } else {
          // Fill a cube for missing fullres flags only once.
          if (itsFullResCube.empty()) {
            itsFullResCube.resize (itsFillNChan*itsFullResNChanAvg,
                                   itsFullResNTimeAvg,
                                   itsNrBl);
            itsFullResCube = True;
          }
          fullResFlags.push_back (itsFullResCube);
        }
      }
      combineFullResFlags (fullResFlags, flags);
      return flags;
    }

    Cube<Complex> MultiMSReader::getData (const String& columnName,
                                          const RefRows& rowNrs)
    {
      Cube<Complex> data(itsNrCorr, itsNrChan, itsNrBl);
      IPosition s(3, 0, 0, 0);
      IPosition e(3, itsNrCorr-1, 0, itsNrBl-1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          uint nchan = itsReaders[i]->nchan();
          e[1] = s[1] + nchan-1;
          data(s,e) = itsReaders[i]->getData (columnName, rowNrs);
        } else {
          e[1] = s[1] + itsFillNChan-1;
          data(s,e) = Complex();
        }
        s[1] = e[1] + 1;
      }
      return data;
    }

    void MultiMSReader::combineFullResFlags (const vector<Cube<bool> >& vec,
                                             Cube<bool>& flags) const
    {
      // The cubes have axes nchan, ntimeavg, nbl.
      IPosition s(3, 0);
      IPosition e(vec[0].shape());
      // Count nr of channels.
      uint nchan = 0;
      for (uint i=0; i<vec.size(); ++i) {
        nchan += vec[i].shape()[0];
      }
      e[0] = nchan;
      flags.resize (e);
      e -= 1;
      for (uint i=0; i<vec.size(); ++i) {
        e[0] = s[0] + vec[i].shape()[0] - 1;
        flags(s,e) = vec[i];
        s[0] = e[0] + 1;
      }
    }

  } //# end namespace
}
