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
#include <Common/ParameterSet.h>
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
                                  const ParameterSet& parset,
                                  const string& prefix)
      : itsFirst    (-1),
        itsNMissing (0),
        itsMSNames  (msNames)
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
      ASSERTSTR (itsFirst>=0, "All input MeasurementSets do not exist");
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
      ASSERTSTR (!itsOrderMS, "Cannot order the MSs if some are missing");
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
          ASSERTSTR (itsReaders[i]->getInfo().nchan() == itsFillNChan,
                     "An MS is missing; the others should have equal nchan");
          // Check if all channels have the same width and are consecutive.
          const Vector<double>& freqs = itsReaders[i]->getInfo().chanFreqs();
          const Vector<double>& width = itsReaders[i]->getInfo().chanWidths();
          ASSERTSTR (freqs[0] > freq  ||  near(freqs[0], freq, 1e-5),
                     "Subbands should be in increasing order of frequency; found "
		     << freqs[0] << ", expected " << freq << " (diff="
		     << freqs[0]-freq << ')');
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
      itsBuffer.setUVW      (buf1.getUVW());
      // Size the buffers.
      if (itsReadVisData) {
        itsBuffer.getData().resize (IPosition(3, itsNrCorr,
                                              itsNrChan, itsNrBl));
      }
      itsBuffer.getFlags().resize (IPosition(3, itsNrCorr,
                                             itsNrChan, itsNrBl));
      // Loop through all readers and get data and flags.
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
          ASSERTSTR (near(itsStartTime, rdinfo.startTime())  &&
                     near(itsLastTime, itsReaders[i]->lastTime())  &&
                     near(itsTimeInterval, rdinfo.timeInterval())  &&
                     itsNrCorr == rdinfo.ncorr()  &&
                     itsNrBl   == rdinfo.nbaselines()  &&
                     itsFullResNChanAvg == itsReaders[i]->nchanAvgFullRes()  &&
                     itsFullResNTimeAvg == itsReaders[i]->ntimeAvgFullRes()  &&
                     allEQ (getInfo().getAnt1(), rdinfo.getAnt1())  &&
                     allEQ (getInfo().getAnt2(), rdinfo.getAnt2()),
                     "Meta data of MS " << itsMSNames[i]
                     << " differs from " << itsMSNames[itsFirst]);
          itsNrChan += rdinfo.nchan();
          itsHasFullResFlags = (itsHasFullResFlags  &&
                                itsReaders[i]->hasFullResFlags());
        }
      }
      // Handle the bands and take care of missing MSs.
      // Sort them if needed.
      handleBands();
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
         << ')' << std::endl;
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
          uint nchan = itsReaders[i]->getInfo().nchan();
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

    /*
    Cube<Complex> MultiMSReader::getData (const String& columnName,
                                          const RefRows& rowNrs)
    {
      Cube<Complex> data(itsNrCorr, itsNrChan, itsNrBl);
      IPosition s(3, 0, 0, 0);
      IPosition e(3, itsNrCorr-1, 0, itsNrBl-1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]) {
          uint nchan = itsReaders[i]->getInfo().nchan();
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
    */

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
