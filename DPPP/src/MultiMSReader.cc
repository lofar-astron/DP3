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
      for (uInt i=0; i<msNames.size(); ++i) {
        itsReaders.push_back (new MSReader (msNames[i], parset, prefix));
        itsSteps.push_back   (DPStep::ShPtr(itsReaders[i]));
        // Add a null step for the reader.
        itsSteps[i]->setNextStep (nullStep);
      }
      // Use the first MS as the standard MS (for meta data)
      // Get meta data and check they are equal for all MSs.
      itsMS              = itsReaders[0]->table();
      itsStartTime       = itsReaders[0]->startTime();
      itsFirstTime       = itsReaders[0]->firstTime();
      itsLastTime        = itsReaders[0]->lastTime();
      itsInterval        = itsReaders[0]->timeInterval();
      itsSelBL           = itsReaders[0]->baselineSelection();
      itsSpw             = itsReaders[0]->spectralWindow();
      itsNrCorr          = itsReaders[0]->ncorr();
      itsNrBl            = itsReaders[0]->nbaselines();
      itsNrChan          = itsReaders[0]->nchan();
      itsStartChan       = itsReaders[0]->startChan();
      itsFullResNChanAvg = itsReaders[0]->nchanAvg();
      itsFullResNTimeAvg = itsReaders[0]->ntimeAvg();
      itsAnt1            = itsReaders[0]->getAnt1();
      itsAnt2            = itsReaders[0]->getAnt2();
      itsAntNames        = itsReaders[0]->antennaNames();
      itsAntPos          = itsReaders[0]->antennaPos();
      itsArrayPos        = itsReaders[0]->arrayPos();
      itsPhaseCenter     = itsReaders[0]->phaseCenter();
      itsDelayCenter     = itsReaders[0]->delayCenter();
      itsTileBeamDir     = itsReaders[0]->tileBeamDir();
      itsHasFullResFlags = itsReaders[0]->hasFullResFlags();
      itsBaseRowNrs      = itsReaders[0]->getBaseRowNrs();
      for (uInt i=1; i<msNames.size(); ++i) {
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
                   << " differs from first one");
        itsNrChan += itsReaders[i]->nchan();
        itsHasFullResFlags = (itsHasFullResFlags  &&
                              itsReaders[i]->hasFullResFlags());
      }
      // If needed, sort the MSs in order of frequency.
      if (itsOrderMS) {
        sortBands();
      }
      // Collect the channel info of all MSs.
      itsChanFreqs.resize  (itsNrChan);
      itsChanWidths.resize (itsNrChan);
      uint inx = 0;
      for (uInt i=0; i<itsReaders.size(); ++i) {
        uint nchan = itsReaders[i]->nchan();
        objcopy (itsChanFreqs.data()  + inx,
                 itsReaders[i]->chanFreqs().data(),  nchan);
        objcopy (itsChanWidths.data() + inx,
                 itsReaders[i]->chanWidths().data(), nchan);
        inx += nchan;
      }      
      // Initialize the flag counters.
      itsFlagCounter.init (itsNrBl, itsNrChan, itsNrCorr);
    }

    MultiMSReader::~MultiMSReader()
    {}

    void MultiMSReader::setReadVisData (bool readVisData)
    {
      itsReadVisData = readVisData;
      for (uint i=0; i<itsReaders.size(); ++i) {
        itsReaders[i]->setReadVisData (readVisData);
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

    void MultiMSReader::getFreqInfo (Vector<double>& freq,
                                     Vector<double>& width,
                                     Vector<double>& effBW,
                                     Vector<double>& resolution) const
    {
      freq.resize (itsNrChan);
      width.resize (itsNrChan);
      effBW.resize (itsNrChan);
      resolution.resize (itsNrChan);
      IPosition s(1,0);
      IPosition e(1,0);
      for (uint i=0; i<itsReaders.size(); ++i) {
        e[0] = s[0] + itsReaders[i]->nchan() - 1;
        Vector<double> subfreq (freq(s,e));
        Vector<double> subwidth (width(s,e));
        Vector<double> subeffBW (effBW(s,e));
        Vector<double> subresolution (resolution(s,e));
        itsReaders[i]->getFreqInfo (subfreq, subwidth,
                                    subeffBW, subresolution);
        s[0] += itsReaders[i]->nchan();
      }
    }

    bool MultiMSReader::process (const DPBuffer& buf)
    {
      IPosition s(3, 0, 0, 0);
      IPosition e(3, itsNrCorr-1, 0, itsNrBl-1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (! itsReaders[i]->process (buf)) {
          return false;   // end of input
        }
        const DPBuffer& msBuf = itsReaders[i]->getBuffer();
        ASSERTSTR (! msBuf.getRowNrs().empty(),
                   "When using multiple MSs, the times in all MSs have to be "
                   "consecutive; this is not the case for MS " << i);
        if (i == 0) {
          itsBuffer.setTime (msBuf.getTime());
          itsBuffer.setRowNrs (msBuf.getRowNrs());
          itsBuffer.setUVW (msBuf.getUVW());
          if (! msBuf.getData().empty()) {
            itsBuffer.getData().resize (IPosition(3, itsNrCorr,
                                                  itsNrChan, itsNrBl));
          }
          if (! msBuf.getFlags().empty()) {
            itsBuffer.getFlags().resize (IPosition(3, itsNrCorr,
                                                   itsNrChan, itsNrBl));
          }
          if (! msBuf.getWeights().empty()) {
            itsBuffer.getWeights().resize (IPosition(3, itsNrCorr,
                                                     itsNrChan, itsNrBl));
          }
        }
        e[1] = s[1] + itsReaders[i]->nchan() - 1;
        if (! msBuf.getData().empty()) {
          itsBuffer.getData()(s,e) = msBuf.getData();
        }
        if (! msBuf.getFlags().empty()) {
          itsBuffer.getFlags()(s,e) = msBuf.getFlags();
        }
        if (! msBuf.getWeights().empty()) {
          itsBuffer.getWeights()(s,e) = msBuf.getWeights();
        }
        s[1] = e[1] + 1;
      }
      if (itsHasFullResFlags) {
        vector<Cube<bool> > fullResFlags;
        fullResFlags.reserve (itsReaders.size());
        for (uint i=0; i<itsReaders.size(); ++i) {
          fullResFlags.push_back (itsReaders[i]->getBuffer().getFullResFlags());
        }
        combineFullResFlags (fullResFlags, itsBuffer.getFullResFlags());
      }
      getNextStep()->process (itsBuffer);
      return true;
    }

    void MultiMSReader::finish()
    {
      for (uint i=0; i<itsReaders.size(); ++i) {
        itsReaders[i]->finish();
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
      os << "  input MSs:      " << itsReaders[0]->msName() << std::endl;
      for (uint i=1; i<itsReaders.size(); ++i) {
        os << "                  " << itsReaders[i]->msName() << std::endl;
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
      bool missing = false;
      for (uint i=0; i<itsReaders.size(); ++i) {
        if (itsReaders[i]->missingData()) {
          if (!missing) {
            os << "      missing in:   ";
            missing = true;
          } else {
            os << "                    ";
          }
          os << itsReaders[i]->msName() << std::endl;
        }
      }
      os << "  autoweight:     " << itsAutoWeight << std::endl;
    }

    void MultiMSReader::showCounts (std::ostream& os) const
    {
      for (uint i=0; i<itsReaders.size(); ++i) {
        itsReaders[i]->showCounts (os);
      }
    }

    void MultiMSReader::showTimings (std::ostream& os, double duration) const
    {
      for (uint i=0; i<itsReaders.size(); ++i) {
        itsReaders[i]->showTimings (os, duration);
      }
    }

    Matrix<double> MultiMSReader::getUVW (const RefRows& rowNrs)
    {
      // All MSs have the same UVWs, so use first one.
      return itsReaders[0]->getUVW (rowNrs);
    }

    Cube<float> MultiMSReader::getWeights (const RefRows& rowNrs,
                                           const DPBuffer& buf)
    {
      Cube<float> weights(itsNrCorr, itsNrChan, itsNrBl);
      IPosition s(3, 0, 0, 0);
      IPosition e(3, itsNrCorr-1, 0, itsNrBl-1);
      for (uint i=0; i<itsReaders.size(); ++i) {
        uint nchan = itsReaders[i]->nchan();
        e[1] = s[1] + nchan-1;
        weights(s,e) = itsReaders[i]->getWeights (rowNrs, buf);
        s[1] += nchan;
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
        fullResFlags.push_back (itsReaders[i]->getFullResFlags (rowNrs));
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
        uint nchan = itsReaders[i]->nchan();
        e[1] = s[1] + nchan-1;
        data(s,e) = itsReaders[i]->getData (columnName, rowNrs);
        s[1] += nchan;
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
        cout<<s<<e<<flags.shape()<<vec[i].shape()<<endl;
        flags(s,e) = vec[i];
        s[0] = e[0] + 1;
      }
    }

  } //# end namespace
}
