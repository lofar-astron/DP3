//# VdsPartDesc.cc: Description of a visibility data set or part thereof
//#
//# Copyright (C) 2007
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
//# $Id: VdsPartDesc.cc 16886 2010-12-08 10:43:17Z diepen $
//#
//# @author Ger van Diepen <diepen AT astron nl>

#include "VdsPartDesc.h"
#include "ParameterHandler.h"
#include "StreamUtil.h"

#include "../Blob/BlobArray.h"

#include <casacore/casa/Quanta/MVTime.h>

#include <ostream>

namespace std {
  using DP3::operator<<;
}
using namespace std;
using namespace casacore;

namespace DP3 { namespace CEP {

  VdsPartDesc::VdsPartDesc (const ParameterSet& parset)
  {
    itsName         = parset.getString ("Name");
    itsFileName     = parset.getString ("FileName", "");
    itsFileSys      = parset.getString ("FileSys", "");
    itsCDescName    = parset.getString ("ClusterDesc", "");
    itsStepTime     = parset.getDouble ("StepTime");
    itsNChan        = parset.getInt32Vector ("NChan", vector<int32_t>());
    itsStartFreqs   = parset.getDoubleVector ("StartFreqs", vector<double>());
    itsEndFreqs     = parset.getDoubleVector ("EndFreqs", vector<double>());
    itsParms        = parset.makeSubset ("Extra.");
    string timeStr;
    Quantity q;
    timeStr = parset.getString ("StartTime");
    if (!MVTime::read (q, timeStr, true))
      throw std::runtime_error("Could not parse StartTime");
    itsStartTime = q.getValue ("s");
    timeStr = parset.getString ("EndTime");
    if (!MVTime::read (q, timeStr, true))
      throw std::runtime_error("Could not parse EndTime");
    itsEndTime = q.getValue ("s");
    itsStartTimes = parset.getDoubleVector ("StartTimesDiff", vector<double>());
    itsEndTimes   = parset.getDoubleVector ("EndTimesDiff",   vector<double>());
    if (itsStartTimes.size() != itsEndTimes.size())
      throw std::runtime_error("StartTimeDiff and EndTimeDiff arrays had different sizes");
    double diff = itsStartTime;
    for (unsigned int i=0; i<itsStartTimes.size(); ++i) {
      itsStartTimes[i] += diff;
      diff += itsStepTime;
      itsEndTimes[i]   += diff;
    }
  }

  void VdsPartDesc::write (std::ostream& os, const std::string& prefix) const
  {
    os << prefix << "Name       = " << itsName << endl;
    if (! itsFileName.empty()) {
      os << prefix << "FileName   = " << itsFileName << endl;
    }
    if (! itsFileSys.empty()) {
      os << prefix << "FileSys    = " << itsFileSys << endl;
    }
    if (! itsCDescName.empty()) {
      os << prefix << "ClusterDesc= " << itsCDescName << endl;
    }
    os << prefix << "StartTime  = "
	<< MVTime::Format(MVTime::YMD,9) << MVTime(itsStartTime/86400) << endl;
    os << prefix << "EndTime    = "
	<< MVTime::Format(MVTime::YMD,9) << MVTime(itsEndTime/86400) << endl;
    os << prefix << "StepTime   = " << itsStepTime << endl;
    if (! itsStartTimes.empty()) {
      os << prefix << "StartTimesDiff=[";
      streamsize oldPrec = os.precision (5);
      double diff = itsStartTime;
      for (unsigned int i=0; i<itsStartTimes.size(); ++i) {
        if (i!=0) os << ',';
        os << itsStartTimes[i] - diff;
        diff += itsStepTime;
      }
      os << ']' << endl;
      os.precision (oldPrec);
    }
    if (! itsEndTimes.empty()) {
      os << prefix << "EndTimesDiff=[";
      streamsize oldPrec = os.precision (5);
      double diff = itsStartTime;
      for (unsigned int i=0; i<itsEndTimes.size(); ++i) {
        if (i!=0) os << ',';
        diff += itsStepTime;
        os << itsEndTimes[i] - diff;
      }
      os << ']' << endl;
      os.precision (oldPrec);
    }
    if (! itsNChan.empty()) {
      os << prefix << "NChan      = " << itsNChan << endl;
      streamsize oldPrec = os.precision (12);
      os << prefix << "StartFreqs = " << itsStartFreqs << endl;
      os << prefix << "EndFreqs   = " << itsEndFreqs << endl;
      os.precision (oldPrec);
    }
    // Prepend the extra parameters with Extra..
    ParameterSet parms;
    parms.adoptCollection (itsParms, prefix+"Extra.");
    parms.writeStream (os);
  }

  void VdsPartDesc::setName (const std::string& name,
                             const std::string& fileSys)
  {
    itsName    = name;
    itsFileSys = fileSys;
  }

  void VdsPartDesc::changeBaseName (const string& newBaseName)
  {
    string::size_type pos = itsName.rfind ('/');
    if (pos == string::npos) {
      itsName = newBaseName;
    } else {
      itsName = itsName.substr (0, pos+1) + newBaseName;
    }
  }

  void VdsPartDesc::setTimes (double startTime, double endTime, double stepTime,
                              const vector<double>& startTimes,
                              const vector<double>& endTimes)
  {
    itsStartTime  = startTime;
    itsEndTime    = endTime;
    itsStepTime   = stepTime;
    itsStartTimes = startTimes;
    itsEndTimes   = endTimes;
  }

  void VdsPartDesc::addBand (int nchan, double startFreq, double endFreq)
  {
    itsNChan.push_back (nchan);
    double step = (endFreq-startFreq)/nchan;
    for (int i=0; i<nchan; ++i) {
      itsStartFreqs.push_back (startFreq);
      startFreq += step;
      itsEndFreqs.push_back (startFreq);
    }
  }

  void VdsPartDesc::addBand (int nchan, const vector<double>& startFreq,
    const vector<double>& endFreq)
  {
    assert (startFreq.size() == endFreq.size());
    assert (int(startFreq.size())==nchan || startFreq.size() == 1);
    itsNChan.push_back (nchan);
    for (unsigned int i=0; i<startFreq.size(); ++i) {
      itsStartFreqs.push_back (startFreq[i]);
      itsEndFreqs.push_back   (endFreq[i]);
    }
  }

  BlobOStream& VdsPartDesc::toBlob (BlobOStream& bs) const
  {
    bs.putStart ("VdsPartDesc", 1);
    bs << itsName << itsFileName << itsFileSys << itsCDescName
       << itsStartTime << itsEndTime << itsStepTime
       << itsStartTimes << itsEndTimes
       << itsNChan << itsStartFreqs << itsEndFreqs
       << itsParms;
    bs.putEnd();
    return bs;
  }

  BlobIStream& VdsPartDesc::fromBlob (BlobIStream& bs)
  {
    bs.getStart ("VdsPartDesc");
    bs >> itsName >> itsFileName >> itsFileSys >> itsCDescName
       >> itsStartTime >> itsEndTime >> itsStepTime
       >> itsStartTimes >> itsEndTimes
       >> itsNChan >> itsStartFreqs >> itsEndFreqs
       >> itsParms;
    bs.getEnd();
    return bs;
  }

}} // end namespaces
