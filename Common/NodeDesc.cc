//# NodeDesc.cc: Description of a node
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
//# $Id: NodeDesc.cc 16886 2010-12-08 10:43:17Z diepen $
//#
//# @author Ger van Diepen <diepen AT astron nl>

#include "NodeDesc.h"

#include "StreamUtil.h"
#include "StringUtil.h"

#include <ostream>

#include <boost/algorithm/string/case_conv.hpp>

namespace std {
  using DP3::operator<<;
}
using namespace std;

namespace DP3 { namespace CEP {

  NodeDesc::NodeDesc (const ParameterSet& parset)
  {
    itsName = parset.getString ("NodeName");
    string type (boost::to_lower_copy(parset.getString ("NodeType", "Compute")));
    if (type == "compute") {
      itsType = Compute;
    } else if (type == "storage") {
      itsType = Storage;
    } else if (type == "head") {
      itsType = Head;
    } else {
      itsType = Any;
    }
    itsMounts  = parset.getStringVector ("NodeMountPoints", true);
    itsFileSys = parset.getStringVector ("NodeFileSys", itsMounts, true);
    assert (itsFileSys.size() == itsMounts.size());
    for (uint i=0; i<itsMounts.size(); ++i) {
      assert (itsFileSys[i].size() > 0);
      assert (itsMounts[i].size() > 0  &&  itsMounts[i][0] == '/');
    }
  }

  void NodeDesc::addFileSys (const string& fsName,
			     const string& mountPoint)
  {
    assert (fsName.size() > 0);
    string mp(mountPoint);
    if (mp.size() > 5  &&  mp.substr(0,5) == "/auto") {
      mp = mp.substr(5);
    }
    assert (mp.size() > 0  &&  mp[0] == '/');
    itsFileSys.push_back (fsName);
    itsMounts.push_back (mp);
  }

  string NodeDesc::findFileSys (const string& fileName) const
  {
    // The file name must be absolute.
    assert (fileName.size() > 1  &&  fileName[0] == '/');
    // Determine the max nr of parts in the mount point.
    // Remember the root filesys (a single /).
    int nrp = 0;
    int rootfs = -1;
    for (uint i=0; i<itsMounts.size(); ++i) {
      int nr=0;
      const string& str = itsMounts[i];
      // A single / counts as no part.
      if (str.size() == 1) {
	rootfs = i;
      } else {
	for (uint j=0; j<str.size(); ++j) {
	  if (str[j] == '/') {
	    ++nr;
	  }
	}
      }
      if (nr > nrp) {
	nrp = nr;
      }
    }
    // Find the slashes in the file name for each part.
    vector<int> pos(nrp, -1);
    int nr = 0;
    for (uint i=1; i<fileName.size() && nr<nrp; ++i) {
      if (fileName[i] == '/') {
	pos[nr++] = i;
      }
    }
    // Now compare if it matches the file name.
    // Start with the longest possible string.
    for (int p=nr-1; p>=0; --p) {
      string filePart = fileName.substr(0,pos[p]);
      for (uint i=0; i<itsMounts.size(); ++i) {
	if (filePart == itsMounts[i]) {
	  return itsFileSys[i];
	}
      }
    }
    // No match, so return root file system if there.
    // Otherwise return empty string.
    if (rootfs >= 0) {
      return itsFileSys[rootfs];
    }
    return "";
  }

  void NodeDesc::write (ostream& os, const string& prefix) const
  {
    string type = "Any";
    if (itsType == Compute) {
      type = "Compute";
    } else if (itsType == Storage) {
      type = "Storage";
    } else if (itsType == Head) {
      type = "Head";
    }
    os << prefix << "NodeName = " << itsName << endl;
    os << prefix << "NodeType = " << type << endl;
    os << prefix << "NodeFileSys     = " << itsFileSys << endl;
    os << prefix << "NodeMountPoints = " << itsMounts << endl;
  }

}} // end namespaces
