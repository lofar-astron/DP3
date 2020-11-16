// VdsDesc.cc: Describe an entire visibility data set
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen <diepen AT astron nl>

#include "VdsDesc.h"

#include "StreamUtil.h"

#include <ostream>

namespace std {
using DP3::operator<<;
}
using namespace std;

namespace DP3 {
namespace CEP {

VdsDesc::VdsDesc(const VdsPartDesc& desc) : itsDesc(desc) {}

VdsDesc::VdsDesc(const string& parsetName) { init(ParameterSet(parsetName)); }

void VdsDesc::init(const ParameterSet& parset) {
  itsDesc = VdsPartDesc(parset);
  int npart = parset.getInt32("NParts");
  for (int i = 0; i < npart; ++i) {
    ostringstream prefix;
    prefix << "Part" << i << '.';
    ParameterSet subset = parset.makeSubset(prefix.str());
    itsParts.push_back(VdsPartDesc(subset));
  }
}

void VdsDesc::write(ostream& os) const {
  itsDesc.write(os, "");
  os << "NParts = " << itsParts.size() << endl;
  for (unsigned i = 0; i < itsParts.size(); ++i) {
    ostringstream prefix;
    prefix << "Part" << i << '.';
    itsParts[i].write(os, prefix.str());
  }
}

//   int VdsDesc::findPart (const string& fileSystem,
// 			     const vector<int>& done) const
//   {
//     for (unsigned i=0; i<itsParts.size(); ++i) {
//       if (done[i] < 0  &&  itsParts[i].getFileSys() == fileSystem) {
// 	return i;
//       }
//     }
//     return -1;
//   }

}  // namespace CEP
}  // namespace DP3
