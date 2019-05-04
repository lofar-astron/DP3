//# ParmSet.cc: Set of parameters to be used
//#
//# Copyright (C) 2008
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
//# $Id: ParmSet.cc 16977 2010-12-20 08:40:36Z diepen $

#include "ParmSet.h"
#include "ParmDB.h"

namespace DP3 {
namespace BBS {

  ParmSet::ParmSet()
  {}

  ParmId ParmSet::addParm (ParmDB& parmdb, const string& name)
  {
    // See if the parm has already been added.
    std::map<std::string,int>::const_iterator pos = itsNames.find(name);
    if (pos != itsNames.end()) {
      return itsParms[pos->second].getParmId();
    }
    // A new parameter name, so add it.
    // First find its name id.
    int nameId = parmdb.getNameId (name);
    // Assign a unique parm id.
    ParmId parmId = itsParms.size();
    itsParms.push_back (ParmKey(&parmdb, name, nameId, parmId));
    itsNames.insert (make_pair(name, parmId));
    // If needed, add its ParmDB to the list of used ParmDBs.
    unsigned int i;
    for (i=0; i<itsDBs.size(); ++i) {
      if (&parmdb == itsDBs[i]) {
        break;
      }
    }
    if (i == itsDBs.size()) {
      itsDBs.push_back (&parmdb);
    }
    return parmId;
  }

  ParmId ParmSet::find (const std::string& name) const
  {
    std::map<string,int>::const_iterator pos = itsNames.find(name);
    if (pos == itsNames.end())
			throw std::runtime_error("Parm " + name
               + " not found in ParmSet");
    return itsParms[pos->second].getParmId();
  }

  void ParmSet::getValues (std::vector<ParmValueSet>& vsets,
                           const Box& workDomain) const
  {
    if (vsets.size() == itsParms.size()) {
      return;      // nothing to do
    }
    assert (vsets.size() < itsParms.size());
    unsigned int todo = itsParms.size() - vsets.size();
    unsigned int start = vsets.size();
    vsets.resize (itsParms.size());
    std::vector<unsigned int>   nameIds;
    std::vector<ParmId> parmIds;
    nameIds.reserve (todo);
    parmIds.reserve (todo);
    // Get values for all new parameters.
    // Do it in order of ParmDB to query as little as possible.
    for (unsigned int i=0; i<itsDBs.size(); ++i) {
      ParmDB* pdb = itsDBs[i];
      for (unsigned int j=start; j<itsParms.size(); ++j) {
        if (itsParms[j].getParmDBPtr() == pdb) {
          int nameid = itsParms[j].getNameId();
          ParmId parmid = itsParms[j].getParmId();
          if (nameid < 0) {
            // Use default value.
            vsets[parmid] = pdb->getDefValue (itsParms[parmid].getName(),
                                              ParmValue());
            // If it is a polynomial with a domain, rescale it to the
            // work domain (if needed).
            if (vsets[parmid].getType() == ParmValue::Polc) {
              rescale (vsets[parmid], workDomain);
            }
          } else {
            nameIds.push_back (nameid);
            parmIds.push_back (parmid);
          }
        }
      }
      if (nameIds.size() > 0) {
        pdb->getValues (vsets, nameIds, parmIds, workDomain);
        todo -= nameIds.size();
        if (todo == 0) break;
        nameIds.clear();
        parmIds.clear();
      }
    }
  }

  void ParmSet::write (unsigned int parmId, ParmValueSet& pvset)
  {
    ParmDB* pdb = const_cast<ParmDB*>(itsParms[parmId].getParmDBPtr());
    pdb->putValues (itsParms[parmId].getName(),
                    itsParms[parmId].getNameId(),
                    pvset);
  }

  void ParmSet::clear()
  {
    itsDBs.clear();
    itsParms.clear();
    itsNames.clear();
  }

  void ParmSet::rescale (ParmValueSet& pset, const Box& newDomain) const
  {
    ParmValue defParm = pset.getDefParmValue();
    // If rescale is done, reset the ParmValueSet.
    if (defParm.rescale (newDomain.lowerX(), newDomain.upperX(),
                         newDomain.lowerY(), newDomain.upperY(),
                         pset.getScaleDomain())) {
      ParmValueSet psetNew(defParm, ParmValue::Polc,
                           pset.getPerturbation(), pset.getPertRel(),
                           newDomain);
      psetNew.setSolvableMask (pset.getSolvableMask());
      pset = psetNew;
    }
  }

} //# end namespace BBS
} //# end namspace LOFAR
