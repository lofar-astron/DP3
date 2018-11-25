//# SourceDBUtil.cc: Helper functions to extract patch and source information
//# from a SourceDB.
//#
//# Copyright (C) 2012
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

#include "SourceDBUtil.h"

#include "Exceptions.h"
#include "PointSource.h"
#include "GaussianSource.h"

#include "../ParmDB/SourceDB.h"

#include <sstream>
#include <set>

namespace DP3
{
namespace DPPP
{
using BBS::SourceDB;
using BBS::SourceData;
using BBS::SourceInfo;


std::vector<Patch::ConstPtr> makePatches(SourceDB &sourceDB,
                                    const std::vector<string> &patchNames,
                                    uint nModel)
{
  // Create a component list for each patch name.
  std::vector<std::vector<ModelComponent::Ptr> > componentsList(nModel);

  // Loop over all sources.
  sourceDB.lock();
  sourceDB.rewind();
  SourceData src;
  while (! sourceDB.atEnd()) {
    sourceDB.getNextSource (src);
    // Use the source if its patch matches a patch name.
    for (uint i=0; i<nModel; ++i) {
      if (src.getPatchName() == patchNames[i]) {
        // Fetch position.
        assert (src.getInfo().getRefType() == "J2000");
        Position position;
        position[0] = src.getRa();
        position[1] = src.getDec();

        // Fetch stokes vector.
        Stokes stokes;
        stokes.I = src.getI();
        stokes.V = src.getV();
        if(!src.getInfo().getUseRotationMeasure())
        {
          stokes.Q = src.getQ();
          stokes.U = src.getU();
        }

        PointSource::Ptr source;
        switch(src.getInfo().getType())
        {
        case SourceInfo::POINT:
            {
                source = PointSource::Ptr(new PointSource(position, stokes));
            }
            break;

        case SourceInfo::GAUSSIAN:
            {
                GaussianSource::Ptr gauss(new GaussianSource(position, stokes));

                const double deg2rad = (casacore::C::pi / 180.0);
                gauss->setPositionAngle(src.getOrientation() * deg2rad);

                const double arcsec2rad = (casacore::C::pi / 3600.0) / 180.0;
                gauss->setMajorAxis(src.getMajorAxis() * arcsec2rad);
                gauss->setMinorAxis(src.getMinorAxis() * arcsec2rad);
                source = gauss;
            }
            break;

        default:
            {
                throw Exception("Only point sources and Gaussian sources are"
                    " supported at this time.");
            }
        }

        // Fetch spectral index attributes (if applicable).
        bool isLogarithmic = src.getInfo().getHasLogarithmicSI();
        if (src.getSpectralTerms().size() > 0) {
          source->setSpectralTerms(src.getInfo().getSpectralTermsRefFreq(),
                                   isLogarithmic,
                                   src.getSpectralTerms().begin(),
                                   src.getSpectralTerms().end());
        }

        // Fetch rotation measure attributes (if applicable).
        if(src.getInfo().getUseRotationMeasure())
        {
          source->setRotationMeasure(src.getPolarizedFraction(),
            src.getPolarizationAngle(), src.getRotationMeasure());
        }

        componentsList[i].push_back(source);
        break;
      }
    }
  }
  sourceDB.unlock();

  std::vector<Patch::ConstPtr> patchList;
  patchList.reserve (componentsList.size());
  for (uint i=0; i<componentsList.size(); ++i) {
    if (componentsList[i].empty())
      throw Exception("No sources found for patch "
               + patchNames[i]);
    Patch::Ptr ppatch(new Patch(patchNames[i],
                                componentsList[i].begin(),
                                componentsList[i].end()));
    std::vector<BBS::PatchInfo> patchInfo(sourceDB.getPatchInfo(-1, patchNames[i]));
    assert (patchInfo.size() == 1);
    // Set the position and apparent flux of the patch.
    Position patchPosition;
    patchPosition[0] = patchInfo[0].getRa();
    patchPosition[1] = patchInfo[0].getDec();
    ppatch->setPosition (patchPosition);
    ppatch->setBrightness (patchInfo[0].apparentBrightness());
    ///    ppatch->computePosition();
    patchList.push_back (ppatch);
  }
  return patchList;
}

std::vector<std::pair<ModelComponent::ConstPtr,Patch::ConstPtr> >
makeSourceList (const std::vector<Patch::ConstPtr>& patchList) {
  std::vector<Patch::ConstPtr>::const_iterator pIter=patchList.begin();
  std::vector<Patch::ConstPtr>::const_iterator pEnd =patchList.end();

  uint nSources=0;
  for (; pIter!=pEnd; ++pIter) {
    nSources+=(*pIter)->nComponents();
  }

  std::vector<std::pair<ModelComponent::ConstPtr,Patch::ConstPtr> > sourceList;
  sourceList.reserve(nSources);

  pIter=patchList.begin();

  for (; pIter!=pEnd; ++pIter) {
    Patch::const_iterator sIter=(*pIter)->begin();
    Patch::const_iterator sEnd =(*pIter)->end();
    for (; sIter!=sEnd; ++sIter) {
      sourceList.push_back(make_pair(*sIter,*pIter));
    }
  }

  return sourceList;
}


std::vector<Patch::ConstPtr> makeOnePatchPerComponent(
    const std::vector<Patch::ConstPtr>& patchList) {
    size_t numComponents=0;
    std::vector<Patch::ConstPtr>::const_iterator patchIt;

    for (patchIt=patchList.begin();patchIt!=patchList.end();++patchIt) {
        numComponents+=(*patchIt)->nComponents();
    }

    std::vector<Patch::ConstPtr> largePatchList;
    largePatchList.reserve(numComponents);

    for (patchIt=patchList.begin();patchIt!=patchList.end();++patchIt) {
        Patch::const_iterator compIt;

        size_t compNum=0;
        for (compIt=(*patchIt)->begin();compIt!=(*patchIt)->end();++compIt) {
            // convert compNum to string (blegh)
            std::stringstream ss;
            ss<<compNum;

            Patch::Ptr ppatch(new Patch((*patchIt)->name()+"_"+ss.str(),
                                        compIt,
                                        compIt+1));
            ppatch->setPosition((*compIt)->position());
            largePatchList.push_back(ppatch);
            compNum++;
        }
    }

    return largePatchList;
}


std::vector<string> makePatchList(SourceDB &sourceDB, std::vector<string> patterns)
{
    if(patterns.empty())
    {
        patterns.push_back("*");
    }

    std::set<string> patches;
    std::vector<string>::iterator it = patterns.begin();
    while(it != patterns.end())
    {
        if(!it->empty() && (*it)[0] == '@')
        {
            patches.insert(*it);
            it = patterns.erase(it);
        }
        else
        {
            std::vector<string> match(sourceDB.getPatches(-1, *it));
            patches.insert(match.begin(), match.end());
            ++it;
        }
    }

    return std::vector<string>(patches.begin(), patches.end());
}


bool checkPolarized(SourceDB &sourceDB,
                    const std::vector<string> &patchNames,
                    uint nModel)
{
  bool polarized = false;

  // Loop over all sources.
  sourceDB.lock();
  sourceDB.rewind();
  SourceData src;
  while (! sourceDB.atEnd()) {
    sourceDB.getNextSource (src);
    // Use the source if its patch matches a patch name.
    for (uint i=0; i<nModel; ++i) {
      if (src.getPatchName() == patchNames[i]) {
        // Determine whether source is unpolarized.
        if (src.getV() > 0.0 || src.getQ() > 0.0 || src.getU() > 0.0) {
          polarized = true;
          break;
        }
      }
    }
    if (polarized) {
      break;
    }
  }
  sourceDB.unlock();
  return polarized;
}

} //# namespace DPPP
} //# namespace LOFAR
