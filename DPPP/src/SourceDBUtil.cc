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

#include <lofar_config.h>
#include <DPPP/SourceDBUtil.h>
#include <DPPP/PointSource.h>
#include <DPPP/GaussianSource.h>
#include <ParmDB/SourceDB.h>
#include <Common/LofarLogger.h>
#include <Common/lofar_vector.h>

namespace LOFAR
{
namespace DPPP
{
using BBS::SourceDB;
using BBS::SourceData;
using BBS::SourceInfo;


vector<Patch::ConstPtr> makePatches(SourceDB &sourceDB,
                                    const vector<string> &patchNames,
                                    uint nModel)
{
  // Create a component list for each patch name.
  vector<vector<ModelComponent::Ptr> > componentsList(nModel);

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
        ASSERT (src.getInfo().getRefType() == "J2000");
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

                const double deg2rad = (casa::C::pi / 180.0);
                gauss->setPositionAngle(src.getOrientation() * deg2rad);

                const double arcsec2rad = (casa::C::pi / 3600.0) / 180.0;
                gauss->setMajorAxis(src.getMajorAxis() * arcsec2rad);
                gauss->setMinorAxis(src.getMinorAxis() * arcsec2rad);
                source = gauss;
            }
            break;

        default:
            {
                ASSERTSTR(false, "Only point sources and Gaussian sources are"
                    " supported at this time.");
            }
        }

        // Fetch spectral index attributes (if applicable).
        if (src.getSpectralIndex().size() > 0) {
          source->setSpectralIndex(src.getInfo().getSpectralIndexRefFreq(),
                                   src.getSpectralIndex().begin(),
                                   src.getSpectralIndex().end());
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

  vector<Patch::ConstPtr> patchList;
  patchList.reserve (componentsList.size());
  for (uint i=0; i<componentsList.size(); ++i) {
    ASSERTSTR (!componentsList[i].empty(), "No sources found for patch "
               << patchNames[i]);
    patchList.push_back (Patch::Ptr (new Patch(patchNames[i],
                                               componentsList[i].begin(),
                                               componentsList[i].end())));
  }
  return patchList;
}


} //# namespace DPPP
} //# namespace LOFAR
