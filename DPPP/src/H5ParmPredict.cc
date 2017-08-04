//# GainCal.cc: DPPP step class to H5ParmPredict visibilities
//# Copyright (C) 2013
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
//# $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/H5ParmPredict.h>

#include <iostream>
#include <Common/ParameterSet.h>
#include <Common/Timer.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <Common/StreamUtil.h>
#include <Common/StringUtil.h>

using namespace casacore;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    H5ParmPredict::H5ParmPredict (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix):
                          itsInput(input),
                          itsH5ParmName(parset.getString(prefix+"parmdb")),
                          itsDirections(parset.getStringVector(
                              prefix+"directions", vector<string> ()))
    {
      H5Parm h5parm = H5Parm(itsH5ParmName, false);
      H5Parm::SolTab soltab = h5parm.getSolTab(parset.getString(prefix+"correction"));

      vector<string> h5directions = soltab.getStringAxis("dir");

      if (itsDirections.empty()) {
        itsDirections = h5directions;
      } else {
        for (vector<string>::iterator it = itsDirections.begin();
          // Check that all specified directions are in the h5parm
          it != itsDirections.end(); ++it) {
          if (find(h5directions.begin(), h5directions.end(), *it) ==
              h5directions.end()) {
            THROW(Exception, "Direction "<<*it<<" not found in "<<itsH5ParmName);
          }
        }
      }

      for (uint i=0; i<itsDirections.size(); ++i) {
        string directionStr = itsDirections[i];
        vector<string> directionVec; // each direction should be like '[patch1,patch2]'
        ASSERT(directionStr.size()>2 && directionStr[0]=='[' &&
               directionStr[directionStr.size()-1]==']');
        directionVec = StringUtil::tokenize(directionStr.substr(1, directionStr.size()-2), ",");
        Predict* predictStep = new Predict(input, parset, prefix, directionVec, false);
        predictStep->setApplyCal(input, parset, prefix);
        itsPredictSteps.push_back(Predict::ShPtr(predictStep));
        if (i>0) {
          itsPredictSteps[i-1]->setNextStep(itsPredictSteps[i]);
        }
      }

      itsResultStep=new ResultStep();
      itsPredictSteps[itsPredictSteps.size()-1]->setNextStep(DPStep::ShPtr(itsResultStep));
    }

    H5ParmPredict::~H5ParmPredict()
    {}

    void H5ParmPredict::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();

      for (auto predictstep: itsPredictSteps) {
        predictstep->updateInfo(infoIn);
      }
    }

    void H5ParmPredict::show (std::ostream& os) const
    {
      os << "H5ParmPredict " << itsName << endl;
      os << "  H5Parm:     " << itsH5ParmName << endl;
      os << "  directions: " << itsDirections << endl;
    }

    void H5ParmPredict::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " H5ParmPredict " << itsName << endl;
    }

    bool H5ParmPredict::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      itsBuffer.copy (bufin);
      itsInput->fetchUVW(bufin, itsBuffer, itsTimer);
      itsInput->fetchWeights(bufin, itsBuffer, itsTimer);

      itsPredictSteps[0]->process(itsBuffer);
      itsBuffer = itsResultStep->get();

      itsTimer.stop();
      getNextStep()->process(itsBuffer);
      return false;
    }


    void H5ParmPredict::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } //# end namespace
}
