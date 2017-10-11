//# GainCal.cc: DPPP step class to ApplyCal visibilities
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
#include <DPPP/ApplyCal.h>

#include <iostream>
#include <Common/ParameterSet.h>
#include <Common/ParameterValue.h>
#include <Common/Timer.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using namespace casacore;

namespace LOFAR {
  namespace DPPP {

    ApplyCal::ApplyCal (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix,
                      bool substep,
                      string predictDirection
                      )
      : itsIsSubstep(substep)
    {
      vector<string> subStepNames;
      ParameterValue namesPar (parset.getString(prefix + "steps", ""));

      if (namesPar.isVector()) {
        subStepNames = namesPar.getStringVector();
      } else {
        subStepNames.push_back(namesPar.getString());
      }

      vector<string>::const_iterator subStepNameIter;
      for (subStepNameIter = subStepNames.begin();
           subStepNameIter != subStepNames.end();
           ++subStepNameIter) {
        string subStepName = (*subStepNameIter);
        string subStepPrefix;
        if (subStepName.empty()) {
          // No substeps given, use parameters of this step
          subStepPrefix = prefix;
        } else {
          // Substeps given, use named parameters like applycal.applySol.parmdb
          subStepPrefix = prefix + subStepName + ".";
        }
        itsApplyCals.push_back(OneApplyCal::ShPtr(new OneApplyCal(
                input, parset, subStepPrefix, prefix, substep,
                predictDirection)));
      }

      uint numSteps = itsApplyCals.size();
      for (uint step=0; step<numSteps-1; ++step) {
        itsApplyCals[step]->setNextStep(itsApplyCals[step+1]);
      }
    }

    ApplyCal::ApplyCal()
    {}

    ApplyCal::~ApplyCal()
    {}

    void ApplyCal::setNextStep (DPStep::ShPtr nextStep)
    {
      DPStep::setNextStep(itsApplyCals[0]);
      itsApplyCals[itsApplyCals.size()-1]->setNextStep(nextStep);
    }

    void ApplyCal::show(std::ostream& os) const
    {
      // If not a substep, show will be called by DPRun,
      // through the nextStep() mechanism
      if (itsIsSubstep) {
        vector<OneApplyCal::ShPtr>::const_iterator applycalIter;

        for (applycalIter = itsApplyCals.begin();
             applycalIter != itsApplyCals.end();
             applycalIter++) {
          (*applycalIter)->show(os);
        }
      }
    }

    void ApplyCal::showTimings (std::ostream& os, double duration) const
    {
      if (itsIsSubstep) {
        vector<OneApplyCal::ShPtr>::const_iterator iter;
        for (iter = itsApplyCals.begin();
             iter != itsApplyCals.end();
             iter++) {
          (*iter)->showTimings(os, duration);
        }
      }
    }

    bool ApplyCal::process (const DPBuffer& bufin)
    {
      getNextStep()->process(bufin);
      return true;
    }


    void ApplyCal::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }

    void ApplyCal::applyDiag (const DComplex* gainA, const DComplex* gainB,
                              Complex* vis, float* weight, bool* flag,
                              uint bl, uint chan, bool updateWeights,
                              FlagCounter& flagCounter) {
      // If parameter is NaN or inf, do not apply anything and flag the data
      if (! (isFinite(gainA[0].real()) && isFinite(gainA[0].imag()) &&
             isFinite(gainB[0].real()) && isFinite(gainB[0].imag()) &&
             isFinite(gainA[1].real()) && isFinite(gainA[1].imag()) &&
             isFinite(gainB[1].real()) && isFinite(gainB[1].imag())) ) {
        // Only update flagcounter for first correlation
        if (!flag[0]) {
          flagCounter.incrChannel(chan);
          flagCounter.incrBaseline(bl);
        }
        for (uint corr=0; corr<4; ++corr) {
          flag[corr]=true;
        }
        return;
      }

      vis[0] *= gainA[0] * conj(gainB[0]);
      vis[1] *= gainA[0] * conj(gainB[1]);
      vis[2] *= gainA[1] * conj(gainB[0]);
      vis[3] *= gainA[1] * conj(gainB[1]);

      if (updateWeights) {
        weight[0] /= norm(gainA[0]) * norm(gainB[0]);
        weight[1] /= norm(gainA[0]) * norm(gainB[1]);
        weight[2] /= norm(gainA[1]) * norm(gainB[0]);
        weight[3] /= norm(gainA[1]) * norm(gainB[1]);
      }
    }

    void ApplyCal::applyScalar(const DComplex* gainA, const DComplex* gainB,
                              Complex* vis, float* weight, bool* flag,
                              uint bl, uint chan, bool updateWeights,
                              FlagCounter& flagCounter) {
      // If parameter is NaN or inf, do not apply anything and flag the data
      if (! (isFinite(gainA[0].real()) && isFinite(gainA[0].imag()) &&
             isFinite(gainB[0].real()) && isFinite(gainB[0].imag())) ) {
        // Only update flagcounter for first correlation
        if (!flag[0]) {
          flagCounter.incrChannel(chan);
          flagCounter.incrBaseline(bl);
        }
        for (uint corr=0; corr<4; ++corr) {
          flag[corr]=true;
        }
        return;
      }

      vis[0] *= gainA[0] * conj(gainB[0]);
      vis[1] *= gainA[0] * conj(gainB[0]);
      vis[2] *= gainA[0] * conj(gainB[0]);
      vis[3] *= gainA[0] * conj(gainB[0]);

      if (updateWeights) {
        weight[0] /= norm(gainA[0]) * norm(gainB[0]);
        weight[1] /= norm(gainA[0]) * norm(gainB[0]);
        weight[2] /= norm(gainA[0]) * norm(gainB[0]);
        weight[3] /= norm(gainA[0]) * norm(gainB[0]);
      }
    }

    // Inverts complex 2x2 input matrix
    void ApplyCal::invert (DComplex* v, double sigmaMMSE)
    {
      // Add the variance of the nuisance term to the elements on the diagonal.
      const double variance = sigmaMMSE * sigmaMMSE;
      DComplex v0 = v[0] + variance;
      DComplex v3 = v[3] + variance;
      // Compute inverse in the usual way.
      DComplex invDet(1.0 / (v0 * v3 - v[1] * v[2]));
      v[0] = v3 * invDet;
      v[2] = v[2] * -invDet;
      v[1] = v[1] * -invDet;
      v[3] = v0 * invDet;
    }

    void ApplyCal::applyFull (const DComplex* gainA, const DComplex* gainB,
                              Complex* vis, float* weight, bool* flag,
                              uint bl, uint chan, bool doUpdateWeights,
                              FlagCounter& flagCounter) {
      DComplex gainAxvis[4];

      // If parameter is NaN or inf, do not apply anything and flag the data
      bool anyinfnan = false;
      for (uint corr=0; corr<4; ++corr) {
        if (! (isFinite(gainA[corr].real()) && isFinite(gainA[corr].imag()) &&
               isFinite(gainB[corr].real()) && isFinite(gainB[corr].imag())) ) {
          anyinfnan = true;
          break;
        }
      }
      if (anyinfnan) {
        // Only update flag counter for first correlation
        if (!flag[0]) {
          flagCounter.incrChannel(chan);
          flagCounter.incrBaseline(bl);
        }
        for (uint corr=0; corr<4; ++corr) {
          flag[corr]=true;
        }
        return;
      }

      // gainAxvis = gainA * vis
      for (uint row=0;row<2;++row) {
        for (uint col=0;col<2;++col) {
          gainAxvis[2*row+col]=gainA[2*row+0] * DComplex(vis[2*0+col]) +
                               gainA[2*row+1] * DComplex(vis[2*1+col]);
        }
      }

      // vis = gainAxvis * gainB^H
      for (uint row=0;row<2;++row) {
        for (uint col=0;col<2;++col) {
          vis[2*row+col]=gainAxvis[2*row+0] * conj(gainB[2*col+0])+
                         gainAxvis[2*row+1] * conj(gainB[2*col+1]);
        }
      }

      if (doUpdateWeights) {
        applyWeights(gainA, gainB, weight);
      }
    }

    void ApplyCal::applyWeights(const DComplex* gainA,
                                const DComplex* gainB,
                                float* weight) {
      float cov[4], normGainA[4], normGainB[4];
      for (uint i=0;i<4;++i) {
        cov[i]=1./weight[i];
        normGainA[i]=norm(gainA[i]);
        normGainB[i]=norm(gainB[i]);
      }

      weight[0]=cov[0]*(normGainA[0]*normGainB[0])
                     +cov[1]*(normGainA[0]*normGainB[1])
                     +cov[2]*(normGainA[1]*normGainB[0])
                     +cov[3]*(normGainA[1]*normGainB[1]);
      weight[0]=1./weight[0];

      weight[1]=cov[0]*(normGainA[0]*normGainB[2])
                     +cov[1]*(normGainA[0]*normGainB[3])
                     +cov[2]*(normGainA[1]*normGainB[2])
                     +cov[3]*(normGainA[1]*normGainB[3]);
      weight[1]=1./weight[1];

      weight[2]=cov[0]*(normGainA[2]*normGainB[0])
                     +cov[1]*(normGainA[2]*normGainB[1])
                     +cov[2]*(normGainA[3]*normGainB[0])
                     +cov[3]*(normGainA[3]*normGainB[1]);
      weight[2]=1./weight[2];

      weight[3]=cov[0]*(normGainA[2]*normGainB[2])
                     +cov[1]*(normGainA[2]*normGainB[3])
                     +cov[2]*(normGainA[3]*normGainB[2])
                     +cov[3]*(normGainA[3]*normGainB[3]);
      weight[3]=1./weight[3];
    }

  } //# end namespace
}
