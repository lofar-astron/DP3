//# DemixerNew.cc: DPPP step class to subtract A-team sources in adaptive way
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
//# $Id: DemixerNew.cc 24221 2013-03-12 12:24:48Z diepen $
//#
//# @author Ger van Diepen

#include "DemixerNew.h"
#include "DemixWorker.h"
#include "DemixInfo.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "Exceptions.h"

#include "../ParmDB/ParmDB.h"
#include "../ParmDB/ParmValue.h"

#include "../Common/ParallelFor.h"
#include "../Common/ParameterSet.h"
#include "../Common/StreamUtil.h"

#include <casacore/casa/Arrays/ArrayPartMath.h>

#include <iomanip>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    using DP3::operator<<;

    DemixerNew::DemixerNew (DPInput* input,
                            const ParameterSet& parset,
                            const string& prefix)
      : itsInput          (input),
        itsName           (prefix),
        itsDemixInfo      (parset, prefix),
        itsInstrumentName (parset.getString(prefix+"instrumentmodel",
                                            "instrument")),
        itsFilter         (input, itsDemixInfo.selBL()),
        itsNTime          (0),
        itsNTimeOut       (0),
        itsNChunk         (0)
    {
      if (itsInstrumentName.empty())
        throw Exception(
                 "An empty name is given for the instrument model");
    }

    void DemixerNew::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      // Update the info of this object.
      info().setNeedVisData();
      info().setWriteData();
      info().setWriteFlags();
      // Handle possible data selection.
      itsFilter.updateInfo (getInfo());
      // Update itsDemixInfo and info().
      itsDemixInfo.update (itsFilter.getInfo(), info(), getInfo().nThreads());
      // Size the buffers.
      itsBufIn.resize (itsDemixInfo.ntimeChunk() * itsDemixInfo.chunkSize());
      itsBufOut.resize(itsDemixInfo.ntimeChunk() * itsDemixInfo.ntimeOutSubtr());
      itsSolutions.resize(itsDemixInfo.ntimeChunk() * itsDemixInfo.ntimeOut());
      // Create a worker per thread.
      size_t nthread = getInfo().nThreads();
      itsWorkers.reserve (nthread);
      for (size_t i=0; i<nthread; ++i) {
        itsWorkers.emplace_back (itsInput, itsName, itsDemixInfo,
                                           infoIn, i);
      }
    }

    void DemixerNew::show (ostream& os) const
    {
      os << "DemixerNew " << itsName << endl;
      os << "  instrumentmodel:    " << itsInstrumentName << endl;
      itsDemixInfo.show (os);
      if (itsDemixInfo.selBL().hasSelection()) {
        os << "    demixing " << itsFilter.getInfo().nbaselines()
           << " out of " << getInfo().nbaselines() << " baselines   ("
           << itsFilter.getInfo().antennaUsed().size()
           << " out of " << getInfo().antennaUsed().size()
           << " stations)" << std::endl;
      }
    }

    void DemixerNew::showCounts (ostream& os) const
    {
      os << endl << "Statistics for SmartDemixer " << itsName;
      os << endl << "===========================" << endl;
      // Add the statistics of all workers.
      uint nsolves        = 0;
      uint nconverged     = 0;
      uint niter          = 0;
      uint nnodemix       = 0;
      uint nincludeStrong = 0;
      uint nincludeClose  = 0;
      uint nignore        = 0;
      uint ndeproject     = 0;
      // Sum the counts from all workers.
      uint ndir           = itsDemixInfo.ateamList().size();
      Vector<uint> nsources(ndir, 0);
      Vector<uint> nstation(itsDemixInfo.nstation(), 0);
      Matrix<uint> statsrcs(ndir, nstation.size(), 0);
      Matrix<double> amplSubtrMean (itsDemixInfo.nbl(), ndir+1, 0.);
      Matrix<double> amplSubtrM2   (itsDemixInfo.nbl(), ndir+1, 0.);
      Matrix<size_t> amplSubtrNr   (itsDemixInfo.nbl(), ndir+1, 0);
      double* amplmean = amplSubtrMean.data();
      double* amplm2   = amplSubtrM2.data();
      size_t* amplnr   = amplSubtrNr.data();
      for (size_t i=0; i<itsWorkers.size(); ++i) {
        nsolves        += itsWorkers[i].nSolves();
        nconverged     += itsWorkers[i].nConverged();
        niter          += itsWorkers[i].nIterations();
        nnodemix       += itsWorkers[i].nNoDemix();
        nincludeStrong += itsWorkers[i].nIncludeStrongTarget();
        nincludeClose  += itsWorkers[i].nIncludeCloseTarget();
        nignore        += itsWorkers[i].nIgnoreTarget();
        ndeproject     += itsWorkers[i].nDeprojectTarget();
        nsources       += itsWorkers[i].nsourcesDemixed();
        nstation       += itsWorkers[i].nstationsDemixed();
        statsrcs       += itsWorkers[i].statSourceDemixed();
        const double* partmean = itsWorkers[i].amplSubtrMean().data();
        const double* partm2   = itsWorkers[i].amplSubtrM2().data();
        const size_t* partnr   = itsWorkers[i].amplSubtrNr().data();
        for (uint j=0; j<amplSubtrNr.size(); ++j) {
          // Calculate overall mean and m2 for each baseline/direction.
          addMeanM2 (amplmean[j], amplm2[j], amplnr[j],
                     partmean[j], partm2[j], partnr[j]);
        }
      }
      uint ntimes = (nnodemix + nignore + ndeproject +
                     nincludeStrong + nincludeClose);
      // Show statistics.
      showStat (os, nconverged,    nsolves, "Converged solves:      ", "cells");
      os << "  Average nr of iterations in converged solves: "
         << uint(double(niter)/nconverged + 0.5) << endl;
      showStat (os, nnodemix,       ntimes, "No demixing:           ", "times");
      showStat (os, nignore,        ntimes, "Target ignored:        ", "times");
      showStat (os, ndeproject,     ntimes, "Target deprojected:    ", "times");
      showStat (os, nincludeStrong, ntimes, "Strong target included:", "times");
      showStat (os, nincludeClose,  ntimes, "Close target included: ", "times");
      // Show how often a source/station is demixed.
      os << endl << "Percentage of times a station/source is demixed:" << endl;
      os << std::setw(15) << " ";
      for (size_t dr=0; dr<ndir; ++dr) {
        os << std::setw(8) << itsDemixInfo.ateamList()[dr]->name().substr(0,8);
      }
      os << " Overall" << endl;
      for (size_t st=0; st<nstation.size(); ++st) {
        os << std::setw(4) << st << ' ';
        // Show 10 characters of the source names; append with blanks as needed.
        uint inx = itsFilter.getInfo().antennaUsed()[st];
        string nm = itsFilter.getInfo().antennaNames()[inx];
        os << nm.substr(0,10);
        if (nm.size() < 10) {
          os << std::setw(10 - nm.size()) << ' ';
        }
        for (size_t dr=0; dr<ndir; ++dr) {
          os << "  ";
          showPerc1 (os, statsrcs(dr,st) / double(ntimes));
        }
        os << "  ";
        showPerc1 (os, nstation(st) / double(ntimes));
        os << endl;
      }
      os << "     Overall" << std::setw(3) << ' ';
      for (size_t dr=0; dr<ndir; ++dr) {
        os << "  ";
        showPerc1 (os, nsources[dr] / double(ntimes));
      }
      os << endl;

      // Show percentage subtracted.
      if (itsDemixInfo.doSubtract()) {
        os << endl << "Mean/stddev percentage of subtracted Stokes I amplitude"
           << " for the middle channel" << endl;
        os << std::setw(8) << ' ';
        for (size_t dr=0; dr<ndir; ++dr) {
          if (nsources[dr] > 0) {
            // Print name a bit right of the center.
            const string& nm = itsDemixInfo.ateamList()[dr]->name().substr(0,13);
            if (nm.size() > 10) {
              cout << std::setw(13) << nm;
            } else {
              int szws = 13 - nm.size();    // whitespace
              os << std::setw(szws/2+1) << ' ';
              os << nm;
              os << std::setw(szws-szws/2-1) << ' ';
            }
          }
        }
        os << std::setw(10) << "Total" << endl;
        os << "baseline";
        for (size_t dr=0; dr<ndir; ++dr) {
          if (nsources[dr] > 0) {
            os << "  mean stddev";
          }
        }
        os << "  mean stddev" << endl;
        vector<double> totsump(ndir+1, 0.);
        vector<double> totm2p (ndir+1, 0.);
        vector<size_t> totnrp (ndir+1, 0);
        for (int bl=0; bl<amplSubtrMean.shape()[0]; ++bl) {
          // Do not show if nothing subtracted for this baseline.
          // Last entry contains the sum over all directions!!
          if (amplSubtrNr(bl,ndir)) {
            os << std::setw(4) << itsDemixInfo.getAnt1()[bl] << '-'
               << std::setw(2) << itsDemixInfo.getAnt2()[bl] << "  ";
            for (uint dr=0; dr<ndir+1; ++dr) {
              if (dr == ndir  ||  nsources[dr] > 0) {
                showPerc1 (os, amplSubtrMean(bl,dr));
                double variance = 0;
                if (amplSubtrNr(bl,dr) > 1) {
                  variance = amplSubtrM2(bl,dr) / (amplSubtrNr(bl,dr) - 1);
                }
                showPerc1 (os, sqrt(variance));
                os << ' ';
                addMeanM2 (totsump[dr], totm2p[dr], totnrp[dr],
                           amplSubtrMean(bl,dr), amplSubtrM2(bl,dr),
                           amplSubtrNr(bl,dr));
              }
            }
            os << endl;
          } // end if show
        } // end for bl
        os << "  Total  ";
        for (uint dr=0; dr<ndir+1; ++dr) {
          if (dr == ndir  ||  nsources[dr] > 0) {
            double variance = 0;
            if (totnrp[dr] > 1) {
              variance = totm2p[dr] / (totnrp[dr] - 1);
            }
            showPerc1 (os, totsump[dr]);
            showPerc1 (os, sqrt(variance));
            os << ' ';
          }
        }
        os << endl;
      }
    }
      
    void DemixerNew::addMeanM2 (double& mean, double& m2, size_t& nr,
                                double partmean, double partm2,
                                size_t partnr) const
    {
      // Calculate overall mean and m2 (for stddev) from two partitions.
      // See en.wikipedia.org/wiki/Algorithms_for_calculating_variance
      if (partnr > 0) {
        double delta = partmean - mean;
        double mp    = mean*nr + partmean*partnr;
        double nab   = double(nr) / (nr+partnr) * partnr;
        nr  += partnr;
        mean = mp / nr;
        m2  += partm2 + delta*nab*delta;
      }
    }

    void DemixerNew::showStat (ostream& os, double n, double ntot,
                               const string& str1, const string& str2) const
    {
      os << str1 << ' ';
      showPerc1 (os, ntot==0 ? 0 : n/ntot);
      os << "  ("<< n << ' '<< str2 << " out of " << ntot << ')' << endl;
    }

    void DemixerNew::showPerc1 (ostream& os, float perc) const
    {
      int p = int(1000*perc + 0.5);
      os << std::setw(3) << p/10 << '.' << p%10 << '%';
    }

    void DemixerNew::showTimings (std::ostream& os, double duration) const
    {
      double self  = itsTimer.getElapsed();
      double demix = itsTimerDemix.getElapsed();
      double tottime = 0;
      double coatime = 0;
      double psatime = 0;
      double demtime = 0;
      double pretime = 0;
      double soltime = 0;
      double subtime = 0;
      for (uint i=0; i<itsWorkers.size(); ++i) {
        tottime += itsWorkers[i].getTotalTime();
        coatime += itsWorkers[i].getCoarseTime();
        psatime += itsWorkers[i].getPhaseShiftTime();
        demtime += itsWorkers[i].getDemixTime();
        pretime += itsWorkers[i].getPredictTime();
        soltime += itsWorkers[i].getSolveTime();
        subtime += itsWorkers[i].getSubtractTime();
      }
      os << "  ";
      FlagCounter::showPerc1 (os, self, duration);
      os << " DemixerNew " << itsName << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, demix, self);
      os << " of it spent in demixing the data of which" << endl;
      os << "                ";
      FlagCounter::showPerc1 (os, coatime, tottime);
      os << " in predicting coarse source models" << endl;
      os << "                ";
      FlagCounter::showPerc1 (os, psatime, tottime);
      os << " in phase shifting/averaging data" << endl;
      os << "                ";
      FlagCounter::showPerc1 (os, demtime, tottime);
      os << " in calculating decorrelation factors" << endl;
      os << "                ";
      FlagCounter::showPerc1 (os, pretime, tottime);
      os << " in predicting demix source models" << endl;
      os << "                ";
      FlagCounter::showPerc1 (os, soltime, tottime);
      os << " in solving complex gains" << endl;
      os << "                ";
      FlagCounter::showPerc1 (os, subtime, tottime);
      os << " in subtracting source models" << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerDump.getElapsed(), self);
      os << " of it spent in writing gain solutions to disk" << endl;
    }

    bool DemixerNew::process (const DPBuffer& buf)
    {
      itsTimer.start();
      // Collect sufficient data buffers.
      // Make sure all required data arrays are filled in.
      DPBuffer& newBuf = itsBufIn[itsNTime];
      newBuf.copy (buf);
      itsInput->fetchUVW(buf, newBuf, itsTimer);
      itsInput->fetchWeights(buf, newBuf, itsTimer);
      itsInput->fetchFullResFlags(buf, newBuf, itsTimer);
      // Process the data if entire buffer is filled.
      if (++itsNTime >= itsBufIn.size()) {
        processData();
        itsNTime = 0;
      }
      itsTimer.stop();
      return true;
    }
      
    void DemixerNew::processData()
    {
      itsTimerDemix.start();
      // Last batch might contain fewer time slots.
      uint timeWindowIn  = itsDemixInfo.chunkSize();
      uint timeWindowOut = itsDemixInfo.ntimeOutSubtr();
      uint timeWindowSol = itsDemixInfo.ntimeOut();
      int lastChunk = (itsNTime - 1) / timeWindowIn;
      int lastNTimeIn = itsNTime - lastChunk*timeWindowIn;
      int ntimeOut = ((itsNTime + itsDemixInfo.ntimeAvgSubtr() - 1)
                      / itsDemixInfo.ntimeAvgSubtr());
      int ntimeSol = ((itsNTime + itsDemixInfo.ntimeAvg() - 1)
                      / itsDemixInfo.ntimeAvg());
      ParallelFor<int> loop(getInfo().nThreads());
      loop.Run(0, lastChunk, [&](int i, size_t thread)
      {
        if (i == lastChunk) {
          itsWorkers[thread].process
            (&(itsBufIn[i*timeWindowIn]), lastNTimeIn,
             &(itsBufOut[i*timeWindowOut]),
             &(itsSolutions[i*timeWindowSol]),
             itsNChunk+i);
        } else {
          itsWorkers[thread].process
            (&(itsBufIn[i*timeWindowIn]), timeWindowIn,
             &(itsBufOut[i*timeWindowOut]),
             &(itsSolutions[i*timeWindowSol]),
             itsNChunk+i);
        }
      });
      itsNChunk += lastChunk+1;
      itsTimerDemix.stop();
      // Write the solutions into the instrument ParmDB.
      // Let the next steps process the results.
      // This could be done in parallel.
      ///#pragma omp parallel for num_thread(2)
      for (int i=0; i<2; ++i) {
        if (i == 0) {
          itsTimerDump.start();
          double startTime = (itsBufIn[0].getTime() +
                              0.5 * itsBufIn[0].getExposure());
          writeSolutions (startTime, ntimeSol);
          itsTimerDump.stop();
        } else {
          itsTimer.stop();
          itsTimerNext.start();
          for (int j=0; j<ntimeOut; ++j) {
            getNextStep()->process (itsBufOut[j]);
            itsNTimeOut++;
          }
          itsTimerNext.stop();
          itsTimer.start();
        }
      }
    }

    void DemixerNew::finish()
    {
      cerr << "  " << itsNTime << " time slots to finish in SmartDemixer ..."
           << endl;
      itsTimer.start();
      // Process remaining entries.
      if (itsNTime > 0) {
        processData();
      }
      itsTimer.stop();
      // Let the next steps finish.
      getNextStep()->finish();
    }

    void DemixerNew::writeSolutions (double startTime, int ntime)
    {
      if (itsDemixInfo.verbose() > 12) {
        for (int i=0; i<ntime; ++i) {
          cout << "solution " << i << endl;
          const double* sol = &(itsSolutions[i][0]);
          for (size_t dr=0; dr<itsSolutions[i].size()/(8*itsDemixInfo.nstation()); ++dr) {
            for (size_t st=0; st<itsDemixInfo.nstation(); ++ st) {
              cout << dr<<','<<st<<' ';
              print (cout, sol, sol+8);
              cout << endl;
              sol += 8;
            }
          }
        }
      }
      /// todo: skip solutions that are all 0.
      // Construct solution grid.
      const Vector<double>& freq      = getInfo().chanFreqs();
      const Vector<double>& freqWidth = getInfo().chanWidths();
      BBS::Axis::ShPtr freqAxis(new BBS::RegularAxis(freq[0] - freqWidth[0]
        * 0.5, freqWidth[0], 1));
      BBS::Axis::ShPtr timeAxis(new BBS::RegularAxis
                                (startTime,
                                 itsDemixInfo.timeIntervalAvg(), ntime));
      BBS::Grid solGrid(freqAxis, timeAxis);
      // Create domain grid.
      BBS::Axis::ShPtr tdomAxis(new BBS::RegularAxis
                                (startTime,
                                 itsDemixInfo.timeIntervalAvg() * ntime, 1));
      BBS::Grid domainGrid(freqAxis, tdomAxis);

      // Open the ParmDB at the first write.
      // In that way the instrumentmodel ParmDB can be in the MS directory.
      if (! itsParmDB) {
        itsParmDB = std::shared_ptr<BBS::ParmDB>
          (new BBS::ParmDB(BBS::ParmDBMeta("casa", itsInstrumentName),
                           true));
        itsParmDB->lock();
        // Store the (freq, time) resolution of the solutions.
        vector<double> resolution(2);
        resolution[0] = freqWidth[0];
        resolution[1] = itsDemixInfo.timeIntervalAvg();
        itsParmDB->setDefaultSteps(resolution);
      }
      // Write the solutions per parameter.
      const char* str01[] = {"0:","1:"};
      const char* strri[] = {"Real:","Imag:"};
      Matrix<double> values(1, ntime);
      uint seqnr = 0;
      for (size_t dr=0; dr<itsDemixInfo.ateamList().size()+1; ++dr) {
        for (size_t st=0; st<itsDemixInfo.nstation(); ++st) {
          string suffix(itsDemixInfo.antennaNames()[st]);
          if (dr == itsDemixInfo.ateamList().size()) {
            suffix += ":Target";
          } else {
            suffix += ":" + itsDemixInfo.ateamList()[dr]->name();
          }
          for (int i=0; i<2; ++i) {
            for (int j=0; j<2; j++) {
              for (int k=0; k<2; ++k) {
                string name(string("DirectionalGain:") +
                            str01[i] + str01[j] + strri[k] + suffix);
                // Collect its solutions for all times in a single array.
                for (int ts=0; ts<ntime; ++ts) {
                  values(0, ts) = itsSolutions[ts][seqnr];
                }
                seqnr++;
                BBS::ParmValue::ShPtr pv(new BBS::ParmValue());
                pv->setScalars (solGrid, values);
                BBS::ParmValueSet pvs(domainGrid,
                                      vector<BBS::ParmValue::ShPtr>(1, pv));
                std::map<std::string,int>::const_iterator pit = itsParmIdMap.find(name);
                if (pit == itsParmIdMap.end()) {
                  // First time, so a new nameId will be set.
                  int nameId = -1;
                  itsParmDB->putValues (name, nameId, pvs);
                  itsParmIdMap[name] = nameId;
                } else {
                  // Parm has been put before.
                  int nameId = pit->second;
                  itsParmDB->putValues (name, nameId, pvs);
                }
              }
            }
          }
        }
      }
    }

} //# end namespace DPPP
} //# end namespace LOFAR
