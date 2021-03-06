// DemixerNew.cc: DPPP step class to subtract A-team sources in adaptive way
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "DemixerNew.h"

#include "../base/DemixWorker.h"
#include "../base/DemixInfo.h"
#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"
#include "../base/Exceptions.h"

#include "../parmdb/ParmDB.h"
#include "../parmdb/ParmValue.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <aocommon/parallelfor.h>

#include <casacore/casa/Arrays/ArrayPartMath.h>

#include <iomanip>

using casacore::Matrix;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

namespace dp3 {
namespace steps {

using dp3::common::operator<<;

DemixerNew::DemixerNew(InputStep* input, const common::ParameterSet& parset,
                       const string& prefix)
    : itsInput(input),
      itsName(prefix),
      itsDemixInfo(parset, prefix),
      itsInstrumentName(
          parset.getString(prefix + "instrumentmodel", "instrument")),
      itsFilter(input, itsDemixInfo.selBL()),
      itsNTime(0),
      itsNTimeOut(0),
      itsNChunk(0) {
  if (itsInstrumentName.empty())
    throw Exception("An empty name is given for the instrument model");
}

void DemixerNew::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  // Update the info of this object.
  info().setNeedVisData();
  info().setWriteData();
  info().setWriteFlags();
  // Handle possible data selection.
  itsFilter.updateInfo(getInfo());
  // Update itsDemixInfo and info().
  itsDemixInfo.update(itsFilter.getInfo(), info(), getInfo().nThreads());
  // Size the buffers.
  itsBufIn.resize(itsDemixInfo.ntimeChunk() * itsDemixInfo.chunkSize());
  itsBufOut.resize(itsDemixInfo.ntimeChunk() * itsDemixInfo.ntimeOutSubtr());
  itsSolutions.resize(itsDemixInfo.ntimeChunk() * itsDemixInfo.ntimeOut());
  // Create a worker per thread.
  size_t nthread = getInfo().nThreads();
  itsWorkers.reserve(nthread);
  for (size_t i = 0; i < nthread; ++i) {
    itsWorkers.emplace_back(itsInput, itsName, itsDemixInfo, infoIn, i);
  }
}

void DemixerNew::show(std::ostream& os) const {
  os << "DemixerNew " << itsName << '\n';
  os << "  instrumentmodel:    " << itsInstrumentName << '\n';
  itsDemixInfo.show(os);
  if (itsDemixInfo.selBL().hasSelection()) {
    os << "    demixing " << itsFilter.getInfo().nbaselines() << " out of "
       << getInfo().nbaselines() << " baselines   ("
       << itsFilter.getInfo().antennaUsed().size() << " out of "
       << getInfo().antennaUsed().size() << " stations)" << '\n';
  }
}

void DemixerNew::showCounts(std::ostream& os) const {
  os << '\n' << "Statistics for SmartDemixer " << itsName;
  os << '\n' << "===========================" << '\n';
  // Add the statistics of all workers.
  unsigned int nsolves = 0;
  unsigned int nconverged = 0;
  unsigned int niter = 0;
  unsigned int nnodemix = 0;
  unsigned int nincludeStrong = 0;
  unsigned int nincludeClose = 0;
  unsigned int nignore = 0;
  unsigned int ndeproject = 0;
  // Sum the counts from all workers.
  unsigned int ndir = itsDemixInfo.ateamList().size();
  casacore::Vector<unsigned int> nsources(ndir, 0);
  casacore::Vector<unsigned int> nstation(itsDemixInfo.nstation(), 0);
  Matrix<unsigned int> statsrcs(ndir, nstation.size(), 0);
  Matrix<double> amplSubtrMean(itsDemixInfo.nbl(), ndir + 1, 0.);
  Matrix<double> amplSubtrM2(itsDemixInfo.nbl(), ndir + 1, 0.);
  Matrix<size_t> amplSubtrNr(itsDemixInfo.nbl(), ndir + 1, 0);
  double* amplmean = amplSubtrMean.data();
  double* amplm2 = amplSubtrM2.data();
  size_t* amplnr = amplSubtrNr.data();
  for (size_t i = 0; i < itsWorkers.size(); ++i) {
    nsolves += itsWorkers[i].nSolves();
    nconverged += itsWorkers[i].nConverged();
    niter += itsWorkers[i].nIterations();
    nnodemix += itsWorkers[i].nNoDemix();
    nincludeStrong += itsWorkers[i].nIncludeStrongTarget();
    nincludeClose += itsWorkers[i].nIncludeCloseTarget();
    nignore += itsWorkers[i].nIgnoreTarget();
    ndeproject += itsWorkers[i].nDeprojectTarget();
    nsources += itsWorkers[i].nsourcesDemixed();
    nstation += itsWorkers[i].nstationsDemixed();
    statsrcs += itsWorkers[i].statSourceDemixed();
    const double* partmean = itsWorkers[i].amplSubtrMean().data();
    const double* partm2 = itsWorkers[i].amplSubtrM2().data();
    const size_t* partnr = itsWorkers[i].amplSubtrNr().data();
    for (unsigned int j = 0; j < amplSubtrNr.size(); ++j) {
      // Calculate overall mean and m2 for each baseline/direction.
      addMeanM2(amplmean[j], amplm2[j], amplnr[j], partmean[j], partm2[j],
                partnr[j]);
    }
  }
  unsigned int ntimes =
      (nnodemix + nignore + ndeproject + nincludeStrong + nincludeClose);
  // Show statistics.
  showStat(os, nconverged, nsolves, "Converged solves:      ", "cells");
  os << "  Average nr of iterations in converged solves: "
     << (unsigned int)(double(niter) / nconverged + 0.5) << '\n';
  showStat(os, nnodemix, ntimes, "No demixing:           ", "times");
  showStat(os, nignore, ntimes, "Target ignored:        ", "times");
  showStat(os, ndeproject, ntimes, "Target deprojected:    ", "times");
  showStat(os, nincludeStrong, ntimes, "Strong target included:", "times");
  showStat(os, nincludeClose, ntimes, "Close target included: ", "times");
  // Show how often a source/station is demixed.
  os << '\n' << "Percentage of times a station/source is demixed:" << '\n';
  os << std::setw(15) << " ";
  for (size_t dr = 0; dr < ndir; ++dr) {
    os << std::setw(8) << itsDemixInfo.ateamList()[dr]->name().substr(0, 8);
  }
  os << " Overall" << '\n';
  for (size_t st = 0; st < nstation.size(); ++st) {
    os << std::setw(4) << st << ' ';
    // Show 10 characters of the source names; append with blanks as needed.
    unsigned int inx = itsFilter.getInfo().antennaUsed()[st];
    string nm = itsFilter.getInfo().antennaNames()[inx];
    os << nm.substr(0, 10);
    if (nm.size() < 10) {
      os << std::setw(10 - nm.size()) << ' ';
    }
    for (size_t dr = 0; dr < ndir; ++dr) {
      os << "  ";
      showPerc1(os, statsrcs(dr, st) / double(ntimes));
    }
    os << "  ";
    showPerc1(os, nstation(st) / double(ntimes));
    os << '\n';
  }
  os << "     Overall" << std::setw(3) << ' ';
  for (size_t dr = 0; dr < ndir; ++dr) {
    os << "  ";
    showPerc1(os, nsources[dr] / double(ntimes));
  }
  os << '\n';

  // Show percentage subtracted.
  if (itsDemixInfo.doSubtract()) {
    os << '\n'
       << "Mean/stddev percentage of subtracted Stokes I amplitude"
       << " for the middle channel" << '\n';
    os << std::setw(8) << ' ';
    for (size_t dr = 0; dr < ndir; ++dr) {
      if (nsources[dr] > 0) {
        // Print name a bit right of the center.
        const string& nm = itsDemixInfo.ateamList()[dr]->name().substr(0, 13);
        if (nm.size() > 10) {
          std::cout << std::setw(13) << nm;
        } else {
          int szws = 13 - nm.size();  // whitespace
          os << std::setw(szws / 2 + 1) << ' ';
          os << nm;
          os << std::setw(szws - szws / 2 - 1) << ' ';
        }
      }
    }
    os << std::setw(10) << "Total" << '\n';
    os << "baseline";
    for (size_t dr = 0; dr < ndir; ++dr) {
      if (nsources[dr] > 0) {
        os << "  mean stddev";
      }
    }
    os << "  mean stddev" << '\n';
    std::vector<double> totsump(ndir + 1, 0.);
    std::vector<double> totm2p(ndir + 1, 0.);
    std::vector<size_t> totnrp(ndir + 1, 0);
    for (int bl = 0; bl < amplSubtrMean.shape()[0]; ++bl) {
      // Do not show if nothing subtracted for this baseline.
      // Last entry contains the sum over all directions!!
      if (amplSubtrNr(bl, ndir)) {
        os << std::setw(4) << itsDemixInfo.getAnt1()[bl] << '-' << std::setw(2)
           << itsDemixInfo.getAnt2()[bl] << "  ";
        for (unsigned int dr = 0; dr < ndir + 1; ++dr) {
          if (dr == ndir || nsources[dr] > 0) {
            showPerc1(os, amplSubtrMean(bl, dr));
            double variance = 0;
            if (amplSubtrNr(bl, dr) > 1) {
              variance = amplSubtrM2(bl, dr) / (amplSubtrNr(bl, dr) - 1);
            }
            showPerc1(os, sqrt(variance));
            os << ' ';
            addMeanM2(totsump[dr], totm2p[dr], totnrp[dr],
                      amplSubtrMean(bl, dr), amplSubtrM2(bl, dr),
                      amplSubtrNr(bl, dr));
          }
        }
        os << '\n';
      }  // end if show
    }    // end for bl
    os << "  Total  ";
    for (unsigned int dr = 0; dr < ndir + 1; ++dr) {
      if (dr == ndir || nsources[dr] > 0) {
        double variance = 0;
        if (totnrp[dr] > 1) {
          variance = totm2p[dr] / (totnrp[dr] - 1);
        }
        showPerc1(os, totsump[dr]);
        showPerc1(os, sqrt(variance));
        os << ' ';
      }
    }
    os << '\n';
  }
}

void DemixerNew::addMeanM2(double& mean, double& m2, size_t& nr,
                           double partmean, double partm2,
                           size_t partnr) const {
  // Calculate overall mean and m2 (for stddev) from two partitions.
  // See en.wikipedia.org/wiki/Algorithms_for_calculating_variance
  if (partnr > 0) {
    double delta = partmean - mean;
    double mp = mean * nr + partmean * partnr;
    double nab = double(nr) / (nr + partnr) * partnr;
    nr += partnr;
    mean = mp / nr;
    m2 += partm2 + delta * nab * delta;
  }
}

void DemixerNew::showStat(std::ostream& os, double n, double ntot,
                          const string& str1, const string& str2) const {
  os << str1 << ' ';
  showPerc1(os, ntot == 0 ? 0 : n / ntot);
  os << "  (" << n << ' ' << str2 << " out of " << ntot << ')' << '\n';
}

void DemixerNew::showPerc1(std::ostream& os, float perc) const {
  int p = int(1000 * perc + 0.5);
  os << std::setw(3) << p / 10 << '.' << p % 10 << '%';
}

void DemixerNew::showTimings(std::ostream& os, double duration) const {
  double self = itsTimer.getElapsed();
  double demix = itsTimerDemix.getElapsed();
  double tottime = 0;
  double coatime = 0;
  double psatime = 0;
  double demtime = 0;
  double pretime = 0;
  double soltime = 0;
  double subtime = 0;
  for (unsigned int i = 0; i < itsWorkers.size(); ++i) {
    tottime += itsWorkers[i].getTotalTime();
    coatime += itsWorkers[i].getCoarseTime();
    psatime += itsWorkers[i].getPhaseShiftTime();
    demtime += itsWorkers[i].getDemixTime();
    pretime += itsWorkers[i].getPredictTime();
    soltime += itsWorkers[i].getSolveTime();
    subtime += itsWorkers[i].getSubtractTime();
  }
  os << "  ";
  FlagCounter::showPerc1(os, self, duration);
  os << " DemixerNew " << itsName << '\n';
  os << "          ";
  FlagCounter::showPerc1(os, demix, self);
  os << " of it spent in demixing the data of which" << '\n';
  os << "                ";
  FlagCounter::showPerc1(os, coatime, tottime);
  os << " in predicting coarse source models" << '\n';
  os << "                ";
  FlagCounter::showPerc1(os, psatime, tottime);
  os << " in phase shifting/averaging data" << '\n';
  os << "                ";
  FlagCounter::showPerc1(os, demtime, tottime);
  os << " in calculating decorrelation factors" << '\n';
  os << "                ";
  FlagCounter::showPerc1(os, pretime, tottime);
  os << " in predicting demix source models" << '\n';
  os << "                ";
  FlagCounter::showPerc1(os, soltime, tottime);
  os << " in solving complex gains" << '\n';
  os << "                ";
  FlagCounter::showPerc1(os, subtime, tottime);
  os << " in subtracting source models" << '\n';
  os << "          ";
  FlagCounter::showPerc1(os, itsTimerDump.getElapsed(), self);
  os << " of it spent in writing gain solutions to disk" << '\n';
}

bool DemixerNew::process(const DPBuffer& buf) {
  itsTimer.start();
  // Collect sufficient data buffers.
  // Make sure all required data arrays are filled in.
  DPBuffer& newBuf = itsBufIn[itsNTime];
  newBuf.copy(buf);
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

void DemixerNew::processData() {
  itsTimerDemix.start();
  // Last batch might contain fewer time slots.
  unsigned int timeWindowIn = itsDemixInfo.chunkSize();
  unsigned int timeWindowOut = itsDemixInfo.ntimeOutSubtr();
  unsigned int timeWindowSol = itsDemixInfo.ntimeOut();
  int lastChunk = (itsNTime - 1) / timeWindowIn;
  int lastNTimeIn = itsNTime - lastChunk * timeWindowIn;
  int ntimeOut = ((itsNTime + itsDemixInfo.ntimeAvgSubtr() - 1) /
                  itsDemixInfo.ntimeAvgSubtr());
  int ntimeSol =
      ((itsNTime + itsDemixInfo.ntimeAvg() - 1) / itsDemixInfo.ntimeAvg());
  aocommon::ParallelFor<int> loop(getInfo().nThreads());
  loop.Run(0, lastChunk, [&](int i, size_t thread) {
    if (i == lastChunk) {
      itsWorkers[thread].process(&(itsBufIn[i * timeWindowIn]), lastNTimeIn,
                                 &(itsBufOut[i * timeWindowOut]),
                                 &(itsSolutions[i * timeWindowSol]),
                                 itsNChunk + i);
    } else {
      itsWorkers[thread].process(&(itsBufIn[i * timeWindowIn]), timeWindowIn,
                                 &(itsBufOut[i * timeWindowOut]),
                                 &(itsSolutions[i * timeWindowSol]),
                                 itsNChunk + i);
    }
  });
  itsNChunk += lastChunk + 1;
  itsTimerDemix.stop();
  // Write the solutions into the instrument ParmDB.
  // Let the next steps process the results.
  // This could be done in parallel.
  ///#pragma omp parallel for num_thread(2)
  for (int i = 0; i < 2; ++i) {
    if (i == 0) {
      itsTimerDump.start();
      double startTime =
          (itsBufIn[0].getTime() + 0.5 * itsBufIn[0].getExposure());
      writeSolutions(startTime, ntimeSol);
      itsTimerDump.stop();
    } else {
      itsTimer.stop();
      itsTimerNext.start();
      for (int j = 0; j < ntimeOut; ++j) {
        getNextStep()->process(itsBufOut[j]);
        itsNTimeOut++;
      }
      itsTimerNext.stop();
      itsTimer.start();
    }
  }
}

void DemixerNew::finish() {
  std::cerr << "  " << itsNTime << " time slots to finish in SmartDemixer ..."
            << '\n';
  itsTimer.start();
  // Process remaining entries.
  if (itsNTime > 0) {
    processData();
  }
  itsTimer.stop();
  // Let the next steps finish.
  getNextStep()->finish();
}

void DemixerNew::writeSolutions(double startTime, int ntime) {
  if (itsDemixInfo.verbose() > 12) {
    for (int i = 0; i < ntime; ++i) {
      std::cout << "solution " << i << '\n';
      const double* sol = &(itsSolutions[i][0]);
      for (size_t dr = 0;
           dr < itsSolutions[i].size() / (8 * itsDemixInfo.nstation()); ++dr) {
        for (size_t st = 0; st < itsDemixInfo.nstation(); ++st) {
          std::cout << dr << ',' << st << ' ';
          common::print(std::cout, sol, sol + 8);
          std::cout << '\n';
          sol += 8;
        }
      }
    }
  }
  /// todo: skip solutions that are all 0.
  // Construct solution grid.
  const casacore::Vector<double>& freq = getInfo().chanFreqs();
  const casacore::Vector<double>& freqWidth = getInfo().chanWidths();
  parmdb::Axis::ShPtr freqAxis(
      new parmdb::RegularAxis(freq[0] - freqWidth[0] * 0.5, freqWidth[0], 1));
  parmdb::Axis::ShPtr timeAxis(new parmdb::RegularAxis(
      startTime, itsDemixInfo.timeIntervalAvg(), ntime));
  parmdb::Grid solGrid(freqAxis, timeAxis);
  // Create domain grid.
  parmdb::Axis::ShPtr tdomAxis(new parmdb::RegularAxis(
      startTime, itsDemixInfo.timeIntervalAvg() * ntime, 1));
  parmdb::Grid domainGrid(freqAxis, tdomAxis);

  // Open the ParmDB at the first write.
  // In that way the instrumentmodel ParmDB can be in the MS directory.
  if (!itsParmDB) {
    itsParmDB = std::shared_ptr<parmdb::ParmDB>(new parmdb::ParmDB(
        parmdb::ParmDBMeta("casa", itsInstrumentName), true));
    itsParmDB->lock();
    // Store the (freq, time) resolution of the solutions.
    std::vector<double> resolution(2);
    resolution[0] = freqWidth[0];
    resolution[1] = itsDemixInfo.timeIntervalAvg();
    itsParmDB->setDefaultSteps(resolution);
  }
  // Write the solutions per parameter.
  const char* str01[] = {"0:", "1:"};
  const char* strri[] = {"Real:", "Imag:"};
  Matrix<double> values(1, ntime);
  unsigned int seqnr = 0;
  for (size_t dr = 0; dr < itsDemixInfo.ateamList().size() + 1; ++dr) {
    for (size_t st = 0; st < itsDemixInfo.nstation(); ++st) {
      string suffix(itsDemixInfo.antennaNames()[st]);
      if (dr == itsDemixInfo.ateamList().size()) {
        suffix += ":Target";
      } else {
        suffix += ":" + itsDemixInfo.ateamList()[dr]->name();
      }
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; j++) {
          for (int k = 0; k < 2; ++k) {
            string name(string("DirectionalGain:") + str01[i] + str01[j] +
                        strri[k] + suffix);
            // Collect its solutions for all times in a single array.
            for (int ts = 0; ts < ntime; ++ts) {
              values(0, ts) = itsSolutions[ts][seqnr];
            }
            seqnr++;
            parmdb::ParmValue::ShPtr pv(new parmdb::ParmValue());
            pv->setScalars(solGrid, values);
            parmdb::ParmValueSet pvs(
                domainGrid, std::vector<parmdb::ParmValue::ShPtr>(1, pv));
            std::map<std::string, int>::const_iterator pit =
                itsParmIdMap.find(name);
            if (pit == itsParmIdMap.end()) {
              // First time, so a new nameId will be set.
              int nameId = -1;
              itsParmDB->putValues(name, nameId, pvs);
              itsParmIdMap[name] = nameId;
            } else {
              // Parm has been put before.
              int nameId = pit->second;
              itsParmDB->putValues(name, nameId, pvs);
            }
          }
        }
      }
    }
  }
}

}  // namespace steps
}  // namespace dp3
