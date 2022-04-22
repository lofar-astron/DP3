// OnePredict.cc: DPPP step class that predicts visibilities.
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "OnePredict.h"
#include "ApplyBeam.h"

#include <iostream>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/StreamUtil.h"

#include "../parmdb/ParmDBMeta.h"
#include "../parmdb/PatchInfo.h"
#include "../parmdb/SkymodelToSourceDB.h"

#include "../base/DPInfo.h"
#include "../base/Exceptions.h"
#include "../base/FlagCounter.h"
#include "../base/Simulate.h"
#include "../base/Simulator.h"
#include "../base/Stokes.h"
#include "../base/PointSource.h"
#include "../base/GaussianSource.h"

#include "../parmdb/SourceDB.h"

#include <aocommon/threadpool.h>
#include <aocommon/parallelfor.h>
#include <aocommon/barrier.h>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/tables/Tables/RefRows.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <numeric>
#include <algorithm>

#include <boost/algorithm/string/case_conv.hpp>

using casacore::Cube;
using casacore::MDirection;
using casacore::MEpoch;
using casacore::MVDirection;
using casacore::MVEpoch;
using casacore::Quantum;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

OnePredict::OnePredict(InputStep* input, const common::ParameterSet& parset,
                       const string& prefix,
                       const std::vector<string>& source_patterns)
    : thread_pool_(nullptr), measures_mutex_(nullptr) {
  std::vector<std::string> copied_patterns = source_patterns;
  if (source_patterns.empty()) {
    copied_patterns =
        parset.getStringVector(prefix + "sources", std::vector<std::string>());
  }
  init(input, parset, prefix, copied_patterns);
}

void OnePredict::init(InputStep* input, const common::ParameterSet& parset,
                      const string& prefix,
                      const std::vector<string>& sourcePatterns) {
  input_ = input;
  name_ = prefix;
  source_db_name_ = parset.getString(prefix + "sourcedb");
  correct_freq_smearing_ =
      parset.getBool(prefix + "correctfreqsmearing", false);
  SetOperation(parset.getString(prefix + "operation", "replace"));
  apply_beam_ = parset.getBool(prefix + "usebeammodel", false);
  thread_over_baselines_ = parset.getBool(prefix + "parallelbaselines", false);
  debug_level_ = parset.getInt(prefix + "debuglevel", 0);
  patch_list_.clear();

  // Save directions specifications to pass to applycal
  std::stringstream ss;
  ss << sourcePatterns;
  direction_str_ = ss.str();

  base::SourceDB source_db{source_db_name_, sourcePatterns,
                           base::SourceDB::FilterMode::kPattern};
  try {
    patch_list_ = source_db.MakePatchList();
    if (patch_list_.empty()) {
      throw Exception("Couldn't find patch for direction " + direction_str_);
    }
  } catch (std::exception& exception) {
    throw std::runtime_error(std::string("Something went wrong while reading "
                                         "the source model. The error was: ") +
                             exception.what());
  }

  if (apply_beam_) {
    use_channel_freq_ = parset.getBool(prefix + "usechannelfreq", true);
    one_beam_per_patch_ = parset.getBool(prefix + "onebeamperpatch", false);
    beam_proximity_limit_ =
        parset.getDouble(prefix + "beamproximitylimit", 60.0) *
        (M_PI / (180.0 * 60.0 * 60.0));

    beam_mode_ = everybeam::ParseCorrectionMode(
        parset.getString(prefix + "beammode", "default"));

    string element_model = boost::to_lower_copy(
        parset.getString(prefix + "elementmodel", "hamaker"));
    if (element_model == "hamaker") {
      element_response_model_ = everybeam::ElementResponseModel::kHamaker;
    } else if (element_model == "lobes") {
      element_response_model_ = everybeam::ElementResponseModel::kLOBES;
    } else if (element_model == "oskar") {
      element_response_model_ =
          everybeam::ElementResponseModel::kOSKARSphericalWave;
    } else if (element_model == "oskardipole") {
      element_response_model_ = everybeam::ElementResponseModel::kOSKARDipole;
    } else {
      throw Exception(
          "Elementmodel should be HAMAKER, LOBES, OSKAR or OSKARDIPOLE");
    }

    // By default, a source model has each direction in one patch. Therefore,
    // if one-beam-per-patch is requested, we don't have to do anything.
    if (!one_beam_per_patch_) {
      if (beam_proximity_limit_ > 0.0) {
        // Rework patch list to cluster proximate sources
        patch_list_ =
            clusterProximateSources(patch_list_, beam_proximity_limit_);
      } else {
        // Rework patch list to contain a patch for every source
        patch_list_ = makeOnePatchPerComponent(patch_list_);
      }
    }
  }

  // If called from h5parmpredict, applycal gets set by that step,
  // so must not be read from parset
  if (parset.isDefined(prefix + "applycal.parmdb") ||
      parset.isDefined(prefix + "applycal.steps")) {
    SetApplyCal(input, parset, prefix + "applycal.");
  } else {
    do_apply_cal_ = false;
  }

  source_list_ = makeSourceList(patch_list_);

  // Determine whether any sources are polarized. If not, enable Stokes-I-
  // only mode (note that this mode cannot be used with apply_beam_)
  if (apply_beam_ && beam_mode_ != everybeam::CorrectionMode::kArrayFactor) {
    stokes_i_only_ = false;
  } else {
    stokes_i_only_ = !source_db.CheckPolarized();
  }
}

void OnePredict::SetApplyCal(InputStep* input,
                             const common::ParameterSet& parset,
                             const string& prefix) {
  do_apply_cal_ = true;
  apply_cal_step_ = ApplyCal(input, parset, prefix, true, direction_str_);
  if (operation_ != "replace" &&
      parset.getBool(prefix + "applycal.updateweights", false))
    throw std::invalid_argument(
        "Weights cannot be updated when operation is not replace");
  result_step_ = std::make_shared<ResultStep>();
  apply_cal_step_.setNextStep(result_step_);
}

OnePredict::~OnePredict() {}

void OnePredict::initializeThreadData() {
  const size_t nBl = info().nbaselines();
  const size_t nSt = info().nantenna();
  const size_t nCh = info().nchan();
  const size_t nCr = stokes_i_only_ ? 1 : info().ncorr();
  const size_t nThreads = getInfo().nThreads();

  station_uvw_.resize(3, nSt);

  std::vector<std::array<double, 3>> antenna_pos(info().antennaPos().size());
  for (unsigned int i = 0; i < info().antennaPos().size(); ++i) {
    casacore::Quantum<casacore::Vector<double>> pos =
        info().antennaPos()[i].get("m");
    antenna_pos[i][0] = pos.getValue()[0];
    antenna_pos[i][1] = pos.getValue()[1];
    antenna_pos[i][2] = pos.getValue()[2];
  }

  uvw_split_index_ = base::nsetupSplitUVW(info().nantenna(), info().getAnt1(),
                                          info().getAnt2(), antenna_pos);

  if (!predict_buffer_) {
    predict_buffer_ = std::make_shared<base::PredictBuffer>();
  }
  if (apply_beam_ && predict_buffer_->GetStationList().empty()) {
    telescope_ =
        input_->GetTelescope(element_response_model_, use_channel_freq_);
  }
  predict_buffer_->resize(nThreads, nCr, nCh, nBl, nSt, apply_beam_);
  // Create the Measure ITRF conversion info given the array position.
  // The time and direction are filled in later.
  meas_convertors_.resize(nThreads);
  meas_frame_.resize(nThreads);

  for (size_t thread = 0; thread < nThreads; ++thread) {
    const bool need_meas_converters = moving_phase_ref_ || apply_beam_;
    if (need_meas_converters) {
      // Prepare measures converters
      meas_frame_[thread].set(info().arrayPosCopy());
      meas_frame_[thread].set(
          MEpoch(MVEpoch(info().startTime() / 86400), MEpoch::UTC));
      meas_convertors_[thread].set(
          MDirection::J2000,
          MDirection::Ref(MDirection::ITRF, meas_frame_[thread]));
    }
  }
}

void OnePredict::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
  if (operation_ == "replace")
    info().setBeamCorrectionMode(everybeam::CorrectionMode::kNone);

  const size_t nBl = info().nbaselines();
  for (size_t i = 0; i != nBl; ++i) {
    baselines_.emplace_back(info().getAnt1()[i], info().getAnt2()[i]);
  }

  try {
    MDirection dirJ2000(
        MDirection::Convert(infoIn.phaseCenter(), MDirection::J2000)());
    Quantum<casacore::Vector<double>> angles = dirJ2000.getAngle();
    moving_phase_ref_ = false;
    phase_ref_ =
        base::Direction(angles.getBaseValue()[0], angles.getBaseValue()[1]);
  } catch (casacore::AipsError&) {
    // Phase direction (in J2000) is time dependent
    moving_phase_ref_ = true;
  }

  initializeThreadData();

  if (do_apply_cal_) {
    info() = apply_cal_step_.setInfo(info());
  }
}

base::Direction OnePredict::GetFirstDirection() const {
  return patch_list_.front()->direction();
}

void OnePredict::SetOperation(const std::string& operation) {
  operation_ = operation;
  if (operation_ != "replace" && operation_ != "add" &&
      operation_ != "subtract")
    throw std::invalid_argument(
        "Operation must be 'replace', 'add' or 'subtract'.");
}

void OnePredict::show(std::ostream& os) const {
  os << "OnePredict " << name_ << '\n';
  os << "  sourcedb:           " << source_db_name_ << '\n';
  os << "   number of patches: " << patch_list_.size() << '\n';
  os << "   number of sources: " << source_list_.size() << '\n';
  os << "   all unpolarized:   " << std::boolalpha << stokes_i_only_ << '\n';
  os << "   correct freq smearing: " << std::boolalpha << correct_freq_smearing_
     << '\n';
  os << "  apply beam:         " << std::boolalpha << apply_beam_ << '\n';
  if (apply_beam_) {
    os << "   mode:              " << everybeam::ToString(beam_mode_);
    os << '\n';
    os << "   use channelfreq:   " << std::boolalpha << use_channel_freq_
       << '\n';
    os << "   one beam per patch:" << std::boolalpha << one_beam_per_patch_
       << '\n';
    os << "   beam proximity lim:"
       << (beam_proximity_limit_ * (180.0 * 60.0 * 60.0) / M_PI) << " arcsec\n";
  }
  os << "  operation:          " << operation_ << '\n';
  os << "  threads:            " << getInfo().nThreads() << '\n';
  if (do_apply_cal_) {
    apply_cal_step_.show(os);
  }
}

void OnePredict::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " OnePredict " << name_ << '\n';

  /*
   * The timer_ measures the time in a single thread. Both predict_time_ and
   * apply_beam_time_ are the sum of time in multiple threads. This makes it
   * hard to determine the exact time spent in these phases. Instead it shows
   * the percentage spent in these two parts.
   */
  const int64_t time{predict_time_ + apply_beam_time_};
  os << "          ";
  base::FlagCounter::showPerc1(os, predict_time_, time);
  os << " of it spent in predict" << '\n';

  os << "          ";
  base::FlagCounter::showPerc1(os, apply_beam_time_, time);
  os << " of it spent in apply beam" << '\n';
}

bool OnePredict::process(const DPBuffer& bufin) {
  timer_.start();
  DPBuffer scratch_buffer;
  scratch_buffer.copy(bufin);
  input_->fetchUVW(bufin, scratch_buffer, timer_);
  input_->fetchWeights(bufin, scratch_buffer, timer_);

  // Determine the various sizes.
  // const size_t nDr = patch_list_.size();
  const size_t nSt = info().nantenna();
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nCr = info().ncorr();
  const size_t nBeamValues = stokes_i_only_ ? nBl * nCh : nBl * nCh * nCr;

  base::nsplitUVW(uvw_split_index_, baselines_, scratch_buffer.getUVW(),
                  station_uvw_);

  double time = scratch_buffer.getTime();
  // Set up directions for beam evaluation
  everybeam::vector3r_t refdir, tiledir;

  const bool need_meas_converters = moving_phase_ref_ || apply_beam_;
  if (need_meas_converters) {
    // Because multiple predict steps might be predicting simultaneously, and
    // Casacore is not thread safe, this needs synchronization.
    std::unique_lock<std::mutex> lock;
    if (measures_mutex_ != nullptr)
      lock = std::unique_lock<std::mutex>(*measures_mutex_);
    for (size_t thread = 0; thread != getInfo().nThreads(); ++thread) {
      meas_frame_[thread].resetEpoch(
          MEpoch(MVEpoch(time / 86400), MEpoch::UTC));
      // Do a conversion on all threads
      refdir = dir2Itrf(info().delayCenter(), meas_convertors_[thread]);
      tiledir = dir2Itrf(info().tileBeamDir(), meas_convertors_[thread]);
    }
  }

  if (moving_phase_ref_) {
    // Convert phase reference to J2000
    MDirection dirJ2000(MDirection::Convert(
        info().phaseCenter(),
        MDirection::Ref(MDirection::J2000, meas_frame_[0]))());
    Quantum<casacore::Vector<double>> angles = dirJ2000.getAngle();
    phase_ref_ =
        base::Direction(angles.getBaseValue()[0], angles.getBaseValue()[1]);
  }

  std::unique_ptr<aocommon::ThreadPool> localThreadPool;
  aocommon::ThreadPool* pool = thread_pool_;
  if (pool == nullptr) {
    // If no ThreadPool was specified, we create a temporary one just
    // for executation of this part.
    localThreadPool = std::make_unique<aocommon::ThreadPool>(info().nThreads());
    pool = localThreadPool.get();
  } else {
    if (pool->NThreads() != info().nThreads())
      throw std::runtime_error(
          "Thread pool has inconsistent number of threads!");
  }
  std::vector<base::Simulator> simulators;
  simulators.reserve(pool->NThreads());

  std::vector<std::pair<size_t, size_t>> baseline_range;
  std::vector<Cube<dcomplex>> sim_buffer;
  std::vector<std::vector<std::pair<size_t, size_t>>> baselines_split;
  std::vector<std::pair<size_t, size_t>> station_range;
  const size_t actual_nCr = (stokes_i_only_ ? 1 : nCr);
  size_t n_threads = pool->NThreads();
  if (thread_over_baselines_) {
    // make sure each thread has enough baselines to work
    if (nBl < n_threads) {
      n_threads = nBl;
    }
    baseline_range.resize(n_threads);
    sim_buffer.resize(n_threads);
    baselines_split.resize(n_threads);
    if (apply_beam_) {
      station_range.resize(n_threads);
    }

    const size_t baselines_thread = (nBl + n_threads - 1) / n_threads;
    size_t end_index = 0;
    size_t start_index = 0;
    for (size_t thread_index = 0;
         thread_index != n_threads && end_index < baselines_.size();
         ++thread_index) {
      auto start_bline =
          std::next(baselines_.cbegin(), thread_index * baselines_thread);
      end_index = (thread_index * baselines_thread + baselines_thread >
                           baselines_.size()
                       ? baselines_.size()
                       : thread_index * baselines_thread + baselines_thread);
      baseline_range[thread_index] = std::make_pair(start_index, end_index);
      auto end_bline = std::next(baselines_.cbegin(), end_index);
      baselines_split[thread_index].resize(end_index -
                                           thread_index * baselines_thread);
      sim_buffer[thread_index].resize(
          actual_nCr, nCh, end_index - thread_index * baselines_thread);
      start_index += baselines_split[thread_index].size();
      std::copy(start_bline, end_bline, baselines_split[thread_index].begin());
    }
    // find min,max station indices for this thread
    if (apply_beam_) {
      const size_t stations_thread = (nSt + n_threads - 1) / n_threads;
      for (size_t thread_index = 0; thread_index != n_threads; ++thread_index) {
        const size_t station_start = thread_index * stations_thread;
        const size_t station_end = station_start + stations_thread < nSt
                                       ? station_start + stations_thread
                                       : nSt;
        if (station_start < nSt) {
          station_range[thread_index] =
              std::make_pair(station_start, station_end);
        } else {
          // fill an invalid station range
          // so that station_start<nSt for valid range
          station_range[thread_index] = std::make_pair(nSt + 1, nSt + 1);
        }
      }
    }
    for (size_t thread_index = 0; thread_index != n_threads; ++thread_index) {
      predict_buffer_->GetModel(thread_index) = dcomplex();
      if (apply_beam_)
        predict_buffer_->GetPatchModel(thread_index) = dcomplex();

      // When applying beam, simulate into patch vector
      Cube<dcomplex>& simulatedest = sim_buffer[thread_index];
      simulatedest = dcomplex();
      simulators.emplace_back(phase_ref_, nSt, baselines_split[thread_index],
                              info().chanFreqs(), info().chanWidths(),
                              station_uvw_, simulatedest,
                              correct_freq_smearing_, stokes_i_only_);
    }
  } else {
    for (size_t thread_index = 0; thread_index != pool->NThreads();
         ++thread_index) {
      predict_buffer_->GetModel(thread_index) = dcomplex();
      if (apply_beam_)
        predict_buffer_->GetPatchModel(thread_index) = dcomplex();

      // When applying beam, simulate into patch vector
      Cube<dcomplex>& simulatedest =
          (apply_beam_ ? predict_buffer_->GetPatchModel(thread_index)
                       : predict_buffer_->GetModel(thread_index));
      simulators.emplace_back(phase_ref_, nSt, baselines_, info().chanFreqs(),
                              info().chanWidths(), station_uvw_, simulatedest,
                              correct_freq_smearing_, stokes_i_only_);
    }
  }
  std::vector<std::shared_ptr<const base::Patch>> curPatches(pool->NThreads());

  if (thread_over_baselines_) {
    aocommon::Barrier barrier(n_threads);
    aocommon::ParallelFor<size_t> loop(n_threads);
    loop.Run(0, n_threads, [&](size_t thread_index) {
      const common::ScopedMicroSecondAccumulator<decltype(predict_time_)>
          scoped_time{predict_time_};
      // OnePredict the source model and apply beam when an entire patch is
      // done
      std::shared_ptr<const base::Patch>& curPatch = curPatches[thread_index];
      casacore::Slicer slicer_dest(
          casacore::IPosition(3, 0, 0, baseline_range[thread_index].first),
          casacore::IPosition(3, actual_nCr, nCh,
                              baseline_range[thread_index].second -
                                  baseline_range[thread_index].first));

      for (size_t source_index = 0; source_index < source_list_.size();
           ++source_index) {
        const bool patchIsFinished =
            curPatch != source_list_[source_index].second &&
            curPatch != nullptr;

        if (apply_beam_ && patchIsFinished) {
          // PatchModel <- SimulBuffer
          Cube<dcomplex>& simulatedest =
              predict_buffer_->GetPatchModel(thread_index);
          simulatedest(slicer_dest) = sim_buffer[thread_index];

          // Apply the beam and add PatchModel to Model
          addBeamToData(curPatch, time, thread_index, nBeamValues,
                        predict_buffer_->GetPatchModel(thread_index).data(),
                        baseline_range[thread_index],
                        station_range[thread_index], barrier, stokes_i_only_);
          // Initialize patchmodel to zero for the next patch
          sim_buffer[thread_index] = dcomplex();
        }
        // Depending on apply_beam_, the following call will add to either
        // the Model or the PatchModel of the predict buffer
        simulators[thread_index].simulate(source_list_[source_index].first);

        curPatch = source_list_[source_index].second;
      }
      // catch last source
      if (apply_beam_ && curPatch != nullptr) {
        // PatchModel <- SimulBuffer
        Cube<dcomplex>& simulatedest =
            predict_buffer_->GetPatchModel(thread_index);
        simulatedest(slicer_dest) = sim_buffer[thread_index];

        addBeamToData(curPatch, time, thread_index, nBeamValues,
                      predict_buffer_->GetPatchModel(thread_index).data(),
                      baseline_range[thread_index], station_range[thread_index],
                      barrier, stokes_i_only_);
      }
      if (!apply_beam_) {
        Cube<dcomplex>& simulatedest = predict_buffer_->GetModel(thread_index);
        simulatedest(slicer_dest) = sim_buffer[thread_index];
      }
    });
  } else {
    pool->For(0, source_list_.size(), [&](size_t source_index, size_t thread) {
      const common::ScopedMicroSecondAccumulator<decltype(predict_time_)>
          scoped_time{predict_time_};
      // OnePredict the source model and apply beam when an entire patch is
      // done
      std::shared_ptr<const base::Patch>& curPatch = curPatches[thread];
      const bool patchIsFinished =
          curPatch != source_list_[source_index].second && curPatch != nullptr;
      if (apply_beam_ && patchIsFinished) {
        // Apply the beam and add PatchModel to Model
        addBeamToData(curPatch, time, thread, nBeamValues,
                      predict_buffer_->GetPatchModel(thread).data(),
                      stokes_i_only_);
        // Initialize patchmodel to zero for the next patch
        predict_buffer_->GetPatchModel(thread) = dcomplex();
      }
      // Depending on apply_beam_, the following call will add to either
      // the Model or the PatchModel of the predict buffer
      simulators[thread].simulate(source_list_[source_index].first);

      curPatch = source_list_[source_index].second;
    });
    // Apply beam to the last patch
    if (apply_beam_) {
      pool->For(0, pool->NThreads(), [&](size_t thread, size_t) {
        const common::ScopedMicroSecondAccumulator<decltype(predict_time_)>
            scoped_time{predict_time_};
        if (curPatches[thread] != nullptr) {
          addBeamToData(curPatches[thread], time, thread, nBeamValues,
                        predict_buffer_->GetPatchModel(thread).data(),
                        stokes_i_only_);
        }
      });
    }
  }

  // Add all thread model data to one buffer
  scratch_buffer.getData() = casacore::Complex();
  casacore::Complex* tdata = scratch_buffer.getData().data();
  const size_t nVisibilities = nBl * nCh * nCr;
  for (size_t thread = 0; thread < std::min(pool->NThreads(), n_threads);
       ++thread) {
    if (stokes_i_only_) {
      for (size_t i = 0, j = 0; i < nVisibilities; i += nCr, j++) {
        tdata[i] += predict_buffer_->GetModel(thread).data()[j];
        tdata[i + nCr - 1] += predict_buffer_->GetModel(thread).data()[j];
      }
    } else {
      std::transform(tdata, tdata + nVisibilities,
                     predict_buffer_->GetModel(thread).data(), tdata,
                     std::plus<dcomplex>());
    }
  }

  // Call ApplyCal step
  if (do_apply_cal_) {
    apply_cal_step_.process(scratch_buffer);
    scratch_buffer = result_step_->get();
    tdata = scratch_buffer.getData().data();
  }

  // Put predict result from temp buffer into the 'real' buffer
  if (operation_ == "replace") {
    buffer_ = scratch_buffer;
  } else {
    buffer_.copy(bufin);
    casacore::Complex* data = buffer_.getData().data();
    if (operation_ == "add") {
      std::transform(data, data + nVisibilities, tdata, data,
                     std::plus<dcomplex>());
    } else if (operation_ == "subtract") {
      std::transform(data, data + nVisibilities, tdata, data,
                     std::minus<dcomplex>());
    }
  }

  timer_.stop();
  getNextStep()->process(buffer_);
  return false;
}

everybeam::vector3r_t OnePredict::dir2Itrf(const MDirection& dir,
                                           MDirection::Convert& measConverter) {
  const MDirection& itrfDir = measConverter(dir);
  const casacore::Vector<double>& itrf = itrfDir.getValue().getValue();
  everybeam::vector3r_t vec;
  vec[0] = itrf[0];
  vec[1] = itrf[1];
  vec[2] = itrf[2];
  return vec;
}

void OnePredict::addBeamToData(std::shared_ptr<const base::Patch> patch,
                               double time, size_t thread, size_t nBeamValues,
                               dcomplex* data0, bool stokesIOnly) {
  // Apply beam for a patch, add result to Model
  MDirection dir(MVDirection(patch->direction().ra, patch->direction().dec),
                 MDirection::J2000);
  everybeam::vector3r_t srcdir = dir2Itrf(dir, meas_convertors_[thread]);

  if (stokesIOnly) {
    const common::ScopedMicroSecondAccumulator<decltype(apply_beam_time_)>
        scoped_time{apply_beam_time_};
    ApplyBeam::applyBeamStokesIArrayFactor(
        info(), time, data0, srcdir, telescope_.get(),
        predict_buffer_->GetScalarBeamValues(thread), false, beam_mode_,
        &mutex_);
  } else {
    const common::ScopedMicroSecondAccumulator<decltype(apply_beam_time_)>
        scoped_time{apply_beam_time_};
    float* dummyweight = nullptr;
    ApplyBeam::applyBeam(info(), time, data0, dummyweight, srcdir,
                         telescope_.get(),
                         predict_buffer_->GetFullBeamValues(thread), false,
                         beam_mode_, false, &mutex_);
  }

  // Add temporary buffer to Model
  std::transform(predict_buffer_->GetModel(thread).data(),
                 predict_buffer_->GetModel(thread).data() + nBeamValues, data0,
                 predict_buffer_->GetModel(thread).data(),
                 std::plus<dcomplex>());
}

void OnePredict::addBeamToData(std::shared_ptr<const base::Patch> patch,
                               double time, size_t thread, size_t nBeamValues,
                               dcomplex* data0,
                               const std::pair<size_t, size_t>& baseline_range,
                               const std::pair<size_t, size_t>& station_range,
                               aocommon::Barrier& barrier, bool stokesIOnly) {
  // Apply beam for a patch, add result to Model
  MDirection dir(MVDirection(patch->direction().ra, patch->direction().dec),
                 MDirection::J2000);
  everybeam::vector3r_t srcdir = dir2Itrf(dir, meas_convertors_[thread]);

  // We use a common buffer to calculate beam values
  const size_t common_thread = 0;
  if (stokesIOnly) {
    const common::ScopedMicroSecondAccumulator<decltype(apply_beam_time_)>
        scoped_time{apply_beam_time_};
    ApplyBeam::applyBeamStokesIArrayFactor(
        info(), time, data0, srcdir, telescope_.get(),
        predict_buffer_->GetScalarBeamValues(common_thread), baseline_range,
        station_range, barrier, false, beam_mode_, &mutex_);
  } else {
    const common::ScopedMicroSecondAccumulator<decltype(apply_beam_time_)>
        scoped_time{apply_beam_time_};
    float* dummyweight = nullptr;
    ApplyBeam::applyBeam(
        info(), time, data0, dummyweight, srcdir, telescope_.get(),
        predict_buffer_->GetFullBeamValues(common_thread), baseline_range,
        station_range, barrier, false, beam_mode_, false, &mutex_);
  }

  // Add temporary buffer to Model
  std::transform(predict_buffer_->GetModel(thread).data(),
                 predict_buffer_->GetModel(thread).data() + nBeamValues, data0,
                 predict_buffer_->GetModel(thread).data(),
                 std::plus<dcomplex>());
}

void OnePredict::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
