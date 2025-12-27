// OnePredict.cc: DP3 step class that predicts visibilities.
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "OnePredict.h"
#include "ApplyBeam.h"

#include <algorithm>
#include <cassert>
#include <iostream>

#include <xtensor/xview.hpp>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/StreamUtil.h"

#include "../parmdb/ParmDBMeta.h"
#include "../parmdb/PatchInfo.h"
#include "../parmdb/SkymodelToSourceDB.h"

#include "base/DPInfo.h"
#include "../base/FlagCounter.h"
#include "../base/GaussianSource.h"
#include "../base/PointSource.h"
#include "../base/Simulate.h"
#include "../base/Simulator.h"
#include "../base/Stokes.h"
#include "../base/Telescope.h"

#include "../model/SkyModelCache.h"

#include "../parmdb/SourceDB.h"

#include <aocommon/barrier.h>
#include <aocommon/logger.h>
#include <aocommon/recursivefor.h>
#include <aocommon/staticfor.h>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/OS/File.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/tables/Tables/RefRows.h>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <optional>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/algorithm/string/case_conv.hpp>

using casacore::MDirection;
using casacore::MEpoch;
using casacore::MVDirection;
using casacore::MVEpoch;
using casacore::Quantum;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::PredictModel;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

OnePredict::OnePredict(const common::ParameterSet& parset,
                       const std::string& prefix,
                       const std::vector<std::string>& source_patterns)
    : measures_mutex_(nullptr) {
  if (!source_patterns.empty()) {
    init(parset, prefix, source_patterns);
  } else {
    const std::vector<std::string> parset_patterns =
        parset.getStringVector(prefix + "sources", std::vector<std::string>());
    init(parset, prefix, parset_patterns);
  }
}

void OnePredict::init(const common::ParameterSet& parset,
                      const std::string& prefix,
                      const std::vector<std::string>& sourcePatterns) {
  name_ = prefix;
  source_db_name_ = parset.getString(prefix + "sourcedb");
  correct_time_smearing_ =
      parset.getBool(prefix + "correcttimesmearing", false);
  correct_freq_smearing_ =
      parset.getBool(prefix + "correctfreqsmearing", false);

  SetOperation(parset.getString(prefix + "operation", "replace"));
  output_data_name_ = parset.getString(prefix + "outputmodelname", "");

  apply_beam_ = parset.getBool(prefix + "usebeammodel", false);
  coefficients_path_ = parset.getString(prefix + "coefficients_path", "");
  beam_evaluation_interval_ = parset.getDouble(prefix + "beam_interval", 0.0);
  thread_over_baselines_ = parset.getBool(prefix + "parallelbaselines", false);
  debug_level_ = parset.getInt(prefix + "debuglevel", 0);
  patch_list_.clear();

  // Save directions specifications to pass to applycal
  std::stringstream ss;
  ss << sourcePatterns;
  direction_str_ = ss.str();

  aocommon::Logger::Debug << "Loading " << source_db_name_
                          << " in predict step for direction " << direction_str_
                          << ".\n";
  model::SourceDBWrapper source_db =
      model::SkyModelCache::GetInstance().GetSkyModel(source_db_name_);
  source_db.Filter(sourcePatterns,
                   model::SourceDBWrapper::FilterMode::kPattern);
  try {
    patch_list_ = source_db.MakePatchList();
    if (patch_list_.empty()) {
      throw std::runtime_error("Couldn't find patch for direction " +
                               direction_str_);
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

    beam_mode_ = everybeam::ParseBeamMode(
        parset.getString(prefix + "beammode", "default"));

    std::string element_model =
        parset.getString(prefix + "elementmodel", "default");
    element_response_model_ =
        everybeam::ElementResponseModelFromString(element_model);

    // By default, a source model has each direction in one patch. Therefore,
    // if one-beam-per-patch is requested, we don't have to do anything.
    if (!one_beam_per_patch_) {
      if (beam_proximity_limit_ > 0.0) {
        // Rework patch list to cluster proximate sources
        aocommon::Logger::Debug << "Clustering proximate sources for direction "
                                << direction_str_ << ".\n";
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
    SetApplyCal(parset, prefix + "applycal.");
  }

  source_list_ = makeSourceList(patch_list_);
  SetPatchIndices(patch_list_);

  // Determine whether any sources are polarized. If not, enable
  // Stokes-I-only mode (note that this mode cannot be used with apply_beam_)
  if (apply_beam_ && beam_mode_ != everybeam::BeamMode::kArrayFactor) {
    stokes_i_only_ = false;
  } else {
    stokes_i_only_ = !source_db.CheckPolarized();
  }
  any_orientation_is_absolute_ = source_db.CheckAnyOrientationIsAbsolute();
}

void OnePredict::SetApplyCal(const common::ParameterSet& parset,
                             const std::string& prefix) {
  apply_cal_step_ =
      std::make_shared<ApplyCal>(parset, prefix, true, direction_str_);
  if (operation_ != Operation::kReplace &&
      parset.getBool(prefix + "applycal.updateweights", false))
    throw std::invalid_argument(
        "Weights cannot be updated when operation is not replace");
  result_step_ = std::make_shared<ResultStep>();
  apply_cal_step_->setNextStep(result_step_);
}

OnePredict::~OnePredict() = default;

void OnePredict::initializeThreadData() {
  const size_t n_stations = getInfoOut().nantenna();
  const size_t nThreads = aocommon::ThreadPool::GetInstance().NThreads();

  station_uvw_.resize({n_stations, 3});

  std::vector<std::array<double, 3>> antenna_pos(
      getInfoOut().antennaPos().size());
  for (unsigned int i = 0; i < getInfoOut().antennaPos().size(); ++i) {
    casacore::Quantum<casacore::Vector<double>> pos =
        getInfoOut().antennaPos()[i].get("m");
    antenna_pos[i][0] = pos.getValue()[0];
    antenna_pos[i][1] = pos.getValue()[1];
    antenna_pos[i][2] = pos.getValue()[2];
  }

  uvw_split_index_ =
      base::SetupUvwSplitting(getInfoOut().nantenna(), getInfoOut().getAnt1(),
                              getInfoOut().getAnt2(), antenna_pos);

  if (apply_beam_) {
    telescope_ =
        base::GetTelescope(getInfoOut().msName(), element_response_model_,
                           use_channel_freq_, coefficients_path_);
    predict_buffers_ =
        std::make_shared<std::vector<base::PredictBuffer>>(nThreads);
    // TODO:
    // It would be better to use EveryBeam's function to get all antenna values
    // at once, as this would enable this optimization without work here.
    const bool is_dish_telescope = base::IsHomogeneous(*telescope_);
    const size_t n_channels = getInfoOut().nchan();
    predict_buffers_->front().Resize(
        1, n_channels, (is_dish_telescope ? 1 : n_stations), !stokes_i_only_);
  }

  // Create the Measure ITRF conversion info given the array position.
  // The time and direction are filled in later.
  meas_convertors_.resize(nThreads);
  meas_frame_.resize(nThreads);

  for (size_t thread = 0; thread < nThreads; ++thread) {
    const bool need_meas_converters = moving_phase_ref_ || apply_beam_;
    if (need_meas_converters) {
      // Prepare measures converters
      meas_frame_[thread].set(getInfoOut().arrayPosCopy());
      meas_frame_[thread].set(
          MEpoch(MVEpoch(getInfoOut().startTime() / 86400), MEpoch::UTC));
      meas_convertors_[thread].set(
          MDirection::J2000,
          MDirection::Ref(MDirection::ITRF, meas_frame_[thread]));
    }
  }
}

void OnePredict::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  if (operation_ == Operation::kReplace)
    GetWritableInfoOut().setBeamCorrectionMode(
        static_cast<int>(everybeam::BeamMode::kNone));

  for (size_t bl = 0; bl != getInfoOut().nbaselines(); ++bl) {
    baselines_.emplace_back(getInfoOut().getAnt1()[bl],
                            getInfoOut().getAnt2()[bl]);
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

  /**
   * Length of a sidereal day in seconds.
   */
  constexpr double kSiderealDay = 86164.0905;

  const double angular_speed = 2 * M_PI * infoIn.timeInterval() / kSiderealDay;

  // Time smearing is computed by rotating the u,v,w coordinates around the
  // earth axis. The earth axis is assumed here to point to the position
  // dec=90deg in J2000 coordinates. This ignores precession since the true
  // position is dec=90deg for the current epoch. A more accurate method, but
  // probably unnecessary, would be to compute the actual NCP for the current
  // epoch

  scaled_ncp_uvw_.resize(3);
  // In base/Simulator.cc this element is assumed to be zero. If it is updated
  // here to a more accurate value, please update the usage in Simulator as well
  scaled_ncp_uvw_[0] = 0.0;
  scaled_ncp_uvw_[1] =
      angular_speed * std::cos(infoIn.phaseCenterDirection().dec);
  scaled_ncp_uvw_[2] =
      angular_speed * std::sin(infoIn.phaseCenterDirection().dec);

  initializeThreadData();

  if (apply_cal_step_) {
    apply_cal_step_->setInfo(getInfoOut());
    GetWritableInfoOut() = result_step_->getInfoOut();
  }
}

base::Direction OnePredict::GetFirstDirection() const {
  return patch_list_.front()->Direction();
}

void OnePredict::SetOperation(const std::string& operation) {
  if (operation == "replace") {
    operation_ = Operation::kReplace;
  } else if (operation == "add") {
    operation_ = Operation::kAdd;
  } else if (operation == "subtract") {
    operation_ = Operation::kSubtract;
  } else {
    throw std::invalid_argument(
        "Operation must be 'replace', 'add' or 'subtract'.");
  }
}

void OnePredict::show(std::ostream& os) const {
  os << "OnePredict " << name_ << '\n';
  os << "  sourcedb:                " << source_db_name_ << '\n';
  os << "   number of patches:      " << patch_list_.size() << '\n';
  os << "   patches clustered:      " << std::boolalpha
     << (!one_beam_per_patch_ && (beam_proximity_limit_ > 0.0)) << '\n';
  os << "   number of components:   " << source_list_.size() << '\n';
  os << "   absolute orientation:   " << std::boolalpha
     << any_orientation_is_absolute_ << '\n';
  os << "   all unpolarized:        " << std::boolalpha << stokes_i_only_
     << '\n';
  os << "   correct time smearing:  " << std::boolalpha
     << correct_time_smearing_ << '\n';
  os << "   correct freq smearing:  " << std::boolalpha
     << correct_freq_smearing_ << '\n';
  os << "  apply beam:              " << std::boolalpha << apply_beam_ << '\n';
  if (apply_beam_) {
    os << "   mode:                   " << everybeam::ToString(beam_mode_);
    os << '\n';
    os << "   use channelfreq:        " << std::boolalpha << use_channel_freq_
       << '\n';
    os << "   one beam per patch:     " << std::boolalpha << one_beam_per_patch_
       << '\n';
    os << "   beam proximity limit:   "
       << (beam_proximity_limit_ * (180.0 * 60.0 * 60.0) / M_PI) << " arcsec\n";
    os << "   beam interval:          " << beam_evaluation_interval_ << '\n';
  }
  if (!output_data_name_.empty()) {
    os << "  outputmodelname:         " << output_data_name_ << '\n';
  }
  os << "  operation:               ";
  switch (operation_) {
    case Operation::kReplace:
      os << "replace\n";
      break;
    case Operation::kAdd:
      os << "add\n";
      break;
    case Operation::kSubtract:
      os << "subtract\n";
      break;
  }
  if (apply_cal_step_) {
    apply_cal_step_->show(os);
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

void OnePredict::CopyPredictBufferToData(
    base::DPBuffer::DataType& destination,
    const aocommon::xt::UTensor<std::complex<double>, 3>& buffer) {
  if (stokes_i_only_) {
    // Add the predicted model to the first and last correlation.
    const size_t n_correlations = getInfoOut().ncorr();
    auto dest_view = xt::view(destination, xt::all(), xt::all(),
                              xt::keep(0, n_correlations - 1));
    // Without explicit casts, XTensor does not know what to do.
    dest_view = xt::cast<std::complex<float>>(buffer);
  } else {
    // Without explicit casts, XTensor does not know what to do.
    destination = xt::cast<std::complex<float>>(buffer);
  }
}

bool OnePredict::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  // Determine the various sizes.
  const size_t nSt = getInfoOut().nantenna();
  const size_t nBl = getInfoOut().nbaselines();
  const size_t nCh = getInfoOut().nchan();
  const size_t nCr = getInfoOut().ncorr();
  const size_t nThreads = aocommon::ThreadPool::GetInstance().NThreads();

  base::SplitUvw(uvw_split_index_, baselines_, buffer->GetUvw(), station_uvw_);

  double time = buffer->GetTime();

  size_t n_threads = aocommon::ThreadPool::GetInstance().NThreads();
  const bool need_meas_converters = moving_phase_ref_ || apply_beam_;
  if (need_meas_converters) {
    // Because multiple predict steps might be predicting simultaneously, and
    // Casacore is not thread safe, this needs synchronization.
    std::unique_lock<std::mutex> lock;
    if (measures_mutex_ != nullptr)
      lock = std::unique_lock<std::mutex>(*measures_mutex_);
    for (size_t thread = 0; thread != n_threads; ++thread) {
      meas_frame_[thread].resetEpoch(
          MEpoch(MVEpoch(time / 86400), MEpoch::UTC));
      // The convertors are only thread safe after they have been used once;
      // therefore do these dummy conversions. It's unknown if a single
      // conversion using a plain casacore::MDirection() also works, so keep
      // these two conversions for now.
      meas_convertors_[thread](getInfoOut().delayCenter());
      meas_convertors_[thread](getInfoOut().tileBeamDir());
    }
  }

  if (moving_phase_ref_) {
    // Convert phase reference to J2000
    MDirection dirJ2000(MDirection::Convert(
        getInfoOut().phaseCenter(),
        MDirection::Ref(MDirection::J2000, meas_frame_[0]))());
    Quantum<casacore::Vector<double>> angles = dirJ2000.getAngle();
    phase_ref_ =
        base::Direction(angles.getBaseValue()[0], angles.getBaseValue()[1]);
  }

  std::vector<base::Simulator> simulators;
  simulators.reserve(n_threads);

  std::vector<std::pair<size_t, size_t>> baseline_range;
  std::vector<xt::xtensor<std::complex<double>, 3>> sim_buffer;
  std::vector<std::vector<std::pair<size_t, size_t>>> baselines_split;
  std::vector<std::pair<size_t, size_t>> station_range;

  // Take ownership of the input visibilities if we need them later.
  if (operation_ == Operation::kAdd || operation_ == Operation::kSubtract ||
      !output_data_name_.empty()) {
    input_data_ = buffer->TakeData();
  }

  // Determine destination of the predicted visibilities
  if (!output_data_name_.empty()) {
    buffer->AddData(output_data_name_);
  }
  DPBuffer::DataType& data = buffer->GetData(output_data_name_);
  data.resize({nBl, nCh, nCr});
  data.fill(std::complex<float>(0.0, 0.0));

  const size_t actual_nCr = (stokes_i_only_ ? 1 : nCr);
  if (thread_over_baselines_) {
    std::unique_ptr<PredictModel> model_buffer = std::make_unique<PredictModel>(
        nThreads, stokes_i_only_ ? 1 : getInfoOut().ncorr(), nCh, nBl,
        apply_beam_);

    // Reduce the number of threads if there are not enough baselines.
    n_threads = std::min(n_threads, nBl);

    // All threads process 'baselines_per_thread' baselines.
    // The first 'remaining_baselines' threads process an extra baseline.
    const size_t baselines_per_thread = nBl / n_threads;
    const size_t remaining_baselines = nBl % n_threads;

    baseline_range.resize(n_threads);
    sim_buffer.resize(n_threads);
    baselines_split.resize(n_threads);
    if (apply_beam_) {
      station_range.resize(n_threads);
    }

    // Index of the first baseline for the current thread. The loop below
    // updates this variable in each iteration.
    size_t first_baseline = 0;
    for (size_t thread_index = 0; thread_index != n_threads; ++thread_index) {
      const size_t chunk_size =
          baselines_per_thread + ((thread_index < remaining_baselines) ? 1 : 0);

      baseline_range[thread_index] =
          std::make_pair(first_baseline, first_baseline + chunk_size);
      sim_buffer[thread_index].resize({chunk_size, nCh, actual_nCr});

      baselines_split[thread_index].resize(chunk_size);
      std::copy_n(std::next(baselines_.cbegin(), first_baseline), chunk_size,
                  baselines_split[thread_index].begin());

      first_baseline += chunk_size;  // Update for the next loop iteration.
    }
    // Verify that all baselines are assigned to threads.
    assert(first_baseline == nBl);

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

    aocommon::RecursiveFor::NestedRun(0, n_threads, [&](size_t thread_index) {
      const std::complex<double> zero(0.0, 0.0);
      model_buffer->GetModel(thread_index).fill(zero);
      if (apply_beam_) model_buffer->GetPatchModel(thread_index).fill(zero);
      sim_buffer[thread_index].fill(zero);
    });

    // Keep this loop single threaded, I'm not sure if Simulator constructor
    // is thread safe.
    for (size_t thread_index = 0; thread_index != n_threads; ++thread_index) {
      // When applying beam, simulate into patch vector
      // Create a Casacore view since the Simulator still uses Casacore.
      xt::xtensor<std::complex<double>, 3>& thread_buffer =
          sim_buffer[thread_index];
      const casacore::IPosition shape(3, thread_buffer.shape(2),
                                      thread_buffer.shape(1),
                                      thread_buffer.shape(0));
      casacore::Cube<std::complex<double>> simulatedest(
          shape, thread_buffer.data(), casacore::SHARE);

      simulators.emplace_back(
          phase_ref_, nSt, baselines_split[thread_index],
          getInfoOut().chanFreqs(), getInfoOut().chanWidths(), scaled_ncp_uvw_,
          station_uvw_, simulatedest, correct_time_smearing_,
          correct_freq_smearing_, stokes_i_only_);
    }

    std::vector<std::shared_ptr<const model::Patch>> curPatches(n_threads);

    aocommon::Barrier barrier(n_threads);
    // We need to create local threads here because we need to
    // sync only those using the barrier
    aocommon::RecursiveFor::NestedRun(0, n_threads, [&](size_t thread_index) {
      const common::ScopedMicroSecondAccumulator scoped_time(predict_time_);
      // Predict the source model and apply beam when an entire patch is
      // done
      const model::Patch* curPatch = curPatches[thread_index].get();

      for (size_t source_index = 0; source_index < source_list_.size();
           ++source_index) {
        const bool patchIsFinished =
            curPatch != source_list_[source_index].second.get() &&
            curPatch != nullptr;

        if (apply_beam_ && patchIsFinished) {
          // PatchModel <- SimulBuffer
          aocommon::xt::UTensor<std::complex<double>, 3>& patch_model =
              model_buffer->GetPatchModel(thread_index);
          xt::view(patch_model,
                   xt::range(baseline_range[thread_index].first,
                             baseline_range[thread_index].second),
                   xt::all(), xt::all()) = sim_buffer[thread_index];

          // Apply the beam and add PatchModel to Model
          addBeamToDataRange(
              *curPatch, model_buffer->GetModel(thread_index), time,
              thread_index, patch_model, baseline_range[thread_index],
              station_range[thread_index], barrier, stokes_i_only_);
          // Initialize patchmodel to zero for the next patch
          sim_buffer[thread_index].fill(std::complex<double>(0.0, 0.0));
        }
        // Depending on apply_beam_, the following call will add to either
        // the Model or the PatchModel of the predict buffer
        simulators[thread_index].simulate(source_list_[source_index].first);

        curPatch = source_list_[source_index].second.get();
      }
      // catch last source
      if (apply_beam_ && curPatch != nullptr) {
        // PatchModel <- SimulBuffer
        aocommon::xt::UTensor<std::complex<double>, 3>& patch_model =
            model_buffer->GetPatchModel(thread_index);
        xt::view(patch_model,
                 xt::range(baseline_range[thread_index].first,
                           baseline_range[thread_index].second),
                 xt::all(), xt::all()) = sim_buffer[thread_index];

        addBeamToDataRange(
            *curPatch, model_buffer->GetModel(thread_index), time, thread_index,
            patch_model, baseline_range[thread_index],
            station_range[thread_index], barrier, stokes_i_only_);
      }
      if (!apply_beam_) {
        aocommon::xt::UTensor<std::complex<double>, 3>& model =
            model_buffer->GetModel(thread_index);
        xt::view(model,
                 xt::range(baseline_range[thread_index].first,
                           baseline_range[thread_index].second),
                 xt::all(), xt::all()) = sim_buffer[thread_index];
      }
    });

    // Add all thread model data to one buffer
    for (size_t thread = 1; thread < n_threads; ++thread) {
      // Sum thread model data in their own container (doubles) to prevent
      // rounding errors when writing to the data member of the DPBuffer
      // (floats).
      model_buffer->GetModel(0) += model_buffer->GetModel(thread);
    }

    CopyPredictBufferToData(data, model_buffer->GetModel(0));
  } else {
    PredictWithSourceParallelization(data, time);
  }

  if (apply_cal_step_) {
    apply_cal_step_->process(std::move(buffer));
    buffer = result_step_->take();
  }

  if (operation_ == Operation::kAdd) {
    data += input_data_;
  } else if (operation_ == Operation::kSubtract) {
    data = input_data_ - data;
  }
  if (!output_data_name_.empty()) {
    // Put the input visibilities back to the main buffer when needed.
    buffer->GetData() = std::move(input_data_);
  }

  timer_.stop();

  getNextStep()->process(std::move(buffer));
  return false;
}

void OnePredict::PredictSourceRange(
    aocommon::xt::UTensor<std::complex<double>, 3>& result, size_t start,
    size_t end, size_t thread_index, std::mutex& mutex, double time,
    bool update_beam) {
  const size_t n_stations = getInfoOut().nantenna();
  const size_t n_baselines = getInfoOut().nbaselines();
  const size_t n_channels = getInfoOut().nchan();
  const size_t n_buffer_correlations =
      stokes_i_only_ ? 1 : getInfoOut().ncorr();

  aocommon::xt::UTensor<std::complex<double>, 3> model_data(
      {n_baselines, n_channels, n_buffer_correlations},
      std::complex<double>(0.0, 0.0));
  aocommon::xt::UTensor<std::complex<double>, 3> patch_model_data;

  const size_t start_patch = source_list_[start].second->Index();
  if (apply_beam_) {
    base::PredictBuffer& buffer = (*predict_buffers_)[thread_index];
    if (update_beam) {
      const bool is_homogeneous = base::IsHomogeneous(*telescope_);
      const size_t n_patches =
          source_list_.empty()
              ? 0
              : source_list_[end - 1].second->Index() + 1 - start_patch;
      buffer.Resize(n_patches, n_channels, (is_homogeneous ? 1 : n_stations),
                    !stokes_i_only_);
    }
    patch_model_data.resize({n_baselines, n_channels, n_buffer_correlations});
    patch_model_data.fill(std::complex<double>(0.0, 0.0));
  }

  // Create a Casacore view since the Simulator still uses Casacore.
  // model_data.shape() and patch_model_data_data always have the same shape.
  const casacore::IPosition shape(3, model_data.shape(2), model_data.shape(1),
                                  model_data.shape(0));
  std::complex<double>* simulator_data =
      apply_beam_ ? patch_model_data.data() : model_data.data();
  casacore::Cube<std::complex<double>> casacore_data(shape, simulator_data,
                                                     casacore::SHARE);
  base::Simulator simulator(
      phase_ref_, n_stations, baselines_, getInfoOut().chanFreqs(),
      getInfoOut().chanWidths(), scaled_ncp_uvw_, station_uvw_, casacore_data,
      correct_time_smearing_, correct_freq_smearing_, stokes_i_only_);

  const common::ScopedMicroSecondAccumulator scoped_time(predict_time_);
  const model::Patch* patch = nullptr;

  for (size_t source_index = start; source_index != end; ++source_index) {
    // Predict the source model and apply beam when an entire patch is
    // done
    const bool patch_is_finished =
        patch != source_list_[source_index].second.get() && patch != nullptr;
    if (apply_beam_ && patch_is_finished) {
      // Apply the beam and add PatchModel to Model
      addBeamToData(*patch, patch->Index() - start_patch, model_data, time,
                    update_beam, thread_index, patch_model_data,
                    stokes_i_only_);
      // Initialize patchmodel to zero for the next patch
      patch_model_data.fill(std::complex<double>(0.0, 0.0));
    }
    // Depending on apply_beam_, the following call will add to either
    // the Model or the PatchModel predict buffer
    simulator.simulate(source_list_[source_index].first);

    patch = source_list_[source_index].second.get();
  }

  if (apply_beam_ && patch != nullptr) {
    // Apply beam to the last patch
    const common::ScopedMicroSecondAccumulator scoped_time(predict_time_);
    addBeamToData(*patch, patch->Index() - start_patch, model_data, time,
                  update_beam, thread_index, patch_model_data, stokes_i_only_);
  }

  // Add this thread's data to the global buffer
  std::lock_guard lock(mutex);
  result += model_data;
}

void OnePredict::PredictWithSourceParallelization(
    base::DPBuffer::DataType& destination, double time) {
  const size_t n_baselines = getInfoOut().nbaselines();
  const size_t n_channels = getInfoOut().nchan();
  const size_t buffered_correlations =
      stokes_i_only_ ? 1 : getInfoOut().ncorr();

  aocommon::xt::UTensor<std::complex<double>, 3> global_data(
      {n_baselines, n_channels, buffered_correlations},
      std::complex<double>(0.0, 0.0));

  bool update_beam = false;
  double beam_evaluation_time = time;
  if (apply_beam_) {
    const double time_since_beam_update = std::fabs(time - previous_beam_time_);
    update_beam = time_since_beam_update >= beam_evaluation_interval_;
    if (update_beam) {
      beam_evaluation_time = time + 0.5 * beam_evaluation_interval_;
      previous_beam_time_ = time;
      telescope_->SetTime(beam_evaluation_time);
    }
  }

  std::mutex mutex;
  aocommon::StaticFor<size_t> loop;
  // The way source are split into consecutive subranges
  // is important: it makes sure that a single patch is mostly calculated
  // by a single thread, which limits duplicate beam evaluations.
  // The work per thread needs to be the same for each timestep, because
  // threads store beam values for the patches they use.
  loop.Run(0, source_list_.size(),
           [&](size_t start, size_t end, size_t thread_index) {
             PredictSourceRange(global_data, start, end, thread_index, mutex,
                                beam_evaluation_time, update_beam);
           });

  CopyPredictBufferToData(destination, global_data);
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

void OnePredict::addBeamToData(
    const model::Patch& patch, size_t buffer_index,
    aocommon::xt::UTensor<std::complex<double>, 3>& model_data, double time,
    bool update_beam, size_t thread,
    aocommon::xt::UTensor<std::complex<double>, 3>& data, bool stokesIOnly) {
  // Apply beam for a patch, add result to Model
  base::PredictBuffer& buffer = (*predict_buffers_)[thread];

  const common::ScopedMicroSecondAccumulator scoped_time(apply_beam_time_);
  if (stokesIOnly) {
    std::complex<double>* values = buffer.GetScalarBeamValues(buffer_index);
    if (update_beam) {
      const MDirection dir(
          MVDirection(patch.Direction().ra, patch.Direction().dec),
          MDirection::J2000);
      const everybeam::vector3r_t srcdir =
          dir2Itrf(dir, meas_convertors_[thread]);
      ComputeArrayFactor(getInfoOut(), time, srcdir, telescope_.get(), values,
                         false, &mutex_, {});
    }
    ApplyArrayFactorAndAdd(getInfoOut(), buffer.NStations(), data, model_data,
                           values);
  } else {
    aocommon::MC2x2* values = buffer.GetFullBeamValues(buffer_index);
    if (update_beam) {
      const MDirection dir(
          MVDirection(patch.Direction().ra, patch.Direction().dec),
          MDirection::J2000);
      const everybeam::vector3r_t srcdir =
          dir2Itrf(dir, meas_convertors_[thread]);
      ComputeBeam(getInfoOut(), time, srcdir, telescope_.get(), values, false,
                  beam_mode_, &mutex_, {});
    }
    ApplyBeamToDataAndAdd(getInfoOut(), buffer.NStations(), data, model_data,
                          values);
  }
}

void OnePredict::addBeamToDataRange(
    const model::Patch& patch,
    aocommon::xt::UTensor<std::complex<double>, 3>& model_data, double time,
    size_t thread, aocommon::xt::UTensor<std::complex<double>, 3>& data,
    const std::pair<size_t, size_t>& baseline_range,
    const std::pair<size_t, size_t>& station_range, aocommon::Barrier& barrier,
    bool stokesIOnly) {
  // Apply beam for a patch, add result to Model
  MDirection dir(MVDirection(patch.Direction().ra, patch.Direction().dec),
                 MDirection::J2000);
  const everybeam::vector3r_t srcdir = dir2Itrf(dir, meas_convertors_[thread]);

  // We use a common buffer to calculate beam values
  base::PredictBuffer& buffer = predict_buffers_->front();
  if (stokesIOnly) {
    const common::ScopedMicroSecondAccumulator scoped_time(apply_beam_time_);
    ApplyBeam::ApplyBaselineBasedArrayFactor(
        getInfoOut(), time, data.data(), srcdir, telescope_.get(),
        buffer.GetScalarBeamValues(0), baseline_range, station_range, barrier,
        false, beam_mode_, &mutex_);
  } else {
    const common::ScopedMicroSecondAccumulator scoped_time(apply_beam_time_);
    float* weights = nullptr;
    ApplyBeam::ApplyBaselineBasedBeam(
        getInfoOut(), time, data.data(), weights, srcdir, telescope_.get(),
        buffer.GetFullBeamValues(0), baseline_range, station_range, barrier,
        false, beam_mode_, false, &mutex_);
  }

  // Add temporary buffer to Model
  model_data += data;
}

void OnePredict::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
