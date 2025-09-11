// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "FastPredict.h"
#include "ApplyBeam.h"

#include <algorithm>
#include <cassert>
#include <iostream>

#include <stdexcept>
#include <xtensor/xview.hpp>
#include <xtensor/xadapt.hpp>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/StreamUtil.h"

#include <dp3/base/DPInfo.h>
#include "../base/FlagCounter.h"
#include "../base/GaussianSource.h"
#include "../base/PointSource.h"
#include "../base/Simulate.h"
#include "../base/Stokes.h"
#include "../base/Telescope.h"

#include "../model/SkyModelCache.h"

#include <predict/BeamResponse.h>
#include <predict/GaussianSourceCollection.h>
#include <predict/PointSourceCollection.h>
#include <predict/PredictPlan.h>
#include <predict/PredictPlanExecCPU.h>
#include <predict/PointSource.h>
#include <predict/Predict.h>
#include <predict/Spectrum.h>

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

#include <predict/PredictPlan.h>

#include <algorithm>
#include <cstddef>
#include <sstream>
#include <string>
#include <utility>

#include <boost/algorithm/string/case_conv.hpp>

using casacore::MDirection;
using casacore::MEpoch;
using casacore::MVEpoch;
using casacore::Quantum;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

FastPredict::FastPredict(const common::ParameterSet& parset,
                         const std::string& prefix,
                         const std::vector<std::string>& source_patterns) {
  if (!source_patterns.empty()) {
    Init(parset, prefix, source_patterns);
  } else {
    const std::vector<std::string> parset_patterns =
        parset.getStringVector(prefix + "sources", std::vector<std::string>());
    Init(parset, prefix, parset_patterns);
  }
}

void FastPredict::Init(const common::ParameterSet& parset,
                       const std::string& prefix,
                       const std::vector<std::string>& sourcePatterns) {
  name_ = prefix;
  source_db_name_ = parset.getString(prefix + "sourcedb");
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

  if (thread_over_baselines_) {
    aocommon::Logger::Warn
        << "Thread over baselines is not supported in FastPredict.\n";
  }

  if (apply_beam_) {
    use_channel_freq_ = parset.getBool(prefix + "usechannelfreq", true);
    one_beam_per_patch_ = parset.getBool(prefix + "onebeamperpatch", false);
    beam_proximity_limit_ =
        parset.getDouble(prefix + "beamproximitylimit", 60.0) *
        (M_PI / (180.0 * 60.0 * 60.0));

    beam_mode_ = everybeam::ParseCorrectionMode(
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
  if (apply_beam_ && beam_mode_ != everybeam::CorrectionMode::kArrayFactor) {
    stokes_i_only_ = false;
  } else {
    stokes_i_only_ = !source_db.CheckPolarized();
  }
  any_orientation_is_absolute_ = source_db.CheckAnyOrientationIsAbsolute();
}

void FastPredict::SetApplyCal(const common::ParameterSet& parset,
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

FastPredict::~FastPredict() = default;

void FastPredict::InitializePlan() {
  const size_t n_stations = getInfoOut().nantenna();

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
  }

  // Create the Measure ITRF conversion info given the array position.
  // The time and direction are filled in later.
  const bool need_meas_converters = moving_phase_ref_ || apply_beam_;
  if (need_meas_converters) {
    // Prepare measures converters
    meas_frame_.set(getInfoOut().arrayPosCopy());
    meas_frame_.set(
        MEpoch(MVEpoch(getInfoOut().startTime() / 86400), MEpoch::UTC));
    meas_converter_.set(MDirection::J2000,
                        MDirection::Ref(MDirection::ITRF, meas_frame_));
  }

  // Initialize the predict plan.
  const size_t n_buffer_correlations =
      stokes_i_only_ ? 1 : getInfoOut().ncorr();
  predict_plan_.nbaselines = getInfoOut().nbaselines();
  predict_plan_.nchannels = getInfoOut().nchan();
  predict_plan_.nstokes = n_buffer_correlations;
  predict_plan_.nstations =
      getInfoOut().nantenna();  // If telescope is homogeneous this value will
                                // be set to one later.
  predict_plan_.compute_stokes_I_only = stokes_i_only_;
  predict_plan_.correct_frequency_smearing = correct_freq_smearing_;
  predict_plan_.apply_beam = apply_beam_;
  predict_plan_.reference = predict::Direction{phase_ref_.ra, phase_ref_.dec};

  const auto& freqs = getInfoOut().chanFreqs();
  predict_plan_.frequencies = xt::adapt(freqs);

  predict_plan_.channel_widths = getInfoOut().chanWidths();
  predict_plan_.baselines = baselines_;
}

void FastPredict::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  if (operation_ == Operation::kReplace)
    GetWritableInfoOut().setBeamCorrectionMode(
        static_cast<int>(everybeam::CorrectionMode::kNone));

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

  InitializePlan();

  if (apply_cal_step_) {
    apply_cal_step_->setInfo(getInfoOut());
    GetWritableInfoOut() = result_step_->getInfoOut();
  }
}

base::Direction FastPredict::GetFirstDirection() const {
  return patch_list_.front()->Direction();
}

void FastPredict::SetOperation(const std::string& operation) {
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

void FastPredict::show(std::ostream& os) const {
  os << "FastPredict " << name_ << '\n';
  os << "  sourcedb:                " << source_db_name_ << '\n';
  os << "   number of patches:      " << patch_list_.size() << '\n';
  os << "   patches clustered:      " << std::boolalpha
     << (!one_beam_per_patch_ && (beam_proximity_limit_ > 0.0)) << '\n';
  os << "   number of components:   " << source_list_.size() << '\n';
  os << "   absolute orientation:   " << std::boolalpha
     << any_orientation_is_absolute_ << '\n';
  os << "   all unpolarized:        " << std::boolalpha << stokes_i_only_
     << '\n';
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

void FastPredict::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " FastPredict " << name_ << '\n';

  /*
   * The timer_ measures the time in a single thread. Both predict_time_ and
   * apply_beam_time_ are the sum of time in multiple threads. This makes it
   * hard to determine the exact time spent in these phases. Instead it shows
   * the percentage spent in these two parts.
   */
  const int64_t time{predict_time_};
  os << "          ";
  base::FlagCounter::showPerc1(os, predict_time_, time);
  os << " of it spent in predict" << '\n';
}

void FastPredict::CopyPredictBufferToData(
    base::DPBuffer::DataType& destination,
    const xt::xtensor<double, 4, xt::layout_type::row_major>& buffer) {
  const size_t nstokes = buffer.shape()[1];
  const size_t nbaselines = buffer.shape()[2];
  const size_t nchannels = buffer.shape()[3];

  const auto real_part = xt::view(buffer, 0, xt::all(), xt::all(), xt::all());
  const auto imag_part = xt::view(buffer, 1, xt::all(), xt::all(), xt::all());

  // destination:  n_baselines, n_channels, buffered_correlations

  if (stokes_i_only_) {
    const size_t ncorr_out = getInfoOut().ncorr();
    for (size_t bl = 0; bl < nbaselines; ++bl) {
      for (size_t ch = 0; ch < nchannels; ++ch) {
        // First correlation (index 0)
        destination(bl, ch, 0) =
            std::complex<float>(static_cast<float>(real_part(0, bl, ch)),
                                static_cast<float>(imag_part(0, bl, ch)));

        // Last correlation (index ncorr_out - 1)
        destination(bl, ch, ncorr_out - 1) = std::complex<float>(
            static_cast<float>(real_part(ncorr_out - 1, bl, ch)),
            static_cast<float>(imag_part(ncorr_out - 1, bl, ch)));
      }
    }
  } else {
    for (size_t stoke = 0; stoke < nstokes; ++stoke) {
      for (size_t bl = 0; bl < nbaselines; ++bl) {
        for (size_t ch = 0; ch < nchannels; ++ch) {
          destination(bl, ch, stoke) =
              std::complex<float>(static_cast<float>(real_part(stoke, bl, ch)),
                                  static_cast<float>(imag_part(stoke, bl, ch)));
        }
      }
    }
  }
}

bool FastPredict::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  // Determine the various sizes.
  const size_t nBl = getInfoOut().nbaselines();
  const size_t nCh = getInfoOut().nchan();
  const size_t nCr = getInfoOut().ncorr();

  base::SplitUvw(uvw_split_index_, baselines_, buffer->GetUvw(), station_uvw_);
  predict_plan_.uvw = station_uvw_;

  double time = buffer->GetTime();

  const bool need_meas_converters = moving_phase_ref_ || apply_beam_;
  if (need_meas_converters) {
    meas_converter_(getInfoOut().delayCenter());
    meas_converter_(getInfoOut().tileBeamDir());
    meas_frame_.resetEpoch(MEpoch(MVEpoch(time / 86400), MEpoch::UTC));
  }

  if (moving_phase_ref_) {
    // Convert phase reference to J2000
    MDirection dirJ2000(
        MDirection::Convert(getInfoOut().phaseCenter(),
                            MDirection::Ref(MDirection::J2000, meas_frame_))());
    Quantum<casacore::Vector<double>> angles = dirJ2000.getAngle();
    phase_ref_ =
        base::Direction(angles.getBaseValue()[0], angles.getBaseValue()[1]);
  }

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

  RunPlan(data, time);

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

void FastPredict::RunPlan(base::DPBuffer::DataType& destination, double time) {
  xt::xtensor<double, 4, xt::layout_type::row_major> global_data(
      {2, predict_plan_.nstokes, predict_plan_.nbaselines,
       predict_plan_.nchannels});

  bool update_beam = false;
  double beam_evaluation_time = time;
  if (apply_beam_) {
    const double time_since_beam_update = std::fabs(time - previous_beam_time_);
    update_beam = time_since_beam_update >= beam_evaluation_interval_;
    if (update_beam) {
      beam_evaluation_time = time + 0.5 * beam_evaluation_interval_;
      previous_beam_time_ = time;
      telescope_->SetTime(beam_evaluation_time);

      if (base::IsHomogeneous(*telescope_)) predict_plan_.nstations = 1;
    }
  }

  const size_t num_threads = aocommon::ThreadPool::GetInstance().NThreads();

  predict_plan_exec_ =
      std::make_unique<predict::PredictPlanExecCPU>(predict_plan_, num_threads);

  assert(predict_plan_exec_);

  const size_t n_buffer_correlations =
      stokes_i_only_ ? 1 : getInfoOut().ncorr();
  const size_t n_baselines = getInfoOut().nbaselines();
  const size_t n_channels = getInfoOut().nchan();

  xt::xtensor<double, 4, xt::layout_type::row_major> model_data_new(
      {2, predict_plan_.nstokes, predict_plan_.nbaselines,
       predict_plan_.nchannels});
  model_data_new.fill(0.0);

  xt::xtensor<double, 4, xt::layout_type::row_major> patch_model_data_new;

  if (predict_plan_.apply_beam) {
    patch_model_data_new.resize(
        {2, n_buffer_correlations, n_baselines, n_channels});
    patch_model_data_new.fill(0.0);
  }

  xt::xtensor<double, 4, xt::layout_type::row_major>& simulator_data_new =
      predict_plan_.apply_beam ? patch_model_data_new : model_data_new;
  predict::PointSourceCollection point_sources;
  predict::GaussianSourceCollection gaussian_sources;

  const model::Patch* patch = nullptr;
  const size_t num_sources = source_list_.size();

  for (size_t source_index = 0; source_index != num_sources; ++source_index) {
    const std::shared_ptr<base::ModelComponent>& source_ptr =
        source_list_[source_index].first;
    patch = source_list_[source_index].second.get();
    assert(patch);

    if (const base::PointSource* point_source =
            dynamic_cast<base::PointSource*>(source_ptr.get())) {
      const base::Stokes dp3_stokes = point_source->stokes();
      const predict::Stokes fast_predict_stokes(dp3_stokes.I, dp3_stokes.Q,
                                                dp3_stokes.U, dp3_stokes.V);
      predict::Spectrum spectrum{fast_predict_stokes,
                                 point_source->referenceFreq(),
                                 point_source->polarizationAngle(),
                                 point_source->polarizedFraction(),
                                 point_source->rotationMeasure(),
                                 point_source->hasRotationMeasure(),
                                 point_source->hasLogarithmicSI()};

      spectrum.SetSpectralTerms(point_source->referenceFreq(),
                                point_source->hasLogarithmicSI(),
                                point_source->spectrum());

      const size_t patch_index = patch->Index();

      if (const base::GaussianSource* gaussian_source =
              dynamic_cast<base::GaussianSource*>(source_ptr.get());
          gaussian_source) {
        // It's specifically a GaussianSource
        gaussian_sources.Add(predict::GaussianSource{
            predict::Direction{point_source->direction().ra,
                               point_source->direction().dec},
            spectrum, gaussian_source->getPositionAngle(),
            gaussian_source->getPositionAngleIsAbsolute(),
            gaussian_source->getMinorAxis(), gaussian_source->getMajorAxis(),
            patch_index});
        if (predict_plan_.apply_beam) {
          gaussian_sources.AddBeamDirection(
              patch_index, predict::Direction{patch->Direction().ra,
                                              patch->Direction().dec});
        }
      } else {
        point_sources.Add(predict::PointSource{
            predict::Direction{point_source->direction().ra,
                               point_source->direction().dec},
            spectrum, patch_index});
        if (predict_plan_.apply_beam) {
          point_sources.AddBeamDirection(
              patch_index, predict::Direction{patch->Direction().ra,
                                              patch->Direction().dec});
        }
      }
    }
  }

  if (predict_plan_.apply_beam) {
    point_sources.UpdateBeams();
    gaussian_sources.UpdateBeams();
  }

  {
    const common::ScopedMicroSecondAccumulator scoped_time(predict_time_);
    if (predict_plan_.apply_beam) {
      constexpr size_t field_id = 0;
      predict::BeamResponsePlan beam_response_plan{telescope_.get(), time,
                                                   field_id, beam_mode_, false};

      beam_response_plan.SetFrequencies(predict_plan_.frequencies);
      beam_response_plan.SetBaselines(predict_plan_.baselines);

      predict_.runWithStrategy(*predict_plan_exec_, beam_response_plan,
                               point_sources, gaussian_sources,
                               simulator_data_new, meas_converter_,
                               predict::computation_strategy::XSIMD);
    } else {
      predict_.runWithStrategy(*predict_plan_exec_, point_sources,
                               simulator_data_new,
                               predict::computation_strategy::XSIMD);
      predict_.runWithStrategy(*predict_plan_exec_, gaussian_sources,
                               simulator_data_new,
                               predict::computation_strategy::XSIMD);
    }
  }

  global_data = simulator_data_new;

  CopyPredictBufferToData(destination, global_data);
}

void FastPredict::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
