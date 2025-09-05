// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Mick Veldhuis

#include "DynSpec.h"

#include <iostream>
#include <string>
#include <complex>

#include <xtensor/xview.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xlayout.hpp>
#include <xtensor/xindex_view.hpp>

#include <aocommon/logger.h>
#include <aocommon/polarization.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include <dp3/base/Direction.h>

#include "../base/FlagCounter.h"

#include "../model/SkyModelCache.h"
#include "../model/SourceDBUtil.h"

using aocommon::Logger;

using dp3::base::Direction;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace {
template <typename T>
std::string ToString(T t) {
  std::stringstream ss;
  ss << std::setprecision(20) << t;
  return ss.str();
}
}  // namespace

namespace dp3 {
namespace steps {

DynSpec::DynSpec(const common::ParameterSet& parset, const std::string& prefix)
    : name_(prefix),
      source_file_name_(parset.getString(prefix + "sourcelist")),
      fits_prefix_(parset.getString(prefix + "fitsprefix", "")),
      model_column_(parset.getString(prefix + "subtractmodelcolumn", "")) {
  // Read directions from list of sources.
  model::SourceDBWrapper source_db =
      model::SkyModelCache::GetInstance()
          .GetSkyModel(source_file_name_)
          .Filter(std::vector<std::string>(),
                  model::SourceDBWrapper::FilterMode::kPattern);
  source_list_ = source_db.MakePatchList();
  if (source_list_.empty()) {
    // TODO: should we add a Patch at the MSs phase centre instead?
    throw std::runtime_error(
        "No sources available, while at least one source needs to be provided "
        "using the 'sourcelist' option...");
  }

  // The dynamic spectra are created by first subtracting the foreground sources
  // and subsequently shifting the phase centre. Either by subtracting a model
  // data column or using H5ParmPredict. The former takes precedent over the
  // latter.
  subtract_model_column_ = !model_column_.empty();
  subtract_with_h5parmpredict_ =
      parset.isDefined(prefix + "h5parmpredict.applycal.parmdb");
  subtract_sources_ = subtract_model_column_ || subtract_with_h5parmpredict_;
  if (subtract_model_column_ && subtract_with_h5parmpredict_) {
    throw std::runtime_error(
        "Either use a model column or H5ParmPredict to subtract sources, "
        "both are provided...");
  }

  if (subtract_sources_) {
    if (subtract_model_column_) {
      model_step_ =
          std::make_shared<MsColumnReader>(parset, prefix, model_column_);
    } else {
      model_step_ =
          std::make_shared<H5ParmPredict>(parset, prefix + "h5parmpredict.");
    }
    model_result_ = std::make_shared<ResultStep>();
    model_step_->setNextStep(model_result_);
  }

  phase_shifts_.reserve(source_list_.size());
  results_.reserve(source_list_.size());
  for (std::shared_ptr<model::Patch>& source : source_list_) {
    const Direction& source_direction = source->Direction();
    std::vector<std::string> phase_center = {ToString(source_direction.ra),
                                             ToString(source_direction.dec)};
    auto phase_shift = std::make_shared<PhaseShift>(
        parset, prefix + source->Name() + ".", phase_center);
    auto result_step = std::make_shared<ResultStep>();
    phase_shift->setNextStep(result_step);

    phase_shifts_.push_back(phase_shift);
    results_.push_back(result_step);
  }
}

void DynSpec::updateInfo(const DPInfo& info_in) {
  Step::updateInfo(info_in);

  // Propagate DPInfo to substeps
  if (subtract_sources_) {
    model_step_->setInfo(info_in);
  }
  for (std::shared_ptr<PhaseShift>& substep : phase_shifts_) {
    substep->setInfo(info_in);
  }

  // Initialise spectra_ with appropriate shape. Current assumption:
  // 4 polarizations are available in the MS.
  if (getInfoOut().ncorr() != 4) {
    throw std::runtime_error(
        "DynSpec assumes that the MS contains 4 instrumental polarizations: "
        "XX, XY, YX, and YY! But only found " +
        std::to_string(getInfoOut().ncorr()) + "...");
  }

  // Set the proper shape of the dynamic spectra, assuming 4 polarization axes.
  const size_t n_directions = source_list_.size();
  constexpr size_t n_correlations = 4;
  const std::array<std::size_t, 4> shape{
      getInfoOut().ntime(), getInfoOut().nchan(), n_correlations, n_directions};
  dynamic_spectra_ = DynamicSpectrumTensor(shape);
}

void DynSpec::show(std::ostream& os) const {
  os << "DynSpec " << name_ << '\n';
  os << "  subtract sources: " << std::boolalpha << subtract_sources_ << '\n';
  os << "  number of target sources: " << dynamic_spectra_.shape(3) << '\n';
}

void DynSpec::showTimings(std::ostream& os, double duration) const {
  const double total_duration = total_timer_.getElapsed();
  os << "  ";
  base::FlagCounter::showPerc1(os, total_timer_.getElapsed(), duration);
  os << " DynSpec " << name_ << '\n';
  os << "          ";
  base::FlagCounter::showPerc1(os, substep_timer_.getElapsed(), total_duration);
  os << " of it spent in substeps" << '\n';
  os << "          ";
  base::FlagCounter::showPerc1(os, computation_timer_.getElapsed(),
                               total_duration);
  os << " of it spent computing spectra" << '\n';
  os << "          ";
  base::FlagCounter::showPerc1(os, write_timer_.getElapsed(), total_duration);
  os << " of it spent writing spectra to disk" << '\n';
}

bool DynSpec::process(std::unique_ptr<DPBuffer> buffer) {
  total_timer_.start();

  // If requested, subtract the foreground before shifting the phase center.
  if (subtract_sources_) {
    substep_timer_.start();

    const common::Fields fields = model_step_->getRequiredFields();
    std::unique_ptr<DPBuffer> input_buffer =
        std::make_unique<DPBuffer>(*buffer, fields);

    model_step_->process(std::move(input_buffer));
    std::unique_ptr<DPBuffer> model_buffer = model_result_->take();

    buffer->GetData() -= model_buffer->GetData();

    substep_timer_.stop();
  }

  DPBuffer::WeightsType weights = buffer->GetWeights();
  xt::filtration(weights, xt::isnan(weights)) = 0.0f;

  // For each source, shift the phase centre and average over all baselines
  // (weighting accordingly), and store the spectra for the current time slot.
  const common::Fields overall_fields =
      dp3::base::GetChainRequiredFields(phase_shifts_[0]);
  for (size_t direction_index = 0; direction_index < phase_shifts_.size();
       ++direction_index) {
    substep_timer_.start();
    std::unique_ptr<DPBuffer> substep_buffer =
        std::make_unique<DPBuffer>(*buffer, overall_fields);
    phase_shifts_[direction_index]->process(std::move(substep_buffer));
    std::unique_ptr<DPBuffer> result_buffer = results_[direction_index]->take();
    substep_timer_.stop();

    computation_timer_.start();
    DPBuffer::DataType weighted_data = weights * result_buffer->GetData();
    xt::filtration(weighted_data, buffer->GetFlags()) =
        std::complex<float>(0.0f, 0.0f);
    xt::filtration(weighted_data, xt::isnan(weighted_data)) =
        std::complex<float>(0.0f, 0.0f);

    xt::xtensor<std::complex<float>, 2> baseline_averaged_data =
        xt::mean(weighted_data, {0}) / xt::sum(weights, {0});

    xt::view(dynamic_spectra_, time_index_, xt::all(), xt::all(),
             direction_index) =
        ComputeAbsoluteStokesParameters(baseline_averaged_data);
    computation_timer_.stop();
  }

  if (time_index_ < getInfoOut().ntime()) ++time_index_;

  total_timer_.stop();
  getNextStep()->process(std::move(buffer));
  return false;
}

xt::xtensor<float, 2> DynSpec::ComputeAbsoluteStokesParameters(
    xt::xtensor<std::complex<float>, 2>& baseline_averaged_data) const {
  using namespace std::complex_literals;

  xt::xtensor<std::complex<float>, 2, xt::layout_type::column_major>
      stokes_parameters(baseline_averaged_data.shape());

  // Convert instrumental polarizations to Stokes parameters
  // I = (XX + YY) / 2
  // Q = (XX - YY) / 2
  // U = (XY + YX) / 2
  // V = 1i * (YX - XY) / 2
  const auto xx = xt::view(baseline_averaged_data, xt::all(), 0);
  const auto xy = xt::view(baseline_averaged_data, xt::all(), 1);
  const auto yx = xt::view(baseline_averaged_data, xt::all(), 2);
  const auto yy = xt::view(baseline_averaged_data, xt::all(), 3);

  xt::view(stokes_parameters, xt::all(), 0) = 0.5 * (xx + yy);
  xt::view(stokes_parameters, xt::all(), 1) = 0.5 * (xx - yy);
  xt::view(stokes_parameters, xt::all(), 2) = 0.5 * (xy + yx);
  xt::view(stokes_parameters, xt::all(), 3) = 0.5if * (yx - xy);

  return xt::abs(stokes_parameters);
}

void DynSpec::WriteSpectraToDisk() {
  write_timer_.start();
  common::DynSpecFitsWriter dynspec_writer;
  dynspec_writer.InitializeTimeAxis(getInfoOut().ntime(),
                                    getInfoOut().timeInterval(),
                                    getInfoOut().firstTime());
  dynspec_writer.InitializeFrequencyAxis(getInfoOut().nchan(),
                                         getInfoOut().chanWidths().at(0),
                                         getInfoOut().chanFreqs().at(0));
  const std::string ms_name = getInfoOut().msName();
  const std::string fits_origin = "DP3/DynSpec";
  const std::string fits_origin_comment =
      "Dynamic spectrum extracted from " + ms_name;
  dynspec_writer.SetOrigin(fits_origin, fits_origin_comment);

  // Store the Observation ID if the MS is formatted in accordance with
  // LOFAR-USG-ICD-005, i.e.: <Prefix><Observation ID>_<Optional
  // Descriptors>_<Filetype>.<Extension>
  if (ms_name.front() == 'L') {
    const size_t delimiter_index = ms_name.find("_");
    const std::string observation_id = ms_name.substr(1, delimiter_index);

    dynspec_writer.SetObsId(observation_id);
  }

  const std::string prefix =
      fits_prefix_.empty() ? fits_prefix_ : fits_prefix_ + "-";
  for (size_t direction_index = 0; direction_index < phase_shifts_.size();
       ++direction_index) {
    std::shared_ptr<model::Patch>& source = source_list_[direction_index];
    const std::string& source_name = source->Name();
    dynspec_writer.SetObjectName(source_name);
    dynspec_writer.SetObjectCoordinates(source->Direction().ra,
                                        source->Direction().dec);

    const std::string filename = prefix + source_name + "-dynspec.fits";
    Logger::Info << "Writing dynamic spectra for " + source_name +
                        " to FITS file " + filename + "..."
                 << '\n';
    dynspec_writer.Write<float>(filename,
                                &dynamic_spectra_(0, 0, 0, direction_index));
  }

  write_timer_.stop();
}

void DynSpec::finish() {
  total_timer_.start();
  WriteSpectraToDisk();
  total_timer_.stop();
  getNextStep()->finish();
}

}  // namespace steps
}  // namespace dp3
