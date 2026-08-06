// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "AartfaacReader.h"

#include <filesystem>
#include <iostream>
#include <memory>
#include <utility>

#include <aocommon/constants.h>

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MCEpoch.h>
#include <casacore/measures/Measures/Muvw.h>
#include <casacore/measures/Measures/MPosition.h>

#include <boost/filesystem/operations.hpp>

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/misc/xcomplex.hpp>
#include <xtensor/views/xindex_view.hpp>
#include <xtensor/views/xview.hpp>

#include "aartfaac/File.h"
#include "aartfaac/AntennaConfig.h"
#include "base/SubtableWriter.h"
#include "common/baseline_indices.h"

using casacore::Double;
using casacore::MDirection;
using casacore::MPosition;
using casacore::Muvw;
using casacore::Time;
using dp3::aartfaac::AntennaConfig;
using dp3::aartfaac::File;
using dp3::aartfaac::Timestep;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace {
const double DAY_IN_SECONDS = 3600. * 24.;
const Time UNIX_EPOCH = Time(1970, 1, 1, 0, 0, 0);

casacore::MDirection computePhaseDirection(double centralTime,
                                           MPosition reference) {
  casacore::MEpoch time = casacore::MEpoch(
      casacore::MVEpoch(centralTime / DAY_IN_SECONDS), casacore::MEpoch::UTC);
  casacore::MeasFrame frame(reference, time);

  const MDirection::Ref azel_reference(MDirection::AZEL, frame);
  const MDirection::Ref j2000_reference(MDirection::J2000, frame);
  MDirection azel_zenith(casacore::MVDirection(0.0, 0.0, 1.0), azel_reference);
  const MDirection phaseDirection =
      MDirection::Convert(azel_zenith, j2000_reference)();
  return MDirection::Convert(azel_zenith, j2000_reference)();
}

}  // namespace

namespace dp3::steps {

AartfaacReader::AartfaacReader(const common::ParameterSet& parset,
                               const std::string& prefix)
    : mode_(base::RcuMode::FromNumber(parset.getInt(prefix + "mode"))),
      antenna_configuration_(parset.getString(prefix + "antennafield")),
      filename_(
          parset.getString(prefix + "name", parset.getString("msin", ""))),
      n_times_(parset.getInt(prefix + "ntimes", 0)) {}

AartfaacReader::~AartfaacReader() = default;

void AartfaacReader::InitializeInfo() {
  AntennaConfig antenna_configuration =
      AntennaConfig(antenna_configuration_.c_str());
  antenna_positions_ = antenna_configuration.GetArrayFromMode(mode_);

  DPInfo dataset_info(metadata_.NrPolarizations(), metadata_.NrChannels(),
                      mode_.ToString());

  std::vector<double> channels(metadata_.NrChannels());
  std::vector<double> channels_width(metadata_.NrChannels());

  const double channel_width = metadata_.ChannelBandwidth();
  double start_frequency = input_->Frequency() + channel_width * 0.5;
  std::fill(channels_width.begin(), channels_width.end(), channel_width);
  for (size_t ch = 0; ch < channels.size(); ch++) {
    channels[ch] = start_frequency + channel_width * (0.5 + double(ch));
  }
  dataset_info.setChannels(std::move(channels), std::move(channels_width));
  antenna_names_.resize(metadata_.NrReceivers());
  antenna_diameters_.resize(metadata_.NrReceivers());
  std::fill(antenna_diameters_.begin(), antenna_diameters_.end(), 1);

  std::vector<int> antenna1(metadata_.NrBaselines());
  std::vector<int> antenna2(metadata_.NrBaselines());

  for (unsigned i = 0; i < antenna_names_.size(); i++) {
    antenna_names_[i] = std::string("A12") + "_" + std::to_string(i);
  }

  for (unsigned i = 0; i < antenna1.size(); i++) {
    const auto pair = dp3::common::ComputeBaseline(
        i, antenna_names_.size(), dp3::common::BaselineOrder::kColumnMajor);
    antenna1[i] = pair.first;
    antenna2[i] = pair.second;
  }

  integration_time_ = metadata_.EndTime() - metadata_.StartTime();
  dataset_info.setAntennas(antenna_names_, antenna_diameters_,
                           antenna_positions_, antenna1, antenna2);

  const unsigned int n_times = n_times_ == 0 ? input_->NTimesteps() : n_times_;
  phase_direction_ =
      computePhaseDirection(aartfaac::TimeToCasa(metadata_.StartTime()) +
                                integration_time_ * n_times * 0.5,
                            antenna_positions_[0]);
  dataset_info.setArrayInformation(antenna_positions_[0], phase_direction_,
                                   phase_direction_, phase_direction_);
  uvw_calculator_ = std::make_unique<base::UVWCalculator>(
      phase_direction_, antenna_positions_[0], antenna_positions_);

  dataset_info.setTimeIntervalAndSteps(integration_time_, n_times);
  GetWritableInfoOut() = dataset_info;
}

void AartfaacReader::ApplyDelayCorrection(DPBuffer& buffer, size_t baseline) {
  for (int ch = 0; ch < metadata_.NrChannels(); ch++) {
    const float frequency = getInfoIn().chanFreqs()[ch];
    const float w = buffer.GetUvw()(baseline, 2);
    const float angle = -2.0f * M_PI * frequency * w / aocommon::kSpeedOfLight;
    const float delay_sin = sinf(angle);
    const float delay_cos = cosf(angle);

    for (int pol = 0; pol < metadata_.NrPolarizations(); pol++) {
      auto& visibilities = buffer.GetData();

      visibilities(baseline, ch, pol) = std::complex<float>(
          visibilities(baseline, ch, pol).real() * delay_cos -
              visibilities(baseline, ch, pol).imag() * delay_sin,
          visibilities(baseline, ch, pol).real() * delay_sin +
              visibilities(baseline, ch, pol).imag() * delay_cos);
    }
  }
}

void AartfaacReader::CreateInitialSubtables() {
  out_name_ = (boost::filesystem::temp_directory_path() /
               boost::filesystem::unique_path("TEMP%%%%%%%.MS"))
                  .string();

  GetWritableInfoOut().setMsNames(out_name_, data_col_name_, flag_col_name_,
                                  weight_col_name_);
  std::unique_ptr<base::SubtableWriter> writer =
      std::make_unique<base::SubtableWriter>(out_name_, metadata_.NrChannels());

  start_time_ = metadata_.StartTime();
  const double start_time = aartfaac::TimeToCasa(metadata_.StartTime());
  const double end_time = aartfaac::TimeToCasa(metadata_.EndTime());

  CreateAntennaTable(start_time, *writer);
  CreateSpectralWindowTable(metadata_.NrChannels(),
                            metadata_.FirstChannelFrequency(),
                            metadata_.ChannelBandwidth(), *writer);
  CreateSourceTable(start_time, start_time, *writer);
  CreateFieldTable(start_time, *writer);
  CreateObservationTable(start_time, end_time, *writer);
  CreatePolarizationTable(*writer);
}

bool AartfaacReader::process(std::unique_ptr<DPBuffer> buffer) {
  if (getFieldsToRead().Data()) {
    buffer->GetData().resize({metadata_.NrBaselines(), metadata_.NrChannels(),
                              metadata_.NrPolarizations()});
  }
  if (getFieldsToRead().Flags()) {
    buffer->GetFlags().resize({metadata_.NrBaselines(), metadata_.NrChannels(),
                               metadata_.NrPolarizations()});
  }

  if (getFieldsToRead().Weights()) {
    buffer->GetWeights().resize({metadata_.NrBaselines(),
                                 metadata_.NrChannels(),
                                 metadata_.NrPolarizations()});
  }
  if (getFieldsToRead().Uvw()) {
    buffer->GetUvw().resize({metadata_.NrBaselines(), 3});
  }

  common::NSTimer::StartStop sstime(timer_);

  if ((input_->HasMore()) && (current_time_ < n_times_ || n_times_ == 0)) {
    Timestep time_step = input_->ReadTimestep(buffer->GetData().data());

    integration_time_ = time_step.endTime - time_step.startTime;
    casacore::MEpoch time_epoch = casacore::MEpoch(
        casacore::MVEpoch(time_step.startTime / DAY_IN_SECONDS),
        casacore::MEpoch::UTC);
    buffer->SetTime(time_step.startTime);
    buffer->SetExposure(integration_time_);
    buffer->GetWeights().fill(integration_time_ * metadata_.ChannelBandwidth());
    buffer->GetFlags().fill(false);
    for (size_t i = 0; i < metadata_.NrBaselines(); i++) {
      for (size_t k = 0; k < metadata_.NrChannels(); k++) {
        auto mask = !xt::isfinite(
            xt::sum(xt::real(xt::view(buffer->GetData(), i, k, xt::all()))));
        xt::view(buffer->GetFlags(), i, k, xt::all()) = mask;
      }
    }

    xt::filtration(buffer->GetData(), buffer->GetFlags()) =
        std::complex<float>(0.f, 0.f);

    // Calculate UVWs apply geometric delay and compute the weights
    const std::vector<int>& antenna1 = getInfoIn().getAnt1();
    const std::vector<int>& antenna2 = getInfoIn().getAnt2();
    for (size_t i = 0; i < metadata_.NrBaselines(); i++) {
      xt::view(buffer->GetUvw(), i, xt::all()) =
          xt::adapt(uvw_calculator_->getUVW(antenna2[i], antenna1[i],
                                            time_step.startTime));
      ApplyDelayCorrection(*buffer, i);
    }
    getNextStep()->process(std::move(buffer));
    current_time_++;
    return true;
  }
  end_time_ = input_->Header().EndTime();
  // TODO this should not be done every process() call
  GetWritableInfoOut().setTimes(start_time_, end_time_, integration_time_);
  return false;
}

void AartfaacReader::OpenFile() {
  if (std::filesystem::exists(filename_)) {
    input_ = std::make_unique<File>(filename_.c_str(), mode_);
  } else {
    throw std::invalid_argument("No such file: " + filename_);
  }
  metadata_ = input_->Header();
}

void AartfaacReader::CreateSpectralWindowTable(int nchannels,
                                               double startFrequency,
                                               double channelWidth,
                                               base::SubtableWriter& writer) {
  std::ostringstream band_name;
  band_name << "AARTF_BAND_" << round(1e-6 * input_->Frequency() * 10.0) / 10.0;
  std::vector<base::SubtableWriter::ChannelInfo> channels(nchannels);
  for (size_t ch = 0; ch != channels.size(); ++ch) {
    base::SubtableWriter::ChannelInfo& channel = channels[ch];
    channel.channel_frequency = startFrequency + channelWidth * ch;
    channel.channel_width = channelWidth;
    channel.effective_bandwidth = channelWidth;
    channel.resolution = channelWidth;
  }
  writer.WriteBandInfo(band_name.str(), channels, startFrequency, channelWidth,
                       false);
}

void AartfaacReader::CreateAntennaTable(double time,
                                        base::SubtableWriter& writer) {
  // TODO add the possiblity to select the antennas
  std::vector<base::SubtableWriter::AntennaInfo> antennas(
      antenna_positions_.size());
  for (size_t ant = 0; ant != antennas.size(); ++ant) {
    base::SubtableWriter::AntennaInfo& antenna_info = antennas[ant];
    antenna_info.name = antenna_names_[ant];
    antenna_info.station = "AARTFAAC";
    antenna_info.type = "GROUND-BASED";
    // CASA does not support FIXED so we use ALT-AZ instead
    antenna_info.mount = "ALT-AZ";  //
    antenna_info.x = antenna_positions_[ant].getValue()(0);
    antenna_info.y = antenna_positions_[ant].getValue()(1);
    antenna_info.z = antenna_positions_[ant].getValue()(2);
    antenna_info.diameter = antenna_diameters_[ant];
    antenna_info.flag = false;
  }

  writer.WriteAntennas(antennas, axes_, time);
}

void AartfaacReader::CreateSourceTable(double startTime, double endTime,
                                       base::SubtableWriter& writer) {
  base::SubtableWriter::SourceInfo source;

  source.source_id = 0;
  source.time = startTime;
  source.interval = endTime;
  source.spectral_window_id = 0;
  source.num_lines = 0;
  source.name = "AARTFAAC";
  source.calibration_group = 0;
  source.code = "";
  double ra = phase_direction_.getAngle().getValue()[0];
  double dec = phase_direction_.getAngle().getValue()[1];
  source.ra = ra;  // (in radians)
  source.dec = dec;
  source.proper_motion[0] = 0.0;
  source.proper_motion[1] = 0.0;
  writer.WriteSource(source);
}

void AartfaacReader::CreateFieldTable(double time,
                                      base::SubtableWriter& writer) {
  base::SubtableWriter::FieldInfo field;
  field.name = "AARTFAAC";
  field.code = std::string();
  field.time = time;
  field.num_poly = 0;
  double ra = phase_direction_.getAngle().getValue()[0];
  double dec = phase_direction_.getAngle().getValue()[1];
  field.delay_direction_ra = ra;  // (in radians)
  field.delay_direction_dec = dec;
  field.phase_direction_ra = field.delay_direction_ra;
  field.phase_direction_dec = field.delay_direction_dec;
  field.reference_direction_ra = field.delay_direction_ra;
  field.reference_direction_dec = field.delay_direction_dec;
  field.source_id = -1;
  field.flag_row = false;
  writer.WriteField(field);
}

void AartfaacReader::CreateObservationTable(double startTime, double endTime,
                                            base::SubtableWriter& writer) {
  base::SubtableWriter::ObservationInfo observation;
  observation.telescope_name = "AARTFAAC";
  observation.start_time = startTime;
  observation.end_time = endTime;
  observation.observer = "Unknown";
  observation.schedule_type = "AARTFAAC";
  observation.project = "Unknown";
  observation.release_date = 0;
  observation.flag_row = false;

  observation.antenna_type = mode_.AntennaType();
  observation.rcu_mode = mode_.mode;
  observation.flag_window_size = 0;
  writer.WriteObservation(observation);
}

void AartfaacReader::CreatePolarizationTable(base::SubtableWriter& writer) {
  writer.WriteLinearPolarizations(false);
}

std::string AartfaacReader::msName() const { return out_name_; }

void AartfaacReader::updateInfo(const base::DPInfo& infoIn) {
  OpenFile();
  InitializeInfo();
  CreateInitialSubtables();
}

void AartfaacReader::show(std::ostream& os) const { os << "AARTFAACInput\n"; }

void AartfaacReader::finish() { getNextStep()->finish(); }

}  // namespace dp3::steps
