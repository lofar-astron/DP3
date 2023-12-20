// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IDGImager.h"

#include <iostream>

#include <sstream>
#include <base/FlagCounter.h>
#include <base/IDGConfiguration.h>

#ifdef HAVE_IDG
#include <idg-api.h>
#endif
#include <xtensor/xcomplex.hpp>
#include <xtensor/xindex_view.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>
#include <aocommon/staticfor.h>
#include <aocommon/fits/fitswriter.h>
#include <casacore/casa/BasicSL/Constants.h>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3::steps {

IDGImager::IDGImager(const common::ParameterSet& parset,
                     const std::string& prefix)
    : name_(prefix),
      grid_size_{parset.getUint(prefix + "image_size", 1024)},
      dl_{parset.getDouble(prefix + "dl", 0.)},
      dm_{parset.getDouble(prefix + "dm", 0.)},
      image_name_{parset.getString(prefix + "image_name", "image_t%t.fits")} {
  pixel_scale_ = DegToRad(parset.getDouble(prefix + "scale", -1));
  pixel_scale_ = pixel_scale_ > 0 ? pixel_scale_ : M_PI / float(grid_size_);
#ifdef HAVE_IDG
  std::string proxy_type =
      parset.getString(prefix + "proxy_type", "CPU_OPTIMIZED");
  proxy_ = ProxyTypeFromString(proxy_type);
  buffer_set_.reset(idg::api::BufferSet::create(proxy_));
  ConfigureIDG();
#else
  throw std::runtime_error(
      "Cannot use the IDG imager without IDG. Please recompile DP3 with IDG "
      "support.");
#endif
}

void IDGImager::updateInfo(const DPInfo& info_in) { Step::updateInfo(info_in); }

void IDGImager::show(std::ostream& os) const {}

void IDGImager::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  dp3::base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " IDGImager " << name_ << '\n';
}

void IDGImager::ConfigureIDG() {
#ifdef HAVE_IDG
  IdgConfiguration::Read(proxy_, buffer_size_, options_);
#endif
}

void IDGImager::InitBuffers(const double max_w, const double max_baseline) {
#ifdef HAVE_IDG
  size_t kNumberOfTimeSteps = 1;
  buffer_set_->init(grid_size_, pixel_scale_, max_w, dl_, dm_, options_);
  buffer_set_->init_buffers(kNumberOfTimeSteps, {info().chanFreqs()},
                            info().antennaUsed().size(), max_baseline, options_,
                            idg::api::BufferSetType::kGridding);
  init_buffers_ = true;
#endif
}

bool IDGImager::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  const dp3::base::DPBuffer::UvwType uvw = buffer->GetUvw();
  dp3::base::DPBuffer::DataType visibilities = buffer->GetData();
  dp3::base::DPBuffer::WeightsType weights = buffer->GetWeights();
  const dp3::base::DPBuffer::FlagsType flags = buffer->GetFlags();
  const std::vector<int>& antenna1 = info().getAnt1();
  const std::vector<int>& antenna2 = info().getAnt2();

  xt::xtensor<double, 3> image_data({4, grid_size_, grid_size_});

  const double n_visibilities = visibilities.shape(0);

  const std::vector<double>& frequencies = info().chanFreqs();

  // W index
  const ushort kWIndex = 2;
  const float max_w = xt::amax(xt::view(uvw, xt::all(), kWIndex))() *
                      info().refFreq() / casacore::C::c;
  const float max_baseline =
      xt::amax(xt::sqrt(xt::pow(xt::view(uvw, xt::all(), 0), 2) +
                        xt::pow(xt::view(uvw, xt::all(), 1), 2) +
                        xt::pow(xt::view(uvw, xt::all(), 2), 2)))();

  PrepareVisibilities(visibilities, weights, flags);
  if (!init_buffers_) InitBuffers(max_w, max_baseline);

  WeightVisibilities(uvw, weights, frequencies, grid_size_, visibilities);
#ifdef HAVE_IDG
  const size_t kSingleTimeStepGridder = 0;
  idg::api::GridderBuffer& gridder =
      *buffer_set_->get_gridder(kSingleTimeStepGridder);

  buffer_set_->set_image(image_data.data());
  CopyVisibilitiesToBuffer(gridder, antenna1, antenna2, uvw, weights,
                           visibilities);
  gridder.finished();

  buffer_set_->get_image(image_data.data());
#endif

  // Normalize fluxes based on the number of total visibilities
  image_data /= n_visibilities;
  // Flip the image to have declination increasing from
  // bottom to top
  const ushort kYaxis = 1;
  image_data = xt::flip(image_data, kYaxis);

  WriteFITS(FormatName(image_name_, current_image_idx_), info(),
            buffer->GetTime(), dl_, dm_, pixel_scale_, image_data);
  timer_.stop();
  current_image_idx_++;
  getNextStep()->process(std::move(buffer));
  return true;
}

#ifdef HAVE_IDG
idg::api::Type IDGImager::ProxyTypeFromString(const std::string& idg_type) {
  if (idg_type.compare("CPU_OPTIMIZED") == 0) {
    return idg::api::Type::CPU_OPTIMIZED;
  } else if (idg_type.compare("CUDA_GENERIC") == 0) {
    return idg::api::Type::CUDA_GENERIC;
  } else if (idg_type.compare("CPU_REFERENCE") == 0) {
    return idg::api::Type::CPU_REFERENCE;
  } else if (idg_type.compare("HYBRID") == 0) {
    return idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED;
  } else {
    throw std::runtime_error("Unsupported IDG proxy type: " + idg_type);
  }
}
#endif

double IDGImager::DegToRad(double angle) { return angle * M_PI / 180.; }

#ifdef HAVE_IDG
void IDGImager::CopyVisibilitiesToBuffer(
    idg::api::GridderBuffer& buffer, const std::vector<int>& antenna1,
    const std::vector<int>& antenna2, const DPBuffer::UvwType& uvw,
    const DPBuffer::WeightsType& weights,
    const DPBuffer::DataType& visibilities) {
  aocommon::StaticFor<size_t> loop;
  loop.Run(0, visibilities.shape(0), [&](size_t start, size_t end, size_t) {
    for (size_t vis_idx = start; vis_idx < end; vis_idx++) {
      if (antenna1[vis_idx] != antenna2[vis_idx]) {
        buffer.grid_visibilities(0, antenna1[vis_idx], antenna2[vis_idx],
                                 &uvw(vis_idx, 0), &visibilities(vis_idx, 0, 0),
                                 &weights(vis_idx, 0, 0));
      }
    }
  });
}
#endif

void IDGImager::PrepareVisibilities(DPBuffer::DataType& visibilities,
                                    DPBuffer::WeightsType& weights,
                                    const DPBuffer::FlagsType& flags) {
  xt::filtration(visibilities, flags) = std::complex<float>(0.f, 0.f);
  xt::filtration(visibilities, xt::isnan(visibilities)) =
      std::complex<float>(0.f, 0.f);
  xt::filtration(weights, xt::isnan(weights)) = 0.f;
}

void IDGImager::WriteFITS(const std::string file_path,
                          const dp3::base::DPInfo& info, double time, double dl,
                          double dm, double pixel_scale,
                          xt::xtensor<double, 3> image) {
  const size_t n_polarizations = image.shape(0);
  const size_t image_size = image.shape(1);
  aocommon::FitsWriter writer;
  const size_t kNumFrequenciesPerImage = 1;
  writer.SetOrigin("DP3/IDG", " Image Domain Gridder");
  writer.SetImageDimensions(
      image_size, image_size, info.phaseCenter().getAngle().getValue()[0],
      info.phaseCenter().getAngle().getValue()[1], pixel_scale, pixel_scale);
  writer.AddExtraDimension(
      aocommon::FitsWriter::DimensionType::FrequencyDimension,
      kNumFrequenciesPerImage);

  writer.AddExtraDimension(
      aocommon::FitsWriter::DimensionType::PolarizationDimension,
      n_polarizations);

  const double kSecondsInDay = 3600. * 24.;
  writer.SetDate(time / kSecondsInDay);
  writer.SetFrequency(info.refFreq(), info.totalBW());
  writer.SetPhaseCentreShift(dl, dm);
  writer.StartMulti(file_path);

  for (ushort pol = 0u; pol < n_polarizations; pol++) {
    writer.AddToMulti(&image(pol, 0, 0));
  }

  writer.FinishMulti();
}

std::string IDGImager::FormatName(std::string name, int time) {
  std::string kTimeFmt = "%t";
  size_t pos = name.find(kTimeFmt);

  if (pos != std::string::npos) {
    name.replace(pos, kTimeFmt.size(), std::to_string(time));
    return name;
  }
  return name;
}

xt::xtensor<float, 2> IDGImager::ComputeWeightMap(
    size_t grid_size, const DPBuffer::WeightsType& weights,
    const xt::xtensor<size_t, 3>& uv_pixels) {
  const size_t n_baselines = weights.shape(0);
  const size_t n_channels = weights.shape(1);
  const size_t n_correlations = weights.shape(2);
  xt::xtensor<float, 2> weight_map({grid_size, grid_size});

  xt::view(weight_map, xt::all(), xt::all()) = 0.0f;

  for (size_t bl = 0u; bl < n_baselines; ++bl) {
    for (size_t chan = 0u; chan < n_channels; ++chan) {
      const size_t x = uv_pixels(0, bl, chan);
      const size_t y = uv_pixels(1, bl, chan);
      const float weight_xx = weights(bl, chan, 0);
      const float weight_yy = weights(bl, chan, n_correlations - 1);
      const float weight = (weight_xx + weight_yy);
      weight_map(y, x) += weight;
    }
  }

  xt::filtration(weight_map, xt::equal(weight_map, 0.0f)) = 1.0f;
  weight_map = 1.0f / weight_map;
  return weight_map;
}

void IDGImager::ApplyWeights(const xt::xtensor<float, 2>& weight_map,
                             const xt::xtensor<size_t, 3>& uv_pixels,
                             dp3::base::DPBuffer::DataType& visibilities) {
  const size_t n_baselines = visibilities.shape(0);
  const size_t n_channels = visibilities.shape(1);
  const size_t n_correlations = visibilities.shape(2);
  aocommon::StaticFor<size_t> loop;

  for (size_t bl = 0; bl < n_baselines; bl++) {
    for (size_t chan = 0; chan < n_channels; ++chan) {
      for (size_t cor = 0; cor < n_correlations; ++cor) {
        const size_t x = uv_pixels(0, bl, chan);
        const size_t y = uv_pixels(1, bl, chan);
        visibilities(bl, chan, cor) *= weight_map(y, x);
      }
    }
  }
}

void IDGImager::ComputeUVPixels(const dp3::base::DPBuffer::UvwType& uvw,
                                const float uv_max,
                                const std::vector<double>& frequencies,
                                xt::xtensor<size_t, 3>& uv_pixels,
                                size_t grid_size_) {
  const size_t n_baselines = uvw.shape(0);
  const size_t n_channels = frequencies.size();

  for (size_t baseline_idx = 0; baseline_idx < n_baselines; baseline_idx++) {
    for (size_t channel_idx = 0; channel_idx < n_channels; channel_idx++) {
      uv_pixels(0, baseline_idx, channel_idx) =
          roundf((uvw(baseline_idx, 0) / uv_max + 1) * (grid_size_ / 2.f));
      uv_pixels(1, baseline_idx, channel_idx) =
          roundf((uvw(baseline_idx, 1) / uv_max + 1) * (grid_size_ / 2.f));
    }
  }
}

void IDGImager::WeightVisibilities(const DPBuffer::UvwType& uvw,
                                   const DPBuffer::WeightsType& weights,
                                   const std::vector<double> frequencies,
                                   const size_t grid_size,
                                   DPBuffer::DataType& visibilities) {
  const float max_uv =
      xt::amax(xt::sqrt(xt::pow(xt::view(uvw, xt::all(), 0), 2) +
                        xt::pow(xt::view(uvw, xt::all(), 1), 2)))();
  xt::xtensor<size_t, 3> uv_pixels({2, uvw.shape(0), frequencies.size()});
  ComputeUVPixels(uvw, max_uv, frequencies, uv_pixels, grid_size);

  visibilities *= weights;
  const xt::xtensor<float, 2>& weight_map =
      ComputeWeightMap(grid_size, weights, uv_pixels);
  ApplyWeights(weight_map, uv_pixels, visibilities);
}

void IDGImager::finish() { getNextStep()->finish(); }

}  // namespace dp3::steps
