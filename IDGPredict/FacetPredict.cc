// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#include "FacetPredict.h"

#ifdef HAVE_IDG

#include "DS9FacetFile.h"
#include "FitsWriter.h"
#include "IDGConfiguration.h"
#include "../DPPP/DPInput.h"

#include <aocommon/uvector.h>
#include <aocommon/banddata.h>

#include <iostream>

namespace DP3 {
namespace DPPP {

FacetPredict::FacetPredict(DPInput& input,
                           const std::vector<std::string>& fits_model_files,
                           const std::string& ds9_regions_file)
    : buffer_size_(0), input_(input), info_(), ant1_(), ant2_(), timer_() {
  if (fits_model_files.empty()) {
    throw std::runtime_error("No fits files specified for IDG predict");
  }
  readers_.reserve(fits_model_files.size());
  for (const std::string& file : fits_model_files) readers_.emplace_back(file);

  const size_t full_width = readers_.front().ImageWidth();
  const size_t full_height = readers_.front().ImageHeight();
  ref_frequency_ = readers_.front().Frequency();
  pixel_size_x_ = readers_.front().PixelSizeX();
  pixel_size_y_ = readers_.front().PixelSizeY();
  std::vector<aocommon::UVector<double>> models(readers_.size());
  for (size_t img = 0; img != readers_.size(); ++img) {
    if (readers_[img].ImageWidth() != full_width ||
        readers_[img].ImageHeight() != full_height)
      throw std::runtime_error("Image for spectral term " +
                               std::to_string(img) +
                               " has inconsistent dimensions");
    if (readers_[img].PixelSizeX() != pixel_size_x_ ||
        readers_[img].PixelSizeY() != pixel_size_y_)
      throw std::runtime_error("Pixel size of spectral term " +
                               std::to_string(img) +
                               " is inconsistent with first spectral term");
    models[img].resize(full_width * full_height);
    readers_[img].Read(models[img].data());
  }

  DS9FacetFile facet_file(ds9_regions_file);
  std::vector<Facet> facets = facet_file.Read(
      readers_.front().PhaseCentreRA(), readers_.front().PhaseCentreDec(),
      pixel_size_x_, pixel_size_y_, full_width, full_height);
  std::cout << "Read " << facets.size() << " facet definitions.\n";
  if (facets.empty()) {
    throw std::runtime_error("No facet definitions found in " +
                             ds9_regions_file);
  }

  size_t area = 0;
  directions_.reserve(facets.size());
  images_.reserve(facets.size());
  for (const Facet& facet : facets) {
    std::cout << "Facet: Ra,Dec: " << facet.RA() << "," << facet.Dec()
              << " Vertices:";
    for (const Vertex& v : facet) std::cout << " (" << v.x << "," << v.y << ")";
    std::cout << std::endl;

    directions_.emplace_back(facet.RA(), facet.Dec());
    images_.emplace_back();
    FacetImage& image = images_.back();
    constexpr bool kMakeSquare = true;  // only necessary for IDG though
    constexpr double kPadding = 1.0;
    image.CopyFacetPart(facet, models, full_width, full_height, kPadding,
                        kMakeSquare);
    area += image.Width() * image.Height();
  }
  std::cout << "Area covered: " << area / 1024 << " Kpixels^2\n";
}

void FacetPredict::updateInfo(const DP3::DPPP::DPInfo& info) {
  if (info.ncorr() != 4) {
    throw std::invalid_argument("FacetPredict only supports 4 correlations.");
  }

  // TODO, when FacetPredict is a DPStep:
  // Remove info_ member, and call DPStep::updateInfo(info); here.
  info_ = info;

  // Generate ant1_ and ant2_, which contain mapped antenna indices.
  ant1_.clear();
  ant1_.reserve(info.getAnt1().size());
  for (size_t ant : info.getAnt1()) {
    assert(info.antennaMap()[ant] >= 0);
    ant1_.push_back(info.antennaMap()[ant]);
  }

  ant2_.clear();
  ant2_.reserve(info.getAnt2().size());
  for (size_t ant : info.getAnt2()) {
    assert(info.antennaMap()[ant] >= 0);
    ant2_.push_back(info.antennaMap()[ant]);
  }

  // Without a factor 1.5 (instead of 1.0, below), some baselines did not
  // get visibilities from the CPU-optimized IDG version.
  // TODO (AST-223): Determine the logic why and how max_baseline_, which is in
  // meters, depends on the inverse of the pixel size, which is in radians.
  max_baseline_ = 1.0 / std::min(pixel_size_x_, pixel_size_y_);
  max_w_ = max_baseline_ * 0.1;
  std::cout << "Predicting baselines up to " << max_baseline_
            << " wavelengths.\n";
}

bool FacetPredict::IsStarted() const { return !buffersets_.empty(); }

void FacetPredict::StartIDG(bool saveFacets) {
  buffersets_.clear();
  idg::api::Type proxy_type = idg::api::Type::CPU_OPTIMIZED;
  int buffersize = 0;
  size_t n_terms = readers_.size();

  idg::api::options_type options;
  IdgConfiguration::Read(proxy_type, buffersize, options);
  std::vector<aocommon::UVector<double>> data(n_terms);
  meta_data_.clear();
  FitsReader& reader = readers_.front();
  for (FacetImage& img : images_) {
    // dl and dm indicate the relative difference of the center pixel of the
    // facet to the center pixel of the fits image.
    // dl is positive if rA of the facet is larger than the rA of the image.
    // dm is positive if dec of the facet is larger than the dec of the image.
    // Note that rA decreases if the x coordinate increases, while dec
    // increases of the y coordinates increases.
    const double dl = (int(reader.ImageWidth() / 2) -
                       (img.OffsetX() + int(img.Width() / 2))) *
                      pixel_size_x_;
    const double dm = (int(reader.ImageHeight() / 2) -
                       (img.OffsetY() + int(img.Height() / 2))) *
                      pixel_size_y_;
    const double dp = sqrt(1.0 - dl * dl - dm * dm) - 1.0;
    std::cout << "Initializing gridder " << buffersets_.size() << " ("
              << img.Width() << " x " << img.Height() << ", +" << img.OffsetX()
              << "," << img.OffsetY() << ", dl=" << dl * 180.0 / M_PI
              << " deg, dm=" << dm * 180.0 / M_PI << " deg)\n"
              << std::endl;

    // TODO make full polarization
    for (size_t term = 0; term != n_terms; ++term) {
      data[term].assign(img.Width() * img.Height() * 4, 0.0);
      std::copy(img.Data(term), img.Data(term) + img.Width() * img.Height(),
                data[term].data());

      buffersets_.emplace_back(idg::api::BufferSet::create(proxy_type));
      idg::api::BufferSet& bs = *buffersets_.back();
      options["padded_size"] = size_t(1.2 * img.Width());
      // options["max_threads"] = int(1);
      bs.init(img.Width(), pixel_size_x_, max_w_ + 1.0, dl, dm, dp, options);
      bs.set_image(data[term].data());
      bs.init_buffers(0, {info_.chanFreqs()}, info_.antennaUsed().size(),
                      max_baseline_, options,
                      idg::api::BufferSetType::kBulkDegridding);
    }

    if (saveFacets) {
      FitsWriter writer;
      writer.SetImageDimensions(img.Width(), img.Height(),
                                reader.PhaseCentreRA(), reader.PhaseCentreDec(),
                                pixel_size_x_, pixel_size_y_);
      writer.SetPhaseCentreShift(dl, dm);
      writer.Write("facet" + std::to_string(meta_data_.size()) + ".fits",
                   img.Data(0));
    }

    meta_data_.emplace_back(dl, dm, dp);
  }
}

std::vector<DPBuffer> FacetPredict::Predict(const size_t direction) {
  const size_t n_terms = readers_.size();
  const size_t term_size =
      buffers_.size() * info_.nbaselines() * info_.nchan() * info_.ncorr();

  // term_data indexing: [term-1][timestep][baseline][channel][correlation].
  aocommon::UVector<std::complex<float>> term_data((n_terms - 1) * term_size);

  // Initialize a vector with uvw pointers. Raise max_w_ if needed.
  std::vector<const double*> uvws = InitializeUVWs();

  // Initialize output buffers and fill them with IDG predictions. */
  std::vector<DPBuffer> result =
      ComputeVisibilities(direction, uvws, term_data.data());

  // Correct phase shift for the computed visibilities.
  // Apply polynomial-term corrections and add all to values of 'term 0'
  CorrectVisibilities(uvws, result, term_data.data(), direction);

  return result;
}

std::vector<const double*> FacetPredict::InitializeUVWs() {
  std::vector<const double*> uvws;
  uvws.reserve(buffers_.size());

  double old_max_w = max_w_;
  const double max_baseline2 = max_baseline_ * max_baseline_;

  for (DPBuffer& buffer : buffers_) {
    uvws.push_back(input_.fetchUVW(buffer, buffer, timer_).data());

    const double* uvw = uvws.back();
    for (std::size_t bl = 0; bl < info_.nbaselines(); ++bl) {
      if (uvw[2] > max_w_) {
        double uvwr2 = uvw[0] * uvw[0] + uvw[1] * uvw[1] + uvw[2] * uvw[2];
        if (uvwr2 <= max_baseline2) {
          max_w_ = uvw[2];
        }
      }
      uvw += 3;
    }
  }

  if (max_w_ > old_max_w) {
    // Increase max_w by 10% for preventing repeated IDG initialization calls.
    max_w_ *= 1.1;
    std::cout << "Increasing maximum w to " << max_w_ << '\n';
    StartIDG(false);
  }

  return uvws;
}

std::vector<DPBuffer> FacetPredict::ComputeVisibilities(
    const size_t direction, const std::vector<const double*>& uvws,
    std::complex<float>* term_data) const {
  const size_t n_timesteps = buffers_.size();
  const size_t n_terms = readers_.size();
  const size_t timestep_size =
      info_.nbaselines() * info_.nchan() * info_.ncorr();
  const size_t term_size = n_timesteps * timestep_size;

  std::vector<DPBuffer> result;
  result.reserve(buffers_.size());
  for (const DPBuffer& buffer : buffers_) {
    result.emplace_back(buffer);
    result.back().setData(casacore::Cube<casacore::Complex>(
        info_.ncorr(), info_.nchan(), info_.nbaselines()));
  }

  // IDG uses a flipped coordinate system
  static const double kUVWFactors[3] = {1.0, -1.0, -1.0};

  std::vector<std::complex<float>*> data_ptrs;
  data_ptrs.reserve(n_timesteps);

  // Compute visibilities for all terms and timesteps.
  for (size_t term = 0; term < n_terms; ++term) {
    data_ptrs.clear();
    if (term == 0) {  // data_ptrs point directly to the buffers.
      for (size_t t = 0; t < n_timesteps; ++t) {
        data_ptrs.push_back(result[t].getData().data());
      }
    } else {  // Use term_data as auxiliary buffers for the other terms.
      std::complex<float>* data = term_data + (term - 1) * term_size;
      for (size_t t = 0; t < n_timesteps; ++t) {
        data_ptrs.push_back(data);
        data += timestep_size;
      }
    }

    idg::api::BufferSet& bs = *buffersets_[direction * n_terms + term];
    // Since we initialize the BufferSet with a single band, the index is 0.
    constexpr int kDegridderIndex = 0;
    bs.get_bulk_degridder(kDegridderIndex)
        ->compute_visibilities(ant1_, ant2_, uvws, data_ptrs, kUVWFactors);
  }

  return result;
}

double FacetPredict::ComputePhaseShiftFactor(const double* uvw,
                                             const size_t direction) const {
  // The phase shift angle is:
  // 2*pi * (u*dl - v*dm - w*dp) * (frequency / speed_of_light)
  // The v and w terms are negated since IDG uses a flipped coordinate system.
  // The computed phase shift factor contains all parts except the frequency.
  constexpr double two_pi_div_c = 2.0 * M_PI / aocommon::c();
  return (uvw[0] * meta_data_[direction].dl -
          uvw[1] * meta_data_[direction].dm -
          uvw[2] * meta_data_[direction].dp) *
         two_pi_div_c;
}

void FacetPredict::CorrectVisibilities(const std::vector<const double*>& uvws,
                                       std::vector<DPBuffer>& result,
                                       const std::complex<float>* term_data,
                                       const size_t direction) const {
  const size_t n_timesteps = buffers_.size();
  const size_t n_terms = readers_.size();
  const size_t timestep_size =
      info_.nbaselines() * info_.nchan() * info_.ncorr();
  const size_t term_size = n_timesteps * timestep_size;

  // The polynomial term corrections use the "polynomial spectrum" definition,
  // equal to the one e.g. used by WSClean in component outputs (see
  // https://sourceforge.net/p/wsclean/wiki/ComponentList/ ) and in text
  // files when 'logarithmic SI' is false. The definition is:
  //   S(nu) = term0 + term1 (nu/refnu - 1) + term2 (nu/refnu - 1)^2 + ...

  // Create pointers to the visibilities for the extra terms.
  std::vector<const std::complex<float>*> term_data_ptrs;
  term_data_ptrs.reserve(n_terms - 1);
  for (size_t term = 0; term < n_terms - 1; ++term) {
    term_data_ptrs.push_back(term_data + term * term_size);
  }

  // Precompute polynomial frequency factors for all channels.
  std::vector<double> freq_factors;
  if (n_terms > 1) {
    freq_factors.reserve(info_.nchan());
    for (double freq : info_.chanFreqs()) {
      freq_factors.push_back(freq / ref_frequency_ - 1.0);
    }
  }

  for (size_t t = 0; t < n_timesteps; ++t) {
    const double* uvw = uvws[t];
    std::complex<float>* result_data = result[t].getData().data();
    for (size_t bl = 0; bl < info_.nbaselines(); ++bl) {
      const double phase_factor = ComputePhaseShiftFactor(uvw, direction);
      uvw += 3;

      for (size_t ch = 0; ch < info_.nchan(); ++ch) {
        const double angle = phase_factor * info_.chanFreqs()[ch];
        const std::complex<float> phasor(cos(angle), sin(angle));
        result_data[0] *= phasor;
        result_data[1] *= phasor;
        result_data[2] *= phasor;
        result_data[3] *= phasor;

        double polynomial_factor = 1.0;
        for (auto& term_data_ptr : term_data_ptrs) {
          polynomial_factor *= freq_factors[ch];
          // Merge the phase shift correction and the polynomial factor.
          const std::complex<float> term_phasor =
              phasor * float(polynomial_factor);
          result_data[0] += term_data_ptr[0] * term_phasor;
          result_data[1] += term_data_ptr[1] * term_phasor;
          result_data[2] += term_data_ptr[2] * term_phasor;
          result_data[3] += term_data_ptr[3] * term_phasor;
          term_data_ptr += 4;
        }
        result_data += 4;
      }
    }
  }
}

const std::vector<std::pair<double, double>>& FacetPredict::GetDirections()
    const {
  return directions_;
}

void FacetPredict::SetBufferSize(size_t n_timesteps) {
  buffer_size_ = n_timesteps;
}

}  // namespace DPPP
}  // namespace DP3

#else  // HAVE_IDG

namespace {

void notCompiled() {
  throw std::runtime_error(
      "Facet prediction is not available, because DP3 was not compiled with "
      "IDG support");
}

}  // namespace

namespace DP3 {
namespace DPPP {

FacetPredict::FacetPredict(const std::vector<std::string>&, const std::string&,
                           PredictCallback&&) {
  notCompiled();
}

void FacetPredict::updateInfo(const DP3::DPPP::DPInfo&) { notCompiled(); }

bool FacetPredict::IsStarted() const {
  notCompiled();
  return false;
}

void FacetPredict::StartIDG(bool) { notCompiled(); }

void FacetPredict::RequestPredict(size_t, size_t, size_t, size_t, size_t,
                                  size_t, const double*) {
  notCompiled();
}

const std::vector<std::pair<double, double>>& FacetPredict::GetDirections()
    const {
  notCompiled();
  return std::vector<std::pair<double, double>>();
}

void FacetPredict::Flush(size_t) { notCompiled(); }

void FacetPredict::SetBufferSize(size_t) { notCompiled(); }

}  // namespace DPPP
}  // namespace DP3

#endif  // HAVE_IDG
