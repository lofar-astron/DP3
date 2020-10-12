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

#include <iostream>

namespace DP3 {
namespace DPPP {

FacetPredict::FacetPredict(DPInput& input,
                           const std::vector<std::string>& fits_model_files,
                           const std::string& ds9_regions_file,
                           PredictCallback&& callback)
    : predict_callback_(std::move(callback)),
      padding_(1.0),
      buffer_size_(0),
      input_(input),
      info_(),
      ant1_(),
      ant2_(),
      timer_() {
  if (fits_model_files.empty()) {
    throw std::runtime_error("No fits files specified for IDG predict");
  }
  readers_.reserve(fits_model_files.size());
  for (const std::string& file : fits_model_files) readers_.emplace_back(file);

  full_width_ = readers_.front().ImageWidth();
  full_height_ = readers_.front().ImageHeight();
  ref_frequency_ = readers_.front().Frequency();
  pixel_size_x_ = readers_.front().PixelSizeX();
  pixel_size_y_ = readers_.front().PixelSizeY();
  std::vector<aocommon::UVector<double>> models(readers_.size());
  for (size_t img = 0; img != readers_.size(); ++img) {
    if (readers_[img].ImageWidth() != full_width_ ||
        readers_[img].ImageHeight() != full_height_)
      throw std::runtime_error("Image for spectral term " +
                               std::to_string(img) +
                               " has inconsistent dimensions");
    if (readers_[img].PixelSizeX() != pixel_size_x_ ||
        readers_[img].PixelSizeY() != pixel_size_y_)
      throw std::runtime_error("Pixel size of spectral term " +
                               std::to_string(img) +
                               " is inconsistent with first spectral term");
    models[img].resize(full_width_ * full_height_);
    readers_[img].Read(models[img].data());
  }

  DS9FacetFile facet_file(ds9_regions_file);
  std::vector<Facet> facets = facet_file.Read(
      readers_.front().PhaseCentreRA(), readers_.front().PhaseCentreDec(),
      pixel_size_x_, pixel_size_y_, full_width_, full_height_);
  std::cout << "Read " << facets.size() << " facet definitions.\n";
  if (facets.empty()) {
    throw std::runtime_error("No facet definitions found in " +
                             ds9_regions_file);
  }

  bool make_square = true;  // only necessary for IDG though
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
    image.CopyFacetPart(facet, models, full_width_, full_height_, padding_,
                        make_square);
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
  size_t n_terms = readers_.size();

  size_t max_channels = info_.chanFreqs().size();
  long int pageCount = sysconf(_SC_PHYS_PAGES),
           pageSize = sysconf(_SC_PAGE_SIZE);
  int64_t memory = (int64_t)pageCount * (int64_t)pageSize;
  uint64_t memPerTimestep = idg::api::BufferSet::get_memory_per_timestep(
      info_.antennaUsed().size(), max_channels);
  memPerTimestep *= 2;  // IDG uses two internal buffer
  // Allow the directions together to use 1/4th of the available memory for
  // the vis buffers.
  size_t allocatableTimesteps =
      memory / 4 / images_.size() / n_terms / memPerTimestep;
  // TODO once a-terms are supported, this should include the size required
  // for the a-terms.
  std::cout << "Allocatable timesteps per direction: " << allocatableTimesteps
            << '\n';

  int buffersize = std::max(allocatableTimesteps, size_t(1));
  if (buffer_size_ != 0) {
    buffersize = buffer_size_;
    std::cout << "Buffer size manually set to " << buffersize << " timesteps\n";
  }
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
    const double dm = (img.OffsetY() + int(img.Height() / 2) -
                       int(reader.ImageHeight() / 2)) *
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
      bs.init_buffers(buffersize, {info_.chanFreqs()},
                      info_.antennaUsed().size(), max_baseline_, options,
#ifdef BULK_DEGRIDDING
                      idg::api::BufferSetType::kBulkDegridding);
#else
                      idg::api::BufferSetType::degridding);
#endif
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

    meta_data_.emplace_back();
    FacetMetaData& m = meta_data_.back();
    m.dl = dl;
    m.dm = dm;
    m.dp = dp;
    m.is_initialized = false;
    m.row_id_offset = 0;
  }
}

void FacetPredict::RequestPredict(size_t direction, size_t data_desc_id,
                                  size_t rowId, size_t timeIndex,
                                  size_t antenna1, size_t antenna2,
                                  const double* uvw) {
  size_t n_terms = readers_.size();
  double uvwr2 = uvw[0] * uvw[0] + uvw[1] * uvw[1] + uvw[2] * uvw[2];
  if (uvw[2] > max_w_ && uvwr2 <= max_baseline_ * max_baseline_) {
    Flush(data_desc_id);
    max_w_ = std::max(uvw[2], max_w_ * 1.5);
    std::cout << "Increasing maximum w to " << max_w_ << '\n';
    StartIDG(false);
  }
  for (size_t termIndex = 0; termIndex != n_terms; ++termIndex) {
    idg::api::BufferSet& bs = *buffersets_[direction * n_terms + termIndex];
    FacetMetaData& meta = meta_data_[direction];
    if (!meta.is_initialized) {
      meta.row_id_offset = rowId;
      meta.is_initialized = true;
    }
    size_t localRowId = rowId - meta.row_id_offset;
    if (meta.uvws.size() <= localRowId * 3)
      meta.uvws.resize((localRowId + 1) * 3);
    for (size_t i = 0; i != 3; ++i) meta.uvws[localRowId * 3 + i] = uvw[i];
    double uvwFlipped[3] = {uvw[0], -uvw[1],
                            -uvw[2]};  // IDG uses a flipped coordinate system
    while (bs.get_degridder(data_desc_id)
               ->request_visibilities(rowId, timeIndex, antenna1, antenna2,
                                      uvwFlipped)) {
      ComputePredictionBuffer(data_desc_id, direction);
    }
  }
}

std::vector<DPBuffer> FacetPredict::Predict(size_t data_desc_id,
                                            size_t direction) {
  const size_t n_baselines = ant1_.size();
  const size_t n_timesteps = buffers_.size();
  const size_t n_terms = readers_.size();
  const size_t baseline_size = info_.nchan() * info_.ncorr();
  const size_t timestep_size = n_baselines * baseline_size;
  const size_t term_size = n_timesteps * timestep_size;

  std::vector<DPBuffer> result;
  std::vector<const double*> uvws;
  // term_data indexing: [term-1][timestep][baseline][channel][correlation].
  aocommon::UVector<std::complex<float>> term_data((n_terms - 1) * term_size);
  uvws.reserve(n_timesteps);

  double old_max_w = max_w_;
  const double max_baseline2 = max_baseline_ * max_baseline_;

  for (DPBuffer& buffer : buffers_) {
    uvws.push_back(input_.fetchUVW(buffer, buffer, timer_).data());

    const double* uvw = uvws.back();
    for (std::size_t bl = 0; bl < n_baselines; ++bl) {
      if (uvw[2] > max_w_) {
        double uvwr2 = uvw[0] * uvw[0] + uvw[1] * uvw[1] + uvw[2] * uvw[2];
        if (uvwr2 <= max_baseline2) {
          max_w_ = uvw[2];
        }
      }
      uvw += 3;
    }

    result.emplace_back(buffer);
    result.back().setData(casacore::Cube<casacore::Complex>(
        info_.ncorr(), info_.nchan(), n_baselines));
  }

  if (max_w_ > old_max_w) {
    // Increase max_w by 10% for preventing repeated IDG initialization calls.
    max_w_ *= 1.1;
    std::cout << "Increasing maximum w to " << max_w_ << '\n';
    StartIDG(false);
  }

  // IDG uses a flipped coordinate system
  static const double kUVWFactors[3] = {1.0, -1.0, -1.0};

  // Compute visibilities for all terms and timesteps.
  for (size_t term = 0; term < n_terms; ++term) {
    std::vector<std::complex<float>*> data_ptrs;
    data_ptrs.reserve(n_timesteps);

    if (term == 0) {  // data_ptrs point directly to the buffers.
      for (size_t t = 0; t < n_timesteps; ++t) {
        data_ptrs.push_back(result[t].getData().data());
      }
    } else {  // Use term_data as auxiliary buffers for the other terms.
      std::complex<float>* data = term_data.data() + (term - 1) * term_size;
      for (size_t t = 0; t < n_timesteps; ++t) {
        data_ptrs.push_back(data);
        data += timestep_size;
      }
    }

    idg::api::BufferSet& bs = *buffersets_[direction * n_terms + term];
    bs.get_bulk_degridder(data_desc_id)
        ->compute_visibilities(ant1_, ant2_, uvws, data_ptrs, kUVWFactors);
  }

  // Correct phase shift for the computed visibilities. Use the loop over nterms
  // as the inner loop, since it allows reusing the frequency factor.
  for (size_t t = 0; t < buffers_.size(); ++t) {
    const double* uvw = uvws[t];

    // Create pointers to the data for each term.
    std::vector<std::complex<float>*> term_data_ptrs;
    term_data_ptrs.reserve(n_terms);
    term_data_ptrs.push_back(result[t].getData().data());
    for (size_t term = 1; term < n_terms; ++term) {
      term_data_ptrs.push_back(term_data.data() + (term - 1) * term_size +
                               t * timestep_size);
    }

    for (size_t bl = 0; bl < n_baselines; ++bl) {
      // The phase shift angle is:
      // 2*pi * (u*dl + v*dm + w*dp) * (frequency / speed_of_light)
      // frequency_factor contains all parts except the frequency.
      const double frequency_factor = (uvw[0] * meta_data_[direction].dl +
                                       uvw[1] * meta_data_[direction].dm +
                                       uvw[2] * meta_data_[direction].dp) *
                                      (2.0 * M_PI / 299792458.0);

      for (auto& term_data_ptr : term_data_ptrs) {
        CorrectPhaseShift(term_data_ptr, frequency_factor);
        term_data_ptr += baseline_size;
      }

      uvw += 3;
    }
  }

  // Apply polynomial-term corrections and add all to values of 'term 0'
  // The "polynomial spectrum" definition is used, equal to the one e.g.
  // used by WSClean in component outputs (see
  // https://sourceforge.net/p/wsclean/wiki/ComponentList/ ) and in text
  // files when 'logarithmic SI' is false. The definition is:
  //   S(nu) = term0 + term1 (nu/refnu - 1) + term2 (nu/refnu - 1)^2 + ...
  for (size_t ch = 0; ch != info_.nchan(); ++ch) {
    double frequency = info_.chanFreqs()[ch];
    double freq_factor = frequency / ref_frequency_ - 1.0;
    double polynomial_factor = 1.0;
    const size_t ch_offset = ch * info_.ncorr();

    for (size_t term = 1; term != n_terms; ++term) {
      polynomial_factor *= freq_factor;
      const float polynomial_factor_float = polynomial_factor;
      for (size_t t = 0; t < n_timesteps; ++t) {
        std::complex<float>* values0 = result[t].getData().data() + ch_offset;
        const std::complex<float>* values = term_data.data() +
                                            (term - 1) * term_size +
                                            t * timestep_size + ch_offset;
        for (size_t bl = 0; bl < n_baselines; ++bl) {
          values0[0] += values[0] * polynomial_factor_float;
          values0[1] += values[1] * polynomial_factor_float;
          values0[2] += values[2] * polynomial_factor_float;
          values0[3] += values[3] * polynomial_factor_float;
          values0 += baseline_size;
          values += baseline_size;
        }
      }
    }
  }

  return result;
}

const std::vector<std::pair<double, double>>& FacetPredict::GetDirections()
    const {
  return directions_;
}

void FacetPredict::Flush(size_t data_desc_id) {
  for (size_t direction = 0; direction != directions_.size(); ++direction) {
    ComputePredictionBuffer(data_desc_id, direction);
  }
}

void FacetPredict::SetBufferSize(size_t n_timesteps) {
  buffer_size_ = n_timesteps;
}

void FacetPredict::CorrectPhaseShift(std::complex<float>* values,
                                     double frequency_factor) {
  for (size_t ch = 0; ch != info_.nchan(); ++ch) {
    double angle = frequency_factor * info_.chanFreqs()[ch];
    const std::complex<float> phasor(cos(angle), sin(angle));
    values[0] *= phasor;
    values[1] *= phasor;
    values[2] *= phasor;
    values[3] *= phasor;
    values += 4;
  }
}

void FacetPredict::ComputePredictionBuffer(size_t data_desc_id,
                                           size_t direction) {
  size_t n_terms = readers_.size();
  typedef std::vector<std::pair<size_t, std::complex<float>*>> rowidlist_t;
  std::vector<rowidlist_t> available_row_ids(n_terms);
  for (size_t term = 0; term != n_terms; ++term) {
    idg::api::BufferSet& bs = *buffersets_[direction * n_terms + term];
    available_row_ids[term] = bs.get_degridder(data_desc_id)->compute();
  }

  constexpr double two_pi_div_c = 2.0 * M_PI / 299792458.0;

  for (size_t i = 0; i != available_row_ids[0].size(); ++i) {
    size_t row = available_row_ids[0][i].first;
    size_t localRow = row - meta_data_[direction].row_id_offset;
    const double* uvw = &meta_data_[direction].uvws[localRow * 3];

    // The phase shift angle is:
    // 2*pi * (u*dl + v*dm + w*dp) * (frequency / speed_of_light)
    // frequency_factor contains all parts except the frequency.
    const double frequency_factor =
        (uvw[0] * meta_data_[direction].dl + uvw[1] * meta_data_[direction].dm +
         uvw[2] * meta_data_[direction].dp) *
        two_pi_div_c;

    // Correct the phase shift of the values for this facet
    for (size_t term = 0; term != n_terms; ++term) {
      std::complex<float>* values = available_row_ids[term][i].second;
      CorrectPhaseShift(values, frequency_factor);
    }

    // Apply polynomial-term corrections and add all to values of 'term 0'
    // The "polynomial spectrum" definition is used, equal to the one e.g.
    // used by WSClean in component outputs (see
    // https://sourceforge.net/p/wsclean/wiki/ComponentList/ ) and in text
    // files when 'logarithmic SI' is false. The definition is:
    //   S(nu) = term0 + term1 (nu/refnu - 1) + term2 (nu/refnu - 1)^2 + ...

    for (size_t ch = 0; ch != info_.nchan(); ++ch) {
      double frequency = info_.chanFreqs()[ch];
      double freq_factor = frequency / ref_frequency_ - 1.0;
      double polynomial_factor = 1.0;
      const std::size_t ch_offset = ch * 4;  // info_.ncorr() == 4
      for (size_t term = 1; term != n_terms; ++term) {
        polynomial_factor *= freq_factor;
        const float polynomial_factor_float = polynomial_factor;
        std::complex<float>* values0 =
            available_row_ids[0][i].second + ch_offset;
        const std::complex<float>* values =
            available_row_ids[term][i].second + ch_offset;
        values0[0] += values[0] * polynomial_factor_float;
        values0[1] += values[1] * polynomial_factor_float;
        values0[2] += values[2] * polynomial_factor_float;
        values0[3] += values[3] * polynomial_factor_float;
      }
    }
    predict_callback_(row, direction, data_desc_id,
                      available_row_ids[0][i].second);
  }
  for (size_t term = 0; term != n_terms; ++term) {
    idg::api::BufferSet& bs = *buffersets_[direction * n_terms + term];
    bs.get_degridder(data_desc_id)->finished_reading();
  }
  meta_data_[direction].is_initialized = false;
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
