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

#include <iostream>

FacetPredict::FacetPredict(const std::vector<std::string>& fits_model_files,
                           const std::string& ds9_regions_file,
                           PredictCallback&& callback)
    : predict_callback_(std::move(callback)), padding_(1.0), buffer_size_(0) {
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
                      idg::api::BufferSetType::degridding);
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
      computePredictionBuffer(data_desc_id, direction);
    }
  }
}

const std::vector<std::pair<double, double>>& FacetPredict::GetDirections()
    const {
  return directions_;
}

void FacetPredict::Flush(size_t data_desc_id) {
  for (size_t direction = 0; direction != directions_.size(); ++direction) {
    computePredictionBuffer(data_desc_id, direction);
  }
}

void FacetPredict::SetBufferSize(size_t n_timesteps) {
  buffer_size_ = n_timesteps;
}

void FacetPredict::computePredictionBuffer(size_t data_desc_id,
                                           size_t direction) {
  size_t n_terms = readers_.size();
  typedef std::vector<std::pair<size_t, std::complex<float>*>> rowidlist_t;
  std::vector<rowidlist_t> available_row_ids(n_terms);
  for (size_t term = 0; term != n_terms; ++term) {
    idg::api::BufferSet& bs = *buffersets_[direction * n_terms + term];
    available_row_ids[term] = bs.get_degridder(data_desc_id)->compute();
  }

  const double dl_fact = 2.0 * M_PI * meta_data_[direction].dl;
  const double dm_fact = 2.0 * M_PI * meta_data_[direction].dm;
  const double dp_fact = 2.0 * M_PI * meta_data_[direction].dp;
  for (size_t i = 0; i != available_row_ids[0].size(); ++i) {
    size_t row = available_row_ids[0][i].first;
    size_t localRow = row - meta_data_[direction].row_id_offset;
    const double* uvw = &meta_data_[direction].uvws[localRow * 3];

    // Correct the phase shift of the values for this facet
    for (size_t term = 0; term != n_terms; ++term) {
      std::complex<float>* values = available_row_ids[term][i].second;
      for (size_t ch = 0; ch != info_.nchan(); ++ch) {
        double angle = uvw[0] * dl_fact + uvw[1] * dm_fact + uvw[2] * dp_fact;
        angle *= info_.chanFreqs()[ch] / c();
        const std::complex<float> phasor(cos(angle), sin(angle));
        values[0] *= phasor;
        values[1] *= phasor;
        values[2] *= phasor;
        values[3] *= phasor;
        values += 4;
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

#else  // HAVE_IDG

namespace {

void notCompiled() {
  throw std::runtime_error(
      "Facet prediction is not available, because DP3 was not compiled with "
      "IDG support");
}

}  // namespace

FacetPredict::FacetPredict(const std::vector<std::string>&, const std::string&,
                           PredictCallback&&) {
  notCompiled();
}

void FacetPredict::SetMSInfo(std::vector<std::vector<double>>&& bands,
                             size_t nr_stations) {
  notCompiled();
}

bool FacetPredict::IsStarted() const {
  notCompiled();
  return false;
}

void FacetPredict::StartIDG(bool) { notCompiled(); }

void FacetPredict::RequestPredict(size_t, size_t, size_t, size_t, size_t,
                                  size_t, const double*) {
  notCompiled();
}

std::vector<std::pair<double, double>> FacetPredict::GetDirections(
    size_t) const {
  notCompiled();
  return std::vector<std::pair<double, double>>();
}

void FacetPredict::Flush() { notCompiled(); }

void FacetPredict::SetBufferSize(size_t) { notCompiled(); }

#endif  // HAVE_IDG
