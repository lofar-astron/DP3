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

#include "IDGPredict.h"

#ifdef HAVE_IDG

#include "DS9FacetFile.h"
#include "IDGConfiguration.h"
#include "../Common/Memory.h"
#include "../DPPP/DPInput.h"

#include <aocommon/uvector.h>
#include <aocommon/banddata.h>
#include <aocommon/fits/fitswriter.h>

#include <iostream>

using aocommon::FitsWriter;

namespace DP3 {
namespace DPPP {

IDGPredict::IDGPredict(
    DPInput& input, const ParameterSet& parset, const string& prefix,
    std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<double>>>
        readers,
    std::vector<Facet>& facets, const std::string& ds9_regions_file)
    : name_(prefix),
      readers_(readers.first),
      buffer_size_(0),
      input_(input),
      ant1_(),
      ant2_(),
      timer_(),
      save_facets_(parset.getBool(prefix + "savefacets", false)) {
  const size_t full_width = readers_.front().ImageWidth();
  const size_t full_height = readers_.front().ImageHeight();
  ref_frequency_ = readers_.front().Frequency();
  pixel_size_x_ = readers_.front().PixelSizeX();
  pixel_size_y_ = readers_.front().PixelSizeY();
  if (facets.empty()) {
    facets = GetFacets(ds9_regions_file, readers_.front());
  }
  // TODO improve functionality: if no facets are defined, just process the
  // entire image
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
    image.CopyFacetPart(facet, readers.second, full_width, full_height,
                        kPadding, kMakeSquare);
    area += image.Width() * image.Height();
  }
  std::cout << "Area covered: " << area / 1024 << " Kpixels^2\n";
}

IDGPredict::IDGPredict(DPInput& input, const ParameterSet& parset,
                       const string& prefix)
    : IDGPredict(input, parset, prefix,
                 GetReaders(parset.getStringVector(prefix + "images",
                                                   std::vector<string>())),
                 *(new std::vector<Facet>()),
                 parset.getString(prefix + "regions", "")) {}

std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<double>>>
IDGPredict::GetReaders(const std::vector<std::string>& fits_model_files) {
  if (fits_model_files.empty()) {
    throw std::runtime_error("No fits files specified for IDG predict");
  }
  std::vector<FitsReader> readers;
  readers.reserve(fits_model_files.size());
  for (const std::string& file : fits_model_files) readers.emplace_back(file);

  const size_t full_width = readers.front().ImageWidth();
  const size_t full_height = readers.front().ImageHeight();
  const double pixel_size_x = readers.front().PixelSizeX();
  const double pixel_size_y = readers.front().PixelSizeY();

  std::vector<aocommon::UVector<double>> models(readers.size());
  for (size_t img = 0; img != readers.size(); ++img) {
    if (readers[img].ImageWidth() != full_width ||
        readers[img].ImageHeight() != full_height)
      throw std::runtime_error("Image for spectral term " +
                               std::to_string(img) +
                               " has inconsistent dimensions");
    if (readers[img].PixelSizeX() != pixel_size_x ||
        readers[img].PixelSizeY() != pixel_size_y)
      throw std::runtime_error("Pixel size of spectral term " +
                               std::to_string(img) +
                               " is inconsistent with first spectral term");
    models[img].resize(full_width * full_height);
    readers[img].Read(models[img].data());
  }

  return std::make_pair(readers, models);
}

std::vector<Facet> IDGPredict::GetFacets(const std::string& ds9_regions_file,
                                         const double ra, const double dec,
                                         const double pixel_size_x,
                                         const double pixel_size_y,
                                         const size_t full_width,
                                         const size_t full_height) {
  DS9FacetFile facet_file(ds9_regions_file);
  std::vector<Facet> facets_out = facet_file.Read(
      ra, dec, pixel_size_x, pixel_size_y, full_width, full_height);
  std::cout << "Read " << facets_out.size() << " facet definitions.\n";
  return facets_out;
}

std::vector<Facet> IDGPredict::GetFacets(const std::string& ds9_regions_file,
                                         const FitsReader& reader) {
  return GetFacets(ds9_regions_file, reader.PhaseCentreRA(),
                   reader.PhaseCentreDec(), reader.PixelSizeX(),
                   reader.PixelSizeY(), reader.ImageWidth(),
                   reader.ImageHeight());
}

void IDGPredict::updateInfo(const DP3::DPPP::DPInfo& info) {
  if (info.ncorr() != 4) {
    throw std::invalid_argument("IDGPredict only supports 4 correlations.");
  }

  DPStep::updateInfo(info);

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

  // Determine available size for buffering
  if (buffer_size_ == 0) {
    buffer_size_ = GetAllocatableBuffers(AvailableMemory());
  }

  StartIDG();
}

bool IDGPredict::process(const DPBuffer& buffer) {
  buffers_.emplace_back(buffer);
  input_.fetchUVW(buffer, buffers_.back(), timer_);

  if (buffers_.size() == buffer_size_) {
    flush();
  }

  return false;
}

void IDGPredict::finish() {
  flush();
  getNextStep()->finish();
}

void IDGPredict::flush() {
  if (buffers_.empty()) {
    return;
  }

  std::vector<DPBuffer> front_prediction = Predict(0);

  for (size_t dir = 1; dir < directions_.size(); ++dir) {
    auto predictions = Predict(dir);
    for (size_t i = 0; i < predictions.size(); ++i) {
      front_prediction[i].getData() += predictions[i].getData();
    }
  }

  for (DPBuffer& prediction : front_prediction) {
    getNextStep()->process(prediction);
  }

  buffers_.clear();
}

void IDGPredict::show(std::ostream& os) const {
  os << "IDGPredict " << name_ << '\n';
  os << "  buffer size:    " << buffer_size_ << '\n';
  for (const FitsReader& reader : readers_) {
    os << "  Fits models " << reader.Filename() << '\n';
    os << "      RA:     " << reader.PhaseCentreRA() << '\n';
    os << "      Dec:    " << reader.PhaseCentreDec() << '\n';
  }
}

void IDGPredict::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " IDGPredict " << name_ << std::endl;
}

bool IDGPredict::IsStarted() const { return !buffersets_.empty(); }

void IDGPredict::StartIDG() {
  buffersets_.clear();
  idg::api::Type proxy_type = idg::api::Type::CPU_OPTIMIZED;
  int buffersize = 0;
  size_t n_terms = readers_.size();

  idg::api::options_type options;

  IdgConfiguration::Read(proxy_type, buffersize, options);
  std::vector<aocommon::UVector<double>> data(n_terms);
  meta_data_.clear();

  options["disable_wtiling"] = true;

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
    // TODO revert to
    // const double dm = (int(reader.ImageHeight() / 2) -
    //                   (img.OffsetY() + int(img.Height() / 2))) *
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
      bs.init_buffers(0, {info().chanFreqs()}, info().antennaUsed().size(),
                      max_baseline_, options,
                      idg::api::BufferSetType::kBulkDegridding);
    }

    if (save_facets_) {
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

std::vector<DPBuffer> IDGPredict::Predict(const size_t direction) {
  timer_.start();
  const size_t n_terms = readers_.size();
  const size_t term_size =
      buffers_.size() * info().nbaselines() * info().nchan() * info().ncorr();

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

  timer_.stop();
  return result;
}

std::vector<const double*> IDGPredict::InitializeUVWs() {
  std::vector<const double*> uvws;
  uvws.reserve(buffers_.size());

  double old_max_w = max_w_;
  const double max_baseline2 = max_baseline_ * max_baseline_;

  for (DPBuffer& buffer : buffers_) {
    uvws.push_back(input_.fetchUVW(buffer, buffer, timer_).data());

    const double* uvw = uvws.back();
    for (std::size_t bl = 0; bl < info().nbaselines(); ++bl) {
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
    save_facets_ = false;
    StartIDG();
  }

  return uvws;
}

std::vector<DPBuffer> IDGPredict::ComputeVisibilities(
    const size_t direction, const std::vector<const double*>& uvws,
    std::complex<float>* term_data) {
  const size_t n_timesteps = buffers_.size();
  const size_t n_terms = readers_.size();
  const size_t baseline_size = info().nchan() * info().ncorr();
  const size_t timestep_size = info().nbaselines() * baseline_size;
  const size_t term_size = n_timesteps * timestep_size;

  std::vector<DPBuffer> result;
  result.reserve(buffers_.size());
  for (const DPBuffer& buffer : buffers_) {
    result.emplace_back(buffer);
    result.back().setData(casacore::Cube<casacore::Complex>(
        info().ncorr(), info().nchan(), info().nbaselines()));
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
    std::cout << "Starting bulk degridding" << std::endl;
    bs.get_bulk_degridder(kDegridderIndex)
        ->compute_visibilities(ant1_, ant2_, uvws, data_ptrs, kUVWFactors);
    std::cout << "Finished bulk degridding" << std::endl;
  }

  return result;
}

double IDGPredict::ComputePhaseShiftFactor(const double* uvw,
                                           size_t direction) const {
  // The phase shift angle is:
  // 2*pi * (u*dl - v*dm - w*dp) * (frequency / speed_of_light)
  // TODO revert The v and w terms are negated since IDG uses a flipped
  // coordinate system. The computed phase shift factor contains all parts
  // except the frequency.
  constexpr double two_pi_div_c = 2.0 * M_PI / aocommon::c();
  return (uvw[0] * meta_data_[direction].dl +
          uvw[1] * meta_data_[direction].dm +
          uvw[2] * meta_data_[direction].dp) *
         two_pi_div_c;
}

void IDGPredict::CorrectVisibilities(const std::vector<const double*>& uvws,
                                     std::vector<DPBuffer>& result,
                                     const std::complex<float>* term_data,
                                     size_t direction) {
  const size_t n_timesteps = buffers_.size();
  const size_t n_terms = readers_.size();
  const size_t baseline_size = info().nchan() * info().ncorr();
  const size_t timestep_size = info().nbaselines() * baseline_size;
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
    freq_factors.reserve(info().nchan());
    for (double freq : info().chanFreqs()) {
      freq_factors.push_back(freq / ref_frequency_ - 1.0);
    }
  }

  for (size_t t = 0; t < n_timesteps; ++t) {
    const double* uvw = uvws[t];
    std::complex<float>* result_data = result[t].getData().data();
    for (size_t bl = 0; bl < info().nbaselines(); ++bl) {
      const double phase_factor = ComputePhaseShiftFactor(uvw, direction);
      uvw += 3;

      for (size_t ch = 0; ch < info().nchan(); ++ch) {
        const double angle = phase_factor * info().chanFreqs()[ch];
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

const std::vector<std::pair<double, double>>& IDGPredict::GetDirections()
    const {
  return directions_;
}

void IDGPredict::SetBufferSize(size_t n_timesteps) {
  buffer_size_ = n_timesteps;
}

size_t IDGPredict::GetAllocatableBuffers(size_t memory) {
  size_t n_terms = readers_.size();

  size_t max_channels = info().chanFreqs().size();
  uint64_t memPerTimestep = idg::api::BufferSet::get_memory_per_timestep(
      info().antennaUsed().size(), max_channels);
  memPerTimestep *= 2;  // IDG uses two internal buffers
  memPerTimestep *= 4;  // DP3 can store up to 4 times the MS size in memory.
  // Allow the directions together to use 1/4th of the available memory for
  // the vis buffers. We divide by 3, because 3 copies are being kept in memory.
  size_t allocatableTimesteps =
      memory / images_.size() / n_terms / memPerTimestep / 3;
  // TODO once a-terms are supported, this should include the size required
  // for the a-terms.
  std::cout << "Allocatable timesteps per direction: " << allocatableTimesteps
            << '\n';

  size_t buffersize = std::max(allocatableTimesteps, size_t(1));
  return buffersize;
}

}  // namespace DPPP
}  // namespace DP3

#else  // HAVE_IDG

namespace {

void notCompiled() {
  throw std::runtime_error(
      "IDG prediction is not available, because DP3 was not compiled with "
      "IDG support");
}

}  // namespace

namespace DP3 {
namespace DPPP {

IDGPredict::IDGPredict(const std::vector<std::string>&, const std::string&,
                       PredictCallback&&) {
  notCompiled();
}

void IDGPredict::updateInfo(const DP3::DPPP::DPInfo&) { notCompiled(); }

bool IDGPredict::IsStarted() const {
  notCompiled();
  return false;
}

void IDGPredict::StartIDG() { notCompiled(); }

void IDGPredict::RequestPredict(size_t, size_t, size_t, size_t, size_t, size_t,
                                const double*) {
  notCompiled();
}

const std::vector<std::pair<double, double>>& IDGPredict::GetDirections()
    const {
  notCompiled();
  return std::vector<std::pair<double, double>>();
}

void IDGPredict::Flush(size_t) { notCompiled(); }

void IDGPredict::SetBufferSize(size_t) { notCompiled(); }

}  // namespace DPPP
}  // namespace DP3

#endif  // HAVE_IDG