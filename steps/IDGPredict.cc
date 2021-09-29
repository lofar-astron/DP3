// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IDGPredict.h"
#include "InputStep.h"

#ifdef HAVE_IDG
#include "../base/IDGConfiguration.h"
#endif

#include "../base/ParsetAterms.h"
#include "../common/Memory.h"

#include <aocommon/uvector.h>
#include <aocommon/banddata.h>
#include <aocommon/fits/fitswriter.h>

#include <EveryBeam/aterms/atermconfig.h>

#include <schaapcommon/facets/ds9facetfile.h>

#include <iostream>

using aocommon::FitsReader;
using aocommon::FitsWriter;

using everybeam::aterms::ATermBase;
using everybeam::aterms::ATermConfig;

using schaapcommon::facets::DS9FacetFile;
using schaapcommon::facets::Facet;
using schaapcommon::facets::FacetImage;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

using dp3::common::ParameterSet;

namespace dp3 {
namespace steps {

IDGPredict::IDGPredict(InputStep& input, const ParameterSet& parset,
                       const string& prefix)
    : IDGPredict(input, parset, prefix,
                 GetReaders(parset.getStringVector(prefix + "images",
                                                   std::vector<string>())),
                 std::vector<Facet>(),
                 parset.getString(prefix + "regions", "")) {}

#ifdef HAVE_IDG
IDGPredict::IDGPredict(
    InputStep& input, const ParameterSet& parset, const string& prefix,
    std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<float>>>
        readers,
    std::vector<Facet>&& facets, const std::string& ds9_regions_file)
    : name_(prefix),
      parset_(parset),
      ref_frequency_(readers.first.front().Frequency()),
      pixel_size_x_(readers.first.front().PixelSizeX()),
      pixel_size_y_(readers.first.front().PixelSizeY()),
      readers_(std::move(readers.first)),
      buffer_size_(0),
      input_(input),
      ant1_(),
      ant2_(),
      timer_(),
      max_w_(0.0),
      max_baseline_(0.0),
      directions_(),
      save_facets_(parset.getBool(prefix + "savefacets", false)),
      aterms_(),
      aterm_values_() {
  const size_t full_width = readers_.front().ImageWidth();
  const size_t full_height = readers_.front().ImageHeight();
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
    for (const schaapcommon::facets::Pixel& pixel : facet.GetPixels()) {
      std::cout << " (" << pixel.x << "," << pixel.y << ")";
    }
    std::cout << '\n';

    directions_.emplace_back(facet.RA(), facet.Dec());
    images_.emplace_back(full_width, full_height, readers.second.size());
    FacetImage& image = images_.back();

    // The padding is 1.0, so the trimmed and untrimmed boxes are equal.
    const bool kTrimmed = true;
    image.SetFacet(facet, kTrimmed);
    image.CopyToFacet(readers.second);
    area += image.Width() * image.Height();
  }
  std::cout << "Area covered: " << area / 1024 << " Kpixels^2\n";
}
#endif  // HAVE_IDG

std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<float>>>
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

  std::vector<aocommon::UVector<float>> models(readers.size());
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
  const double kShiftL = 0.0;
  const double kShiftM = 0.0;
  const double kPadding = 1.0;
  const size_t kAlign = 4;
  const bool kMakeSquare = true;  // Only necessary for IDG.

  DS9FacetFile facet_file(ds9_regions_file);
  std::vector<Facet> facets = facet_file.Read();
  for (Facet& facet : facets) {
    facet.CalculatePixels(ra, dec, pixel_size_x, pixel_size_y, full_width,
                          full_height, kShiftL, kShiftM, kPadding, kAlign,
                          kMakeSquare);
  }
  std::cout << "Read " << facets.size() << " facet definitions.\n";
  return facets;
}

std::vector<Facet> IDGPredict::GetFacets(const std::string& ds9_regions_file,
                                         const FitsReader& reader) {
  return GetFacets(ds9_regions_file, reader.PhaseCentreRA(),
                   reader.PhaseCentreDec(), reader.PixelSizeX(),
                   reader.PixelSizeY(), reader.ImageWidth(),
                   reader.ImageHeight());
}

#ifdef HAVE_IDG

void IDGPredict::updateInfo(const dp3::base::DPInfo& info) {
  if (info.ncorr() != 4) {
    throw std::invalid_argument("IDGPredict only supports 4 correlations.");
  }

  Step::updateInfo(info);

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
    buffer_size_ = GetAllocatableBuffers(common::AvailableMemory());
  }

  StartIDG();

  if (!parset_.getStringVector(name_ + "aterms", std::vector<string>())
           .empty()) {
    if (info.nantenna() != input_.getInfo().nantenna()) {
      throw std::runtime_error(
          "Number of antennas not matching the number of antennas in the "
          "original MS. This is as yet unsupported for beam calculations.");
    }

    InitializeATerms();
  }
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
  os << " IDGPredict " << name_ << '\n';
}

bool IDGPredict::IsStarted() const { return !buffersets_.empty(); }

void IDGPredict::StartIDG() {
  buffersets_.clear();
  idg::api::Type proxy_type = idg::api::Type::CPU_OPTIMIZED;
  int buffersize = 0;
  size_t n_terms = readers_.size();

  idg::api::options_type options;

  IdgConfiguration::Read(proxy_type, buffersize, options);
  std::vector<double> image_data;
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
    const double dm = (img.OffsetY() + int(img.Height() / 2) -
                       int(reader.ImageHeight() / 2)) *
                      pixel_size_y_;
    const double dp = sqrt(1.0 - dl * dl - dm * dm) - 1.0;
    std::cout << "Initializing gridder " << buffersets_.size() << " ("
              << img.Width() << " x " << img.Height() << ", +" << img.OffsetX()
              << "," << img.OffsetY() << ", dl=" << dl * 180.0 / M_PI
              << " deg, dm=" << dm * 180.0 / M_PI << " deg)\n"
              << '\n';

    size_t subgrid_size = 0;
    for (size_t term = 0; term != n_terms; ++term) {
      // Reserve space for 4 polarizations but only use the first polarization.
      // TODO Copy data for the other polarizations, too.
      image_data.resize(img.Width() * img.Height() * 4, 0.0);
      std::copy_n(img.Data(term), img.Width() * img.Height(),
                  image_data.data());

      buffersets_.emplace_back(idg::api::BufferSet::create(proxy_type));
      idg::api::BufferSet& bs = *buffersets_.back();
      options["padded_size"] = size_t(1.2 * img.Width());
      // options["max_threads"] = int(1);
      bs.init(img.Width(), pixel_size_x_, max_w_ + 1.0, dl, dm, dp, options);
      bs.set_image(image_data.data());
      bs.init_buffers(0, {info().chanFreqs()}, info().antennaUsed().size(),
                      max_baseline_, options,
                      idg::api::BufferSetType::kBulkDegridding);

      // GetSubgridCount assumes that the subgrid size is equal for all terms.
      if (term == 0) {
        subgrid_size = bs.get_subgridsize();
      } else if (bs.get_subgridsize() != subgrid_size) {
        throw std::runtime_error("IDG subgrid sizes do not match");
      }
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

void IDGPredict::InitializeATerms() {
  const size_t n_terms = readers_.size();
  const size_t n_antennas = getInfo().nantenna();

  everybeam::ATermSettings settings;
  // https://wsclean.readthedocs.io/en/latest/a_term_correction.html?highlight=kernel#kernel-size
  // describes why 16 is a reasonable default value.
  settings.max_support = parset_.getInt("atermkernelsize", 16);

  aterms_.clear();
  aterms_.reserve(directions_.size());
  aterm_values_.clear();
  aterm_values_.reserve(directions_.size());
  for (size_t direction = 0; direction < directions_.size(); ++direction) {
    const idg::api::BufferSet& bs = *buffersets_[direction * n_terms];

    everybeam::coords::CoordinateSystem cs;
    // IDG uses a flipped coordinate cs which is moved by half a pixel:
    cs.dl = -bs.get_subgrid_pixelsize();
    cs.dm = -bs.get_subgrid_pixelsize();
    cs.phase_centre_dl = meta_data_[direction].dl - 0.5 * cs.dl;
    cs.phase_centre_dm = meta_data_[direction].dm + 0.5 * cs.dm;
    cs.width = bs.get_subgridsize();
    cs.height = cs.width;
    cs.ra = directions_[direction].first;
    cs.dec = directions_[direction].second;

    auto* aterm = new ATermConfig(n_antennas, cs, settings);
    aterm->Read(input_.table(), base::ParsetATerms(parset_, name_));
    aterms_.emplace_back(aterm);
    aterm_values_.emplace_back(n_antennas * GetSubgridCount(direction));
  }
}

std::vector<DPBuffer> IDGPredict::Predict(const size_t direction) {
  timer_.start();
  const size_t n_terms = readers_.size();
  const size_t term_size =
      buffers_.size() * info().nbaselines() * info().nchan() * info().ncorr();

  // term_data indexing: [term-1][timestep][baseline][channel][correlation].
  aocommon::UVector<std::complex<float>> term_data((n_terms - 1) * term_size);

  std::vector<const double*> uvws = InitializeUVWs();

  std::vector<DPBuffer> result =
      ComputeVisibilities(direction, uvws, term_data.data());

  CorrectVisibilities(result, term_data.data());

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
    std::complex<float>* term_data) const {
  const size_t n_timesteps = buffers_.size();
  const size_t n_terms = readers_.size();
  const size_t baseline_size = getInfo().nchan() * getInfo().ncorr();
  const size_t timestep_size = getInfo().nbaselines() * baseline_size;
  const size_t term_size = n_timesteps * timestep_size;

  std::vector<DPBuffer> result;
  result.reserve(buffers_.size());
  for (const DPBuffer& buffer : buffers_) {
    result.emplace_back(buffer);
    result.back().setData(casacore::Cube<casacore::Complex>(
        getInfo().ncorr(), getInfo().nchan(), getInfo().nbaselines()));
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

    // Since we initialize the BufferSet with a single band, the index is 0.
    constexpr int kDegridderIndex = 0;
    idg::api::BufferSet& bs = *buffersets_[direction * n_terms + term];
    const idg::api::BulkDegridder* degridder =
        bs.get_bulk_degridder(kDegridderIndex);

    if (aterms_.empty()) {
      degridder->compute_visibilities(ant1_, ant2_, uvws, data_ptrs,
                                      kUVWFactors);
    } else {
      // ATerm computations are as yet only meaningful for the leading term.
      assert(term == 0);

      // Note: IDGPredict currently does not support using different aterms
      // for different time slots.
      const aocommon::UVector<std::complex<float>> aterm_values =
          GetAtermValues(direction);

      degridder->compute_visibilities(ant1_, ant2_, uvws, data_ptrs,
                                      kUVWFactors, aterm_values.data());
    }
  }

  return result;
}

aocommon::UVector<std::complex<float>> IDGPredict::GetAtermValues(
    size_t direction) const {
  assert(direction < aterms_.size());
  ATermBase& aterm = *aterms_[direction];
  aocommon::UVector<std::complex<float>>& aterm_values =
      aterm_values_[direction];

  const double time_centroid =
      0.5 * (buffers_.front().getTime() + buffers_.back().getTime());

  // TODO: field_id should be made generic in case DP3 will be used for
  // other telescopes
  constexpr size_t field_id = 0;

  // Computes aterm values for all stations. On the first invocation,
  // aterm.Calculate should always return true and fill aterm_values.
  // On successive invocations, it may return false, which means it does not
  // touch aterm_values and IDGPredict should reuse the old values in the cache.
  aterm.Calculate(aterm_values.data(), time_centroid, getInfo().refFreq(),
                  field_id, nullptr);

  const std::vector<int>& ant_used = getInfo().antennaUsed();
  if (getInfo().nantenna() == ant_used.size()) {
    return aterm_values;
  } else {
    // Copy the relevant data into a buffer for the used antennas.
    const size_t subgrid_count = GetSubgridCount(direction);
    aocommon::UVector<std::complex<float>> aterm_used(ant_used.size() *
                                                      subgrid_count);
    auto aterm_used_it = aterm_used.begin();
    for (int ant : ant_used) {
      std::copy_n(aterm_values.begin() + ant * subgrid_count, subgrid_count,
                  aterm_used_it);
      aterm_used_it += subgrid_count;
    }
    return aterm_used;
  }
}

void IDGPredict::CorrectVisibilities(std::vector<DPBuffer>& result,
                                     const std::complex<float>* term_data) {
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
  std::vector<float> freq_factors;
  if (n_terms > 1) {
    freq_factors.reserve(info().nchan());
    for (double freq : info().chanFreqs()) {
      freq_factors.push_back(freq / ref_frequency_ - 1.0);
    }
  }

  for (size_t t = 0; t < n_timesteps; ++t) {
    std::complex<float>* result_data = result[t].getData().data();
    for (size_t bl = 0; bl < info().nbaselines(); ++bl) {
      for (size_t ch = 0; ch < info().nchan(); ++ch) {
        float polynomial_factor = 1.0;
        for (auto& term_data_ptr : term_data_ptrs) {
          polynomial_factor *= freq_factors[ch];
          result_data[0] += term_data_ptr[0] * polynomial_factor;
          result_data[1] += term_data_ptr[1] * polynomial_factor;
          result_data[2] += term_data_ptr[2] * polynomial_factor;
          result_data[3] += term_data_ptr[3] * polynomial_factor;
          term_data_ptr += 4;
        }
        result_data += 4;
      }
    }
  }
}

base::Direction IDGPredict::GetFirstDirection() const {
  // We can take the front element of an IDG step since it only contains 1.
  return {directions_.front().first, directions_.front().second};
}

void IDGPredict::SetBufferSize(size_t n_timesteps) {
  buffer_size_ = n_timesteps;
}

size_t IDGPredict::GetBufferSize() const { return buffer_size_; }

size_t IDGPredict::GetAllocatableBuffers(size_t memory) {
  size_t n_terms = readers_.size();

  size_t max_channels = info().chanFreqs().size();
  uint64_t memPerTimestep = idg::api::BufferSet::get_memory_per_timestep(
      info().antennaUsed().size(), max_channels);
  memPerTimestep *= 2;  // IDG uses two internal buffers
  memPerTimestep *= 4;  // DP3 can store up to 4 times the MS size in memory.
  // Allow the directions together to use 1/4th of the available memory for
  // the vis buffers. We divide by 3, because 3 copies are being kept in
  // memory.
  size_t allocatableTimesteps =
      memory / images_.size() / n_terms / memPerTimestep / 3;
  // TODO once a-terms are supported, this should include the size required
  // for the a-terms.
  std::cout << "Allocatable timesteps per direction: " << allocatableTimesteps
            << '\n';

  size_t buffersize = std::max(allocatableTimesteps, size_t(1));
  return buffersize;
}

size_t IDGPredict::GetSubgridCount(size_t direction) const {
  const size_t n_terms = readers_.size();
  assert(direction * n_terms < buffersets_.size());

  // The subgrid size should be equal for all terms, so use the first BufferSet
  // for the given direction.
  const idg::api::BufferSet& bs = *buffersets_[direction * n_terms];
  const size_t subgrid_size = bs.get_subgridsize();
  return subgrid_size * subgrid_size * getInfo().ncorr();
}

#else  // HAVE_IDG

namespace {

void notCompiled() {
  throw std::runtime_error(
      "IDG prediction is not available, because DP3 was not compiled with "
      "IDG support");
}

}  // namespace

IDGPredict::IDGPredict(
    InputStep& input, const ParameterSet& parset, const string& prefix,
    std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<float>>>
        readers,
    std::vector<Facet>&& facets, const std::string& ds9_regions_file) {
  // Do nothing / create dummy object.
}

void IDGPredict::updateInfo(const dp3::base::DPInfo&) { notCompiled(); }

bool IDGPredict::IsStarted() const {
  notCompiled();
  return false;
}

base::Direction IDGPredict::GetFirstDirection() const {
  notCompiled();
  return base::Direction();
}

void IDGPredict::flush() { notCompiled(); }

void IDGPredict::SetBufferSize(size_t) { notCompiled(); }

size_t IDGPredict::GetBufferSize() const {
  notCompiled();
  return 0;
}

bool IDGPredict::process(const DPBuffer&) {
  notCompiled();
  return false;
}

void IDGPredict::finish() { notCompiled(); }

void IDGPredict::show(std::ostream&) const { notCompiled(); }

void IDGPredict::showTimings(std::ostream&, double) const { notCompiled(); }

#endif  // HAVE_IDG

}  // namespace steps
}  // namespace dp3
