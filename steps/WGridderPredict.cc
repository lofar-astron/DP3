// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "WGridderPredict.h"

#include <aocommon/banddata.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/logger.h>
#include <aocommon/threadpool.h>

#include <aocommon/uvector.h>
#include <aocommon/xt/utensor.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <ducc0/wgridder/wgridder.h>
#include <ducc0/fft/fftnd_impl.h>
#include <schaapcommon/facets/ds9facetfile.h>
#include <xtensor/xview.hpp>

#include "base/FlagCounter.h"
#include "common/Memory.h"

using namespace ducc0;

using aocommon::FitsReader;
using aocommon::FitsWriter;

using schaapcommon::facets::DS9FacetFile;
using schaapcommon::facets::Facet;
using schaapcommon::facets::FacetImage;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

using dp3::common::ParameterSet;

namespace {
size_t alignto4(size_t i) { return ((i + 3) / 4) * 4; }
}  // namespace

namespace dp3 {
namespace steps {

WGridderPredict::WGridderPredict(const ParameterSet& parset,
                                 const std::string& prefix)
    : WGridderPredict(parset, prefix,
                      GetReaders(parset.getStringVector(
                          prefix + "images", std::vector<std::string>())),
                      std::vector<Facet>(),
                      parset.getString(prefix + "regions", "")) {}

WGridderPredict::WGridderPredict(const ParameterSet& parset,
                                 const std::string& prefix,
                                 std::vector<FitsReader>&& readers,
                                 std::vector<Facet>&& facets,
                                 const std::string& ds9_regions_file)
    : name_(prefix),
      parset_(parset),
      reference_frequency_(readers.front().Frequency()),
      pixel_size_x_(readers.front().PixelSizeX()),
      pixel_size_y_(readers.front().PixelSizeY()),
      readers_(std::move(readers)),
      buffer_size_(parset.getBool(prefix + "buffersize", 0)),
      timer_(),
      directions_(),
      direction_labels_(),
      save_facets_(parset.getBool(prefix + "savefacets", false)),
      sum_facets_(parset.getBool(prefix + "sumfacets", false)) {
  if (facets.empty()) {
    facets = GetFacets(ds9_regions_file, readers_.front());
  }
  // TODO improve functionality: if no facets are defined, just process the
  // entire image

  const size_t n_facets = facets.size();
  const size_t n_terms = readers_.size();

  if (facets.empty()) {
    throw std::runtime_error("No facet definitions found in " +
                             ds9_regions_file);
  }

  size_t area = 0;
  directions_.reserve(n_facets);
  direction_labels_.reserve(n_facets);
  images_.reserve(n_facets);

  const std::vector<aocommon::Image> model_images = GetModelImages();
  std::vector<const float*> model_pointers;
  for (const aocommon::Image& model_image : model_images) {
    model_pointers.push_back(model_image.Data());
  }

  // Loop over facets to fill directions_, direction_labels_ and images_
  for (size_t facet_idx = 0; facet_idx < n_facets; facet_idx++) {
    const Facet& facet = facets[facet_idx];
    aocommon::Logger::Info << "Facet: '" << facet.DirectionLabel() << "'"
                           << " Ra,Dec: " << facet.RA() << "," << facet.Dec()
                           << " Vertices:";
    for (const schaapcommon::facets::PixelPosition& pixel : facet.GetPixels()) {
      aocommon::Logger::Info << " (" << pixel.x << "," << pixel.y << ")";
    }
    aocommon::Logger::Info << '\n';

    directions_.emplace_back(facet.RA(), facet.Dec());
    if (facet.DirectionLabel() != "") {
      direction_labels_.emplace_back(name_ + facet.DirectionLabel());
    } else {
      direction_labels_.emplace_back(name_ + "direction" +
                                     std::to_string(facet_idx));
    }

    const size_t padded_width = model_images.front().Width();
    const size_t padded_height = model_images.front().Height();

    images_.emplace_back(padded_width, padded_height, n_terms);
    FacetImage& image = images_.back();

    // The padding is 1.0, so the trimmed and untrimmed boxes are equal.
    const bool kTrimmed = true;
    image.SetFacet(facet, kTrimmed);
    image.CopyToFacet(model_pointers);
    area += image.Width() * image.Height();

    if (save_facets_) {
      FitsReader& reader = readers_.front();
      FitsWriter writer;
      writer.SetImageDimensions(image.Width(), image.Height(),
                                reader.PhaseCentreRA(), reader.PhaseCentreDec(),
                                pixel_size_x_, pixel_size_y_);
      const double l_shift = (int(reader.ImageWidth() / 2) -
                              (image.OffsetX() + int(image.Width() / 2))) *
                             pixel_size_x_;
      const double m_shift = (image.OffsetY() + int(image.Height() / 2) -
                              int(reader.ImageHeight() / 2)) *
                             pixel_size_y_;
      writer.SetPhaseCentreShift(l_shift, m_shift);
      writer.Write("facet" + std::to_string(meta_data_.size()) + ".fits",
                   image.Data(0));
    }
  }
  aocommon::Logger::Info << "Area covered: " << area / 1024 << " Kpixels^2\n";
}

std::vector<FitsReader> WGridderPredict::GetReaders(
    const std::vector<std::string>& fits_model_files) {
  if (fits_model_files.empty()) {
    throw std::runtime_error("No fits files specified for wgridder predict");
  }
  std::vector<FitsReader> readers;
  readers.reserve(fits_model_files.size());
  for (const std::string& file : fits_model_files) readers.emplace_back(file);

  return readers;
}

std::vector<aocommon::Image> WGridderPredict::GetModelImages() {
  const size_t unpadded_width = readers_.front().ImageWidth();
  const size_t unpadded_height = readers_.front().ImageHeight();

  // Round width and height upwards to the nearest multiple of 4
  // This Alignment is required by the faceting code in schaapcommon
  const size_t padded_width = alignto4(unpadded_width);
  const size_t padded_height = alignto4(unpadded_height);

  const double pixel_size_x = readers_.front().PixelSizeX();
  const double pixel_size_y = readers_.front().PixelSizeY();

  std::vector<aocommon::Image> models(readers_.size());

  for (size_t img = 0; img != readers_.size(); ++img) {
    if (readers_[img].ImageWidth() != unpadded_width ||
        readers_[img].ImageHeight() != unpadded_height)
      throw std::runtime_error("Image for spectral term " +
                               std::to_string(img) +
                               " has inconsistent dimensions");
    if (readers_[img].PixelSizeX() != pixel_size_x ||
        readers_[img].PixelSizeY() != pixel_size_y)
      throw std::runtime_error("Pixel size of spectral term " +
                               std::to_string(img) +
                               " is inconsistent with first spectral term");

    if ((padded_width == unpadded_width) &&
        (padded_height == unpadded_height)) {
      models[img] = aocommon::Image(unpadded_width, unpadded_height);
      readers_[img].Read(models[img].Data());
    } else {
      aocommon::Image image(unpadded_width, unpadded_height);
      readers_[img].Read(image.Data());
      // Untrim() makes a copy of the image.
      // We could make an in-place version in aocommon::Image.
      models[img] = image.Untrim(padded_width, padded_height);
    }
  }

  return models;
}

std::vector<Facet> WGridderPredict::GetFacets(
    const std::string& ds9_regions_file, const double ra, const double dec,
    const double pixel_size_x, const double pixel_size_y,
    const size_t unpadded_width, const size_t unpadded_height) {
  Facet::InitializationData facet_data(pixel_size_x, pixel_size_y,
                                       unpadded_width, unpadded_height);
  facet_data.phase_centre = schaapcommon::facets::Coord(ra, dec);
  facet_data.align = 4;
  facet_data.make_square = false;

  DS9FacetFile facet_file(ds9_regions_file);
  std::vector<Facet> facets = facet_file.Read(facet_data);
  aocommon::Logger::Info << "Read " << facets.size() << " facet definitions.\n";
  return facets;
}

std::vector<Facet> WGridderPredict::GetFacets(
    const std::string& ds9_regions_file, const FitsReader& reader) {
  const size_t padded_width = alignto4(reader.ImageWidth());
  const size_t padded_height = alignto4(reader.ImageHeight());
  return GetFacets(ds9_regions_file, reader.PhaseCentreRA(),
                   reader.PhaseCentreDec(), reader.PixelSizeX(),
                   reader.PixelSizeY(), padded_width, padded_height);
}

void WGridderPredict::updateInfo(const dp3::base::DPInfo& info_in) {
  if (info_in.ncorr() != 4) {
    throw std::invalid_argument(
        "WGridderPredict only supports 4 correlations.");
  }
  Step::updateInfo(info_in);
  std::map<std::string, dp3::base::Direction>& directions =
      GetWritableInfoOut().GetDirections();
  for (size_t i = 0; i < directions_.size(); i++) {
    directions.insert({direction_labels_[i], directions_[i]});
  }
  // Determine available size for buffering
  if (buffer_size_ == 0) {
    buffer_size_ = GetAllocatableBuffers(common::AvailableMemory());
    aocommon::Logger::Info << "Maximum number of timesteps to accumulate: "
                           << buffer_size_ << '\n';
  }
}

bool WGridderPredict::process(std::unique_ptr<base::DPBuffer> buffer) {
  buffers_.emplace_back(std::move(buffer));
  if (buffers_.size() >= buffer_size_) {
    flush();
  }
  return false;
}

void WGridderPredict::finish() {
  flush();
  getNextStep()->finish();
}

void WGridderPredict::flush() {
  if (buffers_.empty()) {
    return;
  }

  std::vector<base::DPBuffer::DataType*> destinations;
  destinations.reserve(buffers_.size());

  const std::array<size_t, 3> data_size = {
      getInfoOut().nbaselines(), getInfoOut().nchan(), getInfoOut().ncorr()};

  if (sum_facets_) {
    destinations.clear();
    for (size_t i = 0; i < buffers_.size(); i++) {
      destinations.push_back(&buffers_[i]->GetData());
      destinations.back()->resize(data_size);
      destinations.back()->fill(0.0);
    }
  }

  for (size_t dir = 0; dir < directions_.size(); ++dir) {
    if (!sum_facets_) {
      destinations.clear();
      for (size_t i = 0; i < buffers_.size(); i++) {
        buffers_[i]->AddData(direction_labels_[dir]);
        destinations.push_back(&buffers_[i]->GetData(direction_labels_[dir]));
        destinations.back()->resize(data_size);
        destinations.back()->fill(0.0);
      }
    }
    Predict(dir, destinations);
  }
  for (size_t i = 0; i < buffers_.size(); i++) {
    getNextStep()->process(std::move(buffers_[i]));
  }
  buffers_.clear();
}

void WGridderPredict::show(std::ostream& os) const {
  os << "WGridderPredict " << name_ << '\n';
  os << "  buffer size:    " << buffer_size_ << '\n';
  for (const FitsReader& reader : readers_) {
    os << "  Fits models " << reader.Filename() << '\n';
    os << "      RA:     " << reader.PhaseCentreRA() << '\n';
    os << "      Dec:    " << reader.PhaseCentreDec() << '\n';
  }
}

void WGridderPredict::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " WGridderPredict " << name_ << '\n';
}

void WGridderPredict::Predict(
    const size_t direction,
    std::vector<base::DPBuffer::DataType*>& destinations) {
  timer_.start();

  const size_t n_terms = readers_.size();
  const size_t n_timesteps = buffers_.size();
  const size_t n_baselines = getInfoOut().nbaselines();
  const size_t n_channels = getInfoOut().nchan();

  // Concatenate uvw data from the dpbuffers into one uvw buffer
  xt::xtensor<double, 3> uvw{
      xt::xtensor<double, 3>::shape_type{n_timesteps, n_baselines, 3}};
  for (size_t t = 0; t < n_timesteps; t++) {
    // wgridder expects an image that is transposed compared to fits images
    // that wsclean produces. To avoid an actual transpose of the images, the
    // u and v coordinates are swapped and the signs are changed.
    xt::view(uvw, t, xt::all(), 0) =
        -xt::view(buffers_[t]->GetUvw(), xt::all(), 1);
    xt::view(uvw, t, xt::all(), 1) =
        -xt::view(buffers_[t]->GetUvw(), xt::all(), 0);
    xt::view(uvw, t, xt::all(), 2) =
        xt::view(buffers_[t]->GetUvw(), xt::all(), 2);
  }

  // The visibilities predicted by wgridder are for a single polarization.
  // Unlike the visibilities in DPBuffer there is no correlations axis here.
  // The tranformation to full polarization data happens when the predicted
  // visibilities are added to the destination DPBuffer.
  xt::xtensor<std::complex<float>, 3> visibilities{
      xt::xtensor<std::complex<float>, 3>::shape_type{n_timesteps, n_baselines,
                                                      n_channels}};

  const double l_shift =
      (int(readers_[0].ImageWidth() / 2) -
       (images_[direction].OffsetX() + int(images_[direction].Width() / 2))) *
      pixel_size_x_;
  const double m_shift =
      (images_[direction].OffsetY() + int(images_[direction].Height() / 2) -
       int(readers_[0].ImageHeight() / 2)) *
      pixel_size_y_;

  const double* frequency_data = getInfoOut().chanFreqs().data();
  constexpr double sigma_min = 1.1;
  constexpr double sigma_max = 2.0;
  size_t width = images_[direction].Width();
  size_t height = images_[direction].Height();
  size_t nthreads = aocommon::ThreadPool::GetInstance().NThreads();
  double epsilon = 1.0e-4;
  size_t verbosity = 0;

  // View on uvw data in cmav type, as needed by wgridder
  cmav<double, 2> uvw_view(uvw.data(), {n_timesteps * n_baselines, 3});

  bool decreasing_freq =
      (n_channels > 1) && (frequency_data[1] < frequency_data[0]);

  // Create a view on the frequencies in ascending order
  // If the original frequencies are in descending order, a reversed view is
  // created, otherwise a normal view is created. The resulting view is in
  // ascending order, assuming the original frequencies were sorted either is
  // ascending or descending order
  auto frequencies_ascending(
      decreasing_freq
          ? cmav<double, 1>(frequency_data + n_channels - 1, {n_channels}, {-1})
          : cmav<double, 1>(frequency_data, {n_channels}));

  // weights is a required parameter to the wgridder, but it can be empty.
  cmav<float, 2> weights(nullptr, {0, 0});

  // mask is a required parameter to the wgridder, but it can be empty.
  cmav<std::uint8_t, 2> mask(nullptr, {0, 0});

  // Create a view on the visibilities in ascending frequency order
  // If the original freqencies are descending, the channel axis is reversed
  auto visibilities_ascending_frequency(
      decreasing_freq
          ? vmav<std::complex<float>, 2>(
                visibilities.data() + n_channels - 1,
                {n_timesteps * n_baselines, n_channels},
                {ptrdiff_t(n_channels), -1})
          : vmav<std::complex<float>, 2>(
                visibilities.data(), {n_timesteps * n_baselines, n_channels}));

  // Precompute polynomial frequency factors for all channels.
  std::vector<float> frequency_factors;
  frequency_factors.reserve(n_channels);
  for (double frequency : getInfoOut().chanFreqs()) {
    frequency_factors.push_back(frequency / reference_frequency_ - 1.0);
  }

  // For the first term (zero-th order) the polynomal factor is 1.0 for all
  // channels
  std::vector<float> polynomial_factors(n_channels, 1.0);

  // Compute visibilities for all spectral terms and timesteps.
  for (size_t term = 0; term < n_terms; ++term) {
    const float* image_data = images_[direction].Data(term);
    cmav<float, 2> model_image(image_data, {height, width});

    // Call to the wgridder. It computes the visibilities based on the model
    // image and uvw and frequency data. The order of x and y, and l and m are
    // reversed here, to avoid transposing the image as explained in the comment
    // where the uvw data is concatenated.
    dirty2ms<float, float>(uvw_view, frequencies_ascending, model_image,
                           weights, mask, pixel_size_y_, pixel_size_x_, epsilon,
                           true, nthreads, visibilities_ascending_frequency,
                           verbosity, true, false, sigma_min, sigma_max,
                           m_shift, l_shift);

    // Add the predicted visibilities to the destination buffers
    for (size_t t = 0; t < n_timesteps; ++t) {
      for (size_t bl = 0; bl < n_baselines; ++bl) {
        for (size_t ch = 0; ch < n_channels; ++ch) {
          // Add the computed Stokes I visibility to the XX and YY components of
          // the destination
          std::complex<float> weighted_visibility =
              visibilities(t, bl, ch) * polynomial_factors[ch];
          (*destinations[t])(bl, ch, 0) += weighted_visibility;
          (*destinations[t])(bl, ch, 3) += weighted_visibility;
        }
      }
    }
    // Compute the polynomial factors for the next term
    for (size_t ch = 0; ch < n_channels; ++ch) {
      polynomial_factors[ch] *= frequency_factors[ch];
    }
  }
  timer_.stop();
}

base::Direction WGridderPredict::GetFirstDirection() const {
  return {directions_.front().ra, directions_.front().dec};
}

void WGridderPredict::SetBufferSize(size_t n_timesteps) {
  buffer_size_ = n_timesteps;
}

size_t WGridderPredict::GetBufferSize() const { return buffer_size_; }

size_t WGridderPredict::GetAllocatableBuffers(size_t memory) {
  const size_t n_baselines = getInfoOut().nbaselines();
  const size_t n_channels = getInfoOut().nchan();
  const size_t n_correlations = getInfoOut().ncorr();

  // Crude estimate of the memory usage of one timestep
  // Curreently only accounts for the visibilities
  size_t mem_per_timestep =
      sizeof(std::complex<float>) * n_baselines * n_channels * n_correlations;

  // If the data is stored per direction an additional n_directions sets of
  // visibilities will be allocated
  if (!sum_facets_) {
    int n_directions = directions_.size();
    mem_per_timestep *= (n_directions + 1);
  }

  // At least one timestep and use at most 50% of available memory
  size_t n_buffers = std::max(memory / mem_per_timestep / 2, size_t(1));

  return n_buffers;
}

}  // namespace steps
}  // namespace dp3
