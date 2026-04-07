// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_COMMON_DYNSPEC_WRITER_H_
#define DP3_COMMON_DYNSPEC_WRITER_H_

#include <aocommon/fits/fitsbase.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/uvector.h>

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <fitsio.h>

namespace dp3::common {

/// Dynamic Spectrum FITS writer. A specialized version of aocommon::FitsWriter
/// that can only write FITS files with axes: TIME, FREQ, STOKES.
class DynSpecFitsWriter : private aocommon::FitsBase {
 public:
  DynSpecFitsWriter()
      : n_times_(0),
        n_channels_(0),
        n_stokes_(4),
        object_ra_(0.0),
        object_dec_(0.0),
        obs_date_(0.0),
        time_resolution_(0.0),
        frequency_(0.0),
        frequency_resolution_(0.0),
        telescope_name_(),
        observer_(),
        object_name_(),
        origin_("DP3"),
        origin_comment_("DP3/DynSpec - dynamic spectrum extracter"),
        obs_id_() {}

  ~DynSpecFitsWriter() {
    if (file_pointer != nullptr) {
      int status = 0;
      fits_close_file(file_pointer, &status);
    }
  }

  template <typename NumType>
  void Write(const std::string& filename, const NumType* image) {
    int status = 0;
    if (file_pointer != nullptr) {
      fits_close_file(file_pointer, &status);
      file_pointer = nullptr;
    }

    fits_create_file(&file_pointer, (std::string("!") + filename).c_str(),
                     &status);
    checkStatus(status, filename);

    writeHeaders(file_pointer, filename);

    long first_pixel[3] = {1, 1, 1};
    writeImage(file_pointer, filename, image, first_pixel);

    fits_close_file(file_pointer, &status);
    checkStatus(status, filename);
    file_pointer = nullptr;
  }

  void InitializeTimeAxis(size_t n_times, double time_resolution,
                          double obs_start_time) {
    n_times_ = n_times;
    time_resolution_ = time_resolution;
    obs_date_ = obs_start_time;
  }

  void InitializeFrequencyAxis(size_t n_channels, double frequency_resolution,
                               double reference_frequency) {
    n_channels_ = n_channels;
    frequency_resolution_ = frequency_resolution;
    frequency_ = reference_frequency;
  }

  void InitializeStokesAxis(size_t n_stokes_parameters) {
    n_stokes_ = n_stokes_parameters;
  }

  void SetObjectName(const std::string& object_name) {
    object_name_ = object_name;
  }

  void SetObjectCoordinates(double ra, double dec) {
    object_ra_ = ra;
    object_dec_ = dec;
  }

  void SetOrigin(const std::string& origin, const std::string& comment) {
    origin_ = origin;
    origin_comment_ = comment;
  }

  void SetObsId(const std::string& obs_id) { obs_id_ = obs_id; }

 private:
  void writeHeaders(fitsfile*& fptr, const std::string& filename) const {
    // append image HDU
    int bitPixInt = FLOAT_IMG;
    std::vector<long> naxes(3);
    naxes[0] = n_times_;
    naxes[1] = n_channels_;
    naxes[2] = n_stokes_;

    int status = 0;
    fits_create_img(fptr, bitPixInt, naxes.size(), naxes.data(), &status);
    checkStatus(status, filename);
    const double zero = 0.0;
    const double one = 1.0;
    fits_write_key(fptr, TDOUBLE, "BSCALE", (void*)&one, "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "BZERO", (void*)&zero, "", &status);
    checkStatus(status, filename);

    fits_write_key(fptr, TSTRING, "BUNIT", (void*)"JY", "Units are in Jansky",
                   &status);
    checkStatus(status, filename);

    if (!telescope_name_.empty()) {
      fits_write_key(fptr, TSTRING, "TELESCOP", (void*)telescope_name_.c_str(),
                     "", &status);
      checkStatus(status, filename);
    }
    if (!object_name_.empty()) {
      fits_write_key(fptr, TSTRING, "OBJECT", (void*)object_name_.c_str(), "",
                     &status);
      checkStatus(status, filename);
    }
    if (!obs_id_.empty()) {
      fits_write_key(fptr, TSTRING, "OBSID", (void*)obs_id_.c_str(),
                     "Observation ID", &status);
      checkStatus(status, filename);
    }
    fits_write_key(fptr, TSTRING, "ORIGIN", (void*)origin_.c_str(),
                   origin_comment_.c_str(), &status);
    checkStatus(status, filename);

    // First axis is TIME
    fits_write_key(fptr, TSTRING, "CTYPE1", (void*)"TIME", "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", (void*)&one, "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", (void*)&obs_date_, "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CDELT1", (void*)&time_resolution_, "",
                   &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TSTRING, "CUNIT1", (void*)"s", "", &status);
    checkStatus(status, filename);

    // Second axis is FREQUENCY
    fits_write_key(fptr, TSTRING, "CTYPE2", (void*)"FREQ", "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", (void*)&one, "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", (void*)&frequency_, "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CDELT2", (void*)&frequency_resolution_, "",
                   &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TSTRING, "CUNIT2", (void*)"Hz", "", &status);
    checkStatus(status, filename);

    // Third axis is POLARIZATION: I, Q, U, V
    fits_write_key(fptr, TSTRING, "CTYPE3", (void*)"STOKES", "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CRPIX3", (void*)&one, "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CRVAL3", (void*)&one, "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "CDELT3", (void*)&one, "", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TSTRING, "CUNIT3", (void*)"", "", &status);
    checkStatus(status, filename);

    // Store the observation date
    int year, month, day, hour, min, sec, deciSec;
    aocommon::FitsWriter::ModifiedJulianDateToYMD(obs_date_, year, month, day);
    aocommon::FitsWriter::MJDToHMS(obs_date_, hour, min, sec, deciSec);
    char dateStr[40];
    std::sprintf(dateStr, "%d-%02d-%02dT%02d:%02d:%02d.%01d", year, month, day,
                 hour, min, sec, deciSec);
    fits_write_key(fptr, TSTRING, "DATE-OBS", (void*)dateStr, "", &status);
    checkStatus(status, filename);

    // Store the target's right ascension and declination
    double ra_deg = (object_ra_ / M_PI) * 180.0,
           dec_deg = (object_dec_ / M_PI) * 180.0;
    fits_write_key(fptr, TDOUBLE, "OBJ-RA", (void*)&ra_deg,
                   "Target right ascension in degrees", &status);
    checkStatus(status, filename);
    fits_write_key(fptr, TDOUBLE, "OBJ-DEC", (void*)&dec_deg,
                   "Target declination in degrees", &status);
    checkStatus(status, filename);
  }

  void writeImage(fitsfile* fptr, const std::string& filename,
                  const double* image, long* currentPixel) const {
    double null_value = std::numeric_limits<double>::max();
    int status = 0;
    const size_t total_size = n_times_ * n_channels_ * n_stokes_;
    fits_write_pixnull(fptr, TDOUBLE, currentPixel, total_size,
                       const_cast<double*>(image), &null_value, &status);
    checkStatus(status, filename);
  }

  void writeImage(fitsfile* fptr, const std::string& filename,
                  const float* image, long* currentPixel) const {
    float null_value = std::numeric_limits<float>::max();
    int status = 0;
    const size_t total_size = n_times_ * n_channels_ * n_stokes_;
    fits_write_pixnull(fptr, TFLOAT, currentPixel, total_size,
                       const_cast<float*>(image), &null_value, &status);
    checkStatus(status, filename);
  }

  template <typename NumType>
  void writeImage(fitsfile* fptr, const std::string& filename,
                  const NumType* image, long* currentPixel) const {
    double null_value = std::numeric_limits<double>::max();
    int status = 0;
    const size_t total_size = n_times_ * n_channels_ * n_stokes_;
    std::vector<double> copy(total_size);
    for (size_t i = 0; i != total_size; ++i) copy[i] = image[i];
    fits_write_pixnull(fptr, TDOUBLE, currentPixel, total_size, &copy[0],
                       &null_value, &status);
    checkStatus(status, filename);
  }

  fitsfile* file_pointer = nullptr;
  std::size_t n_times_;
  std::size_t n_channels_;
  std::size_t n_stokes_;
  double object_ra_;
  double object_dec_;
  double obs_date_;
  double time_resolution_;
  double frequency_;
  double frequency_resolution_;
  std::string telescope_name_;
  std::string observer_;
  std::string object_name_;
  std::string origin_;
  std::string origin_comment_;
  std::string obs_id_;
};
}  // namespace dp3::common

#endif
