// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MSReorder.h"

#include <iomanip>

#include <boost/filesystem/path.hpp>

namespace dp3 {
namespace reorder {
std::string GetFilenamePrefix(const std::string& ms_path_str,
                              const std::string& temp_dir) {
  boost::filesystem::path prefix_path;
  if (temp_dir.empty())
    prefix_path = ms_path_str;
  else {
    boost::filesystem::path ms_path(ms_path_str);
    prefix_path = boost::filesystem::path(temp_dir) / ms_path.filename();
  }
  const std::string prefix(prefix_path.remove_trailing_separator().string());
  return prefix;
}

std::string GetPartPrefix(const std::string& ms_path_str, size_t part_index,
                          aocommon::PolarizationEnum pol, size_t data_desc_id,
                          const std::string& temp_dir) {
  const std::string prefix = GetFilenamePrefix(ms_path_str, temp_dir);

  std::ostringstream part_prefix;
  part_prefix << prefix << "-part" << std::setw(4) << std::setfill('0')
              << part_index << "-"
              << aocommon::Polarization::TypeToShortString(pol) << "-b"
              << data_desc_id;
  return part_prefix.str();
}

std::string GetMetaFilename(const std::string& ms_path_str,
                            const std::string& temp_dir, size_t data_desc_id) {
  std::string prefix = GetFilenamePrefix(ms_path_str, temp_dir);

  std::ostringstream s;
  s << prefix << "-spw" << data_desc_id << "-parted-meta.tmp";
  return s.str();
}

void ExtractData(std::complex<float>* dest, size_t start_channel,
                 size_t end_channel,
                 const std::set<aocommon::PolarizationEnum>& pols_in,
                 const std::complex<float>* data,
                 aocommon::PolarizationEnum pol_out) {
  const size_t pol_count = pols_in.size();
  const std::complex<float>* in_ptr = data + start_channel * pol_count;
  const size_t selected_channel_count = end_channel - start_channel;

  if (pol_out == aocommon::Polarization::Instrumental) {
    if (pols_in.size() != 4) {
      throw std::runtime_error(
          "This mode requires the four polarizations to be present in the "
          "measurement set");
    }
    for (size_t ch = 0; ch != selected_channel_count * pols_in.size(); ++ch) {
      if (IsCFinite(*in_ptr))
        dest[ch] = *in_ptr;
      else
        dest[ch] = 0;
      ++in_ptr;
    }
  } else if (pol_out == aocommon::Polarization::DiagonalInstrumental) {
    if (pols_in.size() == 4) {
      size_t ch = 0;
      while (ch != selected_channel_count * 2) {
        if (IsCFinite(*in_ptr))
          dest[ch] = *in_ptr;
        else
          dest[ch] = 0;
        in_ptr += 3;  // jump from xx to yy
        ++ch;
        if (IsCFinite(*in_ptr))
          dest[ch] = *in_ptr;
        else
          dest[ch] = 0;
        ++in_ptr;
        ++ch;
      }
    } else if (pols_in.size() == 2) {
      for (size_t ch = 0; ch != selected_channel_count * 2; ++ch) {
        if (IsCFinite(*in_ptr))
          dest[ch] = *in_ptr;
        else
          dest[ch] = 0;
        ++in_ptr;
      }
    } else
      throw std::runtime_error(
          "Diagonal instrument visibilities requested, but this requires 2 or "
          "4 polarizations in the data");
  } else if (size_t pol_index;
             aocommon::Polarization::TypeToIndex(pol_out, pols_in, pol_index)) {
    in_ptr += pol_index;
    for (size_t ch = 0; ch != selected_channel_count; ++ch) {
      if (IsCFinite(*in_ptr))
        dest[ch] = *in_ptr;
      else
        dest[ch] = 0;
      in_ptr += pol_count;
    }
  } else {
    // Copy the right visibilities with conversion if necessary.
    switch (pol_out) {
      case aocommon::Polarization::StokesI: {
        size_t pol_index_a = 0, pol_index_b = 0;
        const bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_in, pol_index_a);
        const bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_in, pol_index_b);
        if (!has_XX || !has_YY) {
          const bool has_RR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RR, pols_in, pol_index_a);
          const bool has_LL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LL, pols_in, pol_index_b);
          if (!has_RR || !has_LL)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes I) from available "
                "polarizations");
        }

        for (size_t ch = 0; ch != selected_channel_count; ++ch) {
          in_ptr += pol_index_a;
          std::complex<float> val = *in_ptr;
          in_ptr += pol_index_b - pol_index_a;

          // I = (XX + YY) / 2
          val = (*in_ptr + val) * 0.5f;

          if (IsCFinite(val))
            dest[ch] = val;
          else
            dest[ch] = 0.0;

          in_ptr += pol_count - pol_index_b;
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        size_t pol_index_a = 0, pol_index_b = 0;
        const bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_in, pol_index_a);
        const bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_in, pol_index_b);
        if (has_XX && has_YY) {
          // Convert to StokesQ from XX and YY
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // Q = (XX - YY)/2
            val = (val - *in_ptr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        } else {
          // Convert to StokesQ from RR and LL
          const bool has_RL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RL, pols_in, pol_index_a);
          const bool has_LR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LR, pols_in, pol_index_b);
          if (!has_RL || !has_LR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes Q) from available "
                "polarizations");
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // Q = (RL + LR)/2
            val = (*in_ptr + val) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        }
      } break;
      case aocommon::Polarization::StokesU: {
        size_t pol_index_a = 0, pol_index_b = 0;
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_in, pol_index_b);
        if (has_XY && has_YX) {
          // Convert to StokesU from XY and YX
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // U = (XY + YX)/2
            val = (val + *in_ptr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        } else {
          // Convert to StokesU from RR and LL
          const bool has_RL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RL, pols_in, pol_index_a);
          const bool has_LR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LR, pols_in, pol_index_b);
          if (!has_RL || !has_LR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes U) from available "
                "polarizations");
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // U = -i (RL - LR)/2
            val = (val - *in_ptr) * 0.5f;
            val = std::complex<float>(val.imag(), -val.real());

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        }
      } break;
      case aocommon::Polarization::StokesV: {
        size_t pol_index_a = 0, pol_index_b = 0;
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_in, pol_index_b);
        if (has_XY && has_YX) {
          // Convert to StokesV from XX and YY
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // V = -i(XY - YX)/2
            val = (val - *in_ptr) * 0.5f;
            val = std::complex<float>(val.imag(), -val.real());

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        } else {
          // Convert to StokesV from RR and LL
          const bool has_RL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RR, pols_in, pol_index_a);
          const bool has_LR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LL, pols_in, pol_index_b);
          if (!has_RL || !has_LR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes V) from available "
                "polarizations");
          for (size_t ch = 0; ch != selected_channel_count; ++ch) {
            in_ptr += pol_index_a;
            std::complex<float> val = *in_ptr;
            in_ptr += pol_index_b - pol_index_a;

            // V = (RR - LL)/2
            val = (val - *in_ptr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            in_ptr += pol_count - pol_index_b;
          }
        }
      } break;
      default:
        throw std::runtime_error(
            "Could not convert ms polarizations to requested polarization");
    }
  }
}

template <>
void ExtractWeights<float>(float* dest, size_t start_channel,
                           size_t end_channel,
                           const std::set<aocommon::PolarizationEnum>& pols_in,
                           const std::complex<float>* data,
                           const float* weights, const bool* flags,
                           aocommon::PolarizationEnum pol_out) {
  const size_t pol_count = pols_in.size();
  const std::complex<float>* data_ptr = data + start_channel * pol_count;
  const float* weight_ptr = weights + start_channel * pol_count;
  const bool* flag_ptr = flags + start_channel * pol_count;
  const size_t selected_channel_count = end_channel - start_channel;

  size_t pol_index;
  if (pol_out == aocommon::Polarization::Instrumental) {
    for (size_t ch = 0; ch != selected_channel_count * pols_in.size(); ++ch) {
      if (!*flag_ptr && IsCFinite(*data_ptr))
        // The factor of 4 is to be consistent with StokesI
        // It is for having conjugate visibilities and because IDG doesn't
        // separately count XX and YY visibilities
        dest[ch] = *weight_ptr * 4.0f;
      else
        dest[ch] = 0.0f;
      data_ptr++;
      weight_ptr++;
      flag_ptr++;
    }
  } else if (pol_out == aocommon::Polarization::DiagonalInstrumental) {
    if (pols_in.size() == 4) {
      size_t ch = 0;
      while (ch != selected_channel_count * 2) {
        if (!*flag_ptr && IsCFinite(*data_ptr))
          // See explanation above for factor of 4
          dest[ch] = *weight_ptr * 4.0f;
        else
          dest[ch] = 0.0f;
        data_ptr += 3;  // jump from xx to yy
        weight_ptr += 3;
        flag_ptr += 3;
        ++ch;
        if (!*flag_ptr && IsCFinite(*data_ptr))
          dest[ch] = *weight_ptr * 4.0f;
        else
          dest[ch] = 0.0f;
        ++data_ptr;
        ++weight_ptr;
        ++flag_ptr;
        ++ch;
      }
    } else if (pols_in.size() == 2) {
      for (size_t ch = 0; ch != selected_channel_count * 2; ++ch) {
        if (!*flag_ptr && IsCFinite(*data_ptr))
          dest[ch] = *weight_ptr * 4.0f;
        else
          dest[ch] = 0.0f;
        ++data_ptr;
        ++weight_ptr;
        ++flag_ptr;
      }
    }
  } else if (aocommon::Polarization::TypeToIndex(pol_out, pols_in, pol_index)) {
    data_ptr += pol_index;
    weight_ptr += pol_index;
    flag_ptr += pol_index;
    for (size_t ch = 0; ch != selected_channel_count; ++ch) {
      if (!*flag_ptr && IsCFinite(*data_ptr))
        dest[ch] = *weight_ptr;
      else
        dest[ch] = 0.0f;
      data_ptr += pol_count;
      weight_ptr += pol_count;
      flag_ptr += pol_count;
    }
  } else {
    size_t pol_index_a = 0, pol_index_b = 0;
    switch (pol_out) {
      case aocommon::Polarization::StokesI: {
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_in, pol_index_b);
        if (!has_XY || !has_YX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              pols_in, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              pols_in, pol_index_b);
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        const bool has_XX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, pols_in, pol_index_a);
        const bool has_YY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, pols_in, pol_index_b);
        if (!has_XX || !has_YY) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              pols_in, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              pols_in, pol_index_b);
        }
      } break;
      case aocommon::Polarization::StokesU: {
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_in, pol_index_b);
        if (!has_XY || !has_YX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              pols_in, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              pols_in, pol_index_b);
        }
      } break;
      case aocommon::Polarization::StokesV: {
        const bool has_XY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, pols_in, pol_index_a);
        const bool has_YX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, pols_in, pol_index_b);
        if (!has_XY || !has_YX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              pols_in, pol_index_a);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              pols_in, pol_index_b);
        }
      } break;
      default:
        throw std::runtime_error(
            "Could not convert ms polarizations to requested polarization");
        break;
    }

    weight_ptr += pol_index_a;
    data_ptr += pol_index_a;
    flag_ptr += pol_index_a;
    for (size_t ch = 0; ch != selected_channel_count; ++ch) {
      float w;
      if (!*flag_ptr && IsCFinite(*data_ptr))
        w = *weight_ptr * 4.0f;
      else
        w = 0.0f;
      data_ptr += pol_index_b - pol_index_a;
      weight_ptr += pol_index_b - pol_index_a;
      flag_ptr += pol_index_b - pol_index_a;
      if (!*flag_ptr && IsCFinite(*data_ptr))
        w = std::min<float>(w, *weight_ptr * 4.0f);
      else
        w = 0.0f;
      dest[ch] = w;
      weight_ptr += pol_count - pol_index_b + pol_index_a;
      data_ptr += pol_count - pol_index_b + pol_index_a;
      flag_ptr += pol_count - pol_index_b + pol_index_a;
    }
  }
}

}  // namespace reorder
}  // namespace dp3
