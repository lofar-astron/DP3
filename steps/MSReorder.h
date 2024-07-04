// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_MSREORDER_H_
#define DP3_STEPS_MSREORDER_H_

#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>

#include <aocommon/polarization.h>

namespace dp3 {
namespace reorder {
#pragma pack(push, 1)
struct MetaRecordBuffer {
  double u;
  double v;
  double w;
  double time;
  uint16_t antenna1;
  uint16_t antenna2;
  uint16_t field_id;
};
struct PartHeaderBuffer {
  uint64_t channel_count;
  uint64_t channel_start;
  uint32_t data_desc_id;
  bool has_model;
};
struct MetaHeaderBuffer {
  double start_time;
  uint64_t selected_row_count;
  uint32_t filename_length;
};
#pragma pack(pop)

template <typename NumType>
static bool IsCFinite(const std::complex<NumType>& c) {
  return std::isfinite(c.real()) && std::isfinite(c.imag());
}

struct MetaHeader {
  double start_time = 0.0;
  uint64_t selected_row_count = 0;
  uint32_t filename_length = 0;
  void Read(std::istream& str) {
    MetaHeaderBuffer meta_header_buffer;
    str.read(reinterpret_cast<char*>(&meta_header_buffer),
             sizeof(MetaHeaderBuffer));
    start_time = meta_header_buffer.start_time;
    selected_row_count = meta_header_buffer.selected_row_count;
    filename_length = meta_header_buffer.filename_length;
  }
  void Write(std::ostream& str) const {
    MetaHeaderBuffer meta_header_buffer{start_time, selected_row_count,
                                        filename_length};
    str.write(reinterpret_cast<const char*>(&meta_header_buffer),
              sizeof(MetaHeaderBuffer));
  }
  static constexpr size_t BINARY_SIZE =
      sizeof(start_time) + sizeof(selected_row_count) + sizeof(filename_length);
  static_assert(BINARY_SIZE == 20);
};

struct MetaRecord {
  double u = 0.0, v = 0.0, w = 0.0, time = 0.0;
  uint16_t antenna1 = 0, antenna2 = 0, field_id = 0;
  static constexpr size_t BINARY_SIZE =
      sizeof(double) * 4 + sizeof(uint16_t) * 3;
  static_assert(BINARY_SIZE == 38);
  void Read(std::istream& str) {
    MetaRecordBuffer meta_record_buffer;
    str.read(reinterpret_cast<char*>(&meta_record_buffer),
             sizeof(MetaRecordBuffer));

    u = meta_record_buffer.u;
    v = meta_record_buffer.v;
    w = meta_record_buffer.w;
    time = meta_record_buffer.time;
    antenna1 = meta_record_buffer.antenna1;
    antenna2 = meta_record_buffer.antenna2;
    field_id = meta_record_buffer.field_id;
  }
  void Write(std::ostream& str) const {
    MetaRecordBuffer meta_record_buffer{u,        v,        w,       time,
                                        antenna1, antenna2, field_id};
    str.write(reinterpret_cast<const char*>(&meta_record_buffer),
              sizeof(MetaRecordBuffer));
  }
};

struct PartHeader {
  uint64_t channel_count = 0;
  uint64_t channel_start = 0;
  uint32_t data_desc_id = 0;
  bool has_model = false;
  static constexpr size_t BINARY_SIZE =
      sizeof(channel_count) + sizeof(channel_start) + sizeof(data_desc_id) +
      sizeof(has_model);
  static_assert(BINARY_SIZE == 21);
  void Read(std::istream& str) {
    PartHeaderBuffer part_header_buffer;
    str.read(reinterpret_cast<char*>(&part_header_buffer),
             sizeof(PartHeaderBuffer));
    channel_count = part_header_buffer.channel_count;
    channel_start = part_header_buffer.channel_start;
    data_desc_id = part_header_buffer.data_desc_id;
    has_model = part_header_buffer.has_model;
  }
  void Write(std::ostream& str) const {
    PartHeaderBuffer part_header_buffer{channel_count, channel_start,
                                        data_desc_id, has_model};
    str.write(reinterpret_cast<const char*>(&part_header_buffer),
              sizeof(PartHeaderBuffer));
  }
};

std::string GetMetaFilename(const std::string& ms_path_str,
                            const std::string& temp_dir, size_t data_desc_id);
std::string GetFilenamePrefix(const std::string& ms_path_str,
                              const std::string& temp_dir);

std::string GetPartPrefix(const std::string& ms_path_str, size_t part_index,
                          aocommon::PolarizationEnum pol, size_t data_desc_id,
                          const std::string& temp_dir);

void ExtractData(std::complex<float>* dest, size_t start_channel,
                 size_t end_channel,
                 const std::set<aocommon::PolarizationEnum>& pols_in,
                 const std::complex<float>* data,
                 aocommon::PolarizationEnum pol_out);

template <typename NumType>
void ExtractWeights(NumType* dest, size_t start_channel, size_t end_channel,
                    const std::set<aocommon::PolarizationEnum>& pols_in,
                    const std::complex<float>* data, const float* weights,
                    const bool* flags, aocommon::PolarizationEnum pol_out);

template <>
void ExtractWeights(float* dest, size_t start_channel, size_t end_channel,
                    const std::set<aocommon::PolarizationEnum>& pols_in,
                    const std::complex<float>* data, const float* weights,
                    const bool* flags, aocommon::PolarizationEnum pol_out);
}  // namespace reorder

}  // namespace dp3

#endif
