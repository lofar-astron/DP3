// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This class is responsible to parse the LOFAR antenna configuration file.
/// Such configuration contains the antenna positions for a specific aartfaac
/// configration but also the orientation matrix of the common reference field

/// This code was taken from
/// https://git.astron.nl/RD/aartfaac-tools/-/raw/master/lib/aartfaac2ms/antennaconfig.h?ref_type=heads
/// and adapted to fit into the DP3 codebase.

#ifndef AARTFAACREADER_ANTENNACONFIG_H_
#define AARTFAACREADER_ANTENNACONFIG_H_

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <locale>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <base/RcuMode.h>
#include <casacore/measures/Measures/MPosition.h>

namespace {

// trim from start (in place)
void LTrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                  [](int ch) { return !std::isspace(ch); }));
}

// trim from end (in place)
void RTrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](int ch) { return !std::isspace(ch); })
              .base(),
          s.end());
}

// trim from both ends (in place)
void Trim(std::string &s) {
  LTrim(s);
  RTrim(s);
}
}  // namespace

namespace dp3::aartfaacreader {

class AntennaConfig {
 public:
  explicit AntennaConfig(const char *filename) : file_(filename) {
    ParseFile();
  }

  std::vector<casacore::MPosition> GetLBAPositions() const {
    return GetPositions("LBA");
  }

  std::vector<casacore::MPosition> GetHBAPositions() const {
    return GetPositions("HBA");
  }

  std::array<double, 9> GetLBAAxes() const {
    return GetAxes("LBA_ROTATION_MATRIX");
  }

  std::array<double, 9> GetHBA0Axes() const {
    return GetAxes("HBA0_ROTATION_MATRIX");
  }

  std::array<double, 9> GetHBA1Axes() const {
    return GetAxes("HBA1_ROTATION_MATRIX");
  }

  std::array<double, 9> GetAxesFromMode(base::RcuMode mode) {
    switch (mode.mode) {
      case base::RcuMode::LBAInner10_90:
      case base::RcuMode::LBAInner30_90:
      case base::RcuMode::LBAOuter10_90:
      case base::RcuMode::LBAOuter30_90:
        return GetLBAAxes();
      case base::RcuMode::HBA110_190:
      case base::RcuMode::HBA170_230:
      case base::RcuMode::HBA210_270:
        return GetHBA0Axes();
        break;
      default:
        throw std::runtime_error("Wrong RCU mode");
    }
  }

  std::vector<casacore::MPosition> GetArrayFromMode(base::RcuMode mode) {
    switch (mode.mode) {
      case base::RcuMode::LBAInner10_90:
      case base::RcuMode::LBAInner30_90:
      case base::RcuMode::LBAOuter10_90:
      case base::RcuMode::LBAOuter30_90:
        return GetLBAPositions();
        break;
      case base::RcuMode::HBA110_190:
      case base::RcuMode::HBA170_230:
      case base::RcuMode::HBA210_270:
        return GetHBAPositions();
      default:
        throw std::runtime_error("Wrong RCU mode");
    }
  }

 private:
  struct Array {
    std::string name;
    std::string band;
    std::vector<double> data;
  };

  void ParseFile() {
    Next();
    struct Array a;
    while (ReadArray(a.name, a.band, a.data)) {
      std::string key;
      if (a.band.empty())
        key = a.name;
      else
        key = a.band + "_" + a.name;
      values_.insert(std::make_pair(key, a));
    }
  }
  const std::vector<double> &GetArray(const std::string &name) const {
    return values_.find(name)->second.data;
  }
  std::vector<casacore::MPosition> GetPositions(
      const std::string &arrayName) const {
    const std::vector<double> &arr = GetArray(arrayName);
    std::vector<casacore::MPosition> position;
    for (size_t index = 0; index < arr.size(); index += 6) {
      position.emplace_back(casacore::MPosition{
          casacore::MVPosition{arr[index], arr[index + 1], arr[index + 2]},
          casacore::MPosition::ITRF});
    }
    return position;
  }

  std::array<double, 9> GetAxes(const std::string &arrayName) const {
    const std::vector<double> &arr = GetArray(arrayName);
    if (arr.size() != 9)
      throw std::runtime_error(
          "The array for coordinate axes in the antenna "
          "config file had an incorrect size");

    std::array<double, 9> axes;
    for (size_t index = 0; index < 9; ++index) axes[index] = arr[index];
    return axes;
  }

  bool Next() {
    if (line_.empty() || line_position_ >= line_.size()) {
      do {
        std::getline(file_, line_);
        if (!file_) {
          token_.clear();
          line_.clear();
          return false;
        }
        Trim(line_);
      } while (line_.empty() || line_[0] == '#');
      line_position_ = 0;
    }
    size_t pos = line_.find_first_of(" \t", line_position_);
    if (pos == std::string::npos) {
      token_ = line_.substr(line_position_);
      line_.clear();
      line_position_ = 0;
      return true;
    } else {
      token_ = line_.substr(line_position_, pos - line_position_);
      line_position_ = pos + 1;
      while (line_position_ < line_.size() &&
             (line_[line_position_] == ' ' || line_[line_position_] == '\t'))
        ++line_position_;
      Trim(token_);
      if (token_.empty())
        return Next();
      else
        return true;
    }
  }

  std::vector<int> ReadDimensions() {
    std::vector<int> dimensions = {std::atoi(token_.c_str())};
    do {
      if (!Next())
        throw std::runtime_error(
            "Antenna config file has bad format: expected dimensions");
      if (token_ == "x") {
        if (!Next())
          throw std::runtime_error(
              "Antenna config file has bad format: "
              "expected another dimension after x");
        int dimension_value = std::atoi(token_.c_str());
        dimensions.push_back(dimension_value);
      } else if (token_ != "[") {
        throw std::runtime_error("Antenna config file has bad format");
      }
    } while (token_ != "[");
    return dimensions;
  }

  const std::vector<double> ReadData(const std::vector<int> &dimensions) {
    int count = std::accumulate(dimensions.begin(), dimensions.end(), 1,
                                [](int a, int b) { return a * b; });
    std::vector<double> values(count);
    for (int i = 0; i != count; ++i) {
      if (!Next()) {
        throw std::runtime_error("Missing numbers");
      }
      values[i] = std::atof(token_.c_str());
    }
    Next();  // move TO ']'
    return values;
  }

  bool ReadArray(std::string &name, std::string &band,
                 std::vector<double> &values) {
    values.clear();
    if (!std::isalpha(token_[0])) return false;
    name = token_;
    Next();

    if (std::isalpha(token_[0])) {
      band = token_;
      Next();
    } else {
      band.clear();
    }

    std::vector<int> dimensions1 = ReadDimensions();
    std::vector<double> data1 = ReadData(dimensions1);
    if (Next() && token_[0] >= '0' && token_[0] <= '9') {
      std::vector<int> dimensions2 = ReadDimensions();
      std::vector<double> data2 = ReadData(dimensions2);
      Next();  // skip ']'

      values = std::move(data2);
      for (size_t i = 0; i != values.size(); ++i)
        values[i] += data1[i % data1.size()];
    } else {
      values = std::move(data1);
    }
    return true;
  }

  std::map<std::string, Array> values_;

  std::ifstream file_;
  std::string line_;
  size_t line_position_;
  std::string token_;
};
}  // namespace dp3::aartfaacreader

#endif  // AARTFAACREADER_ANTENNACONFIG_H_
