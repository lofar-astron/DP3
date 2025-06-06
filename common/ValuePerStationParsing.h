// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_COMMON_VALUE_PER_STATION_PARSING_H_
#define DP3_COMMON_VALUE_PER_STATION_PARSING_H_

#include <regex>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>

#include <aocommon/logger.h>
#include <aocommon/uvector.h>

#include "ParameterValue.h"
#include "StringTools.h"

namespace dp3::common {

/**
 * Helper function for @ref ParseValuePerStation(). This function finds all
 * stations that match a given pattern, and assigns the associated values for
 * these stations to the given value. Example: AssignStationValues({0,0,0}, 3,
 * {"a*"}, {"aa", "ab", "b"}) will set the result values to {3,3,0}.
 */
template <typename ValueType>
void AssignStationValues(std::span<ValueType> result, ValueType value,
                         const std::string pattern_string,
                         std::span<const std::string> station_names) {
  // This function is based on StationAdder::GetMatchingStations(). Because we
  // don't need to have the indices of the stations, that function is slightly
  // adapted here.
  aocommon::UVector<bool> is_selected(station_names.size(), false);
  const std::vector<std::string> patterns(
      ParameterValue(pattern_string).getStringVector());
  for (const std::string& pattern : patterns) {
    int n = 0;
    const auto SetSelected = [&](bool new_value, const std::string& pattern) {
      const std::regex pattern_regex(pattern);
      for (size_t i = 0; i < station_names.size(); ++i) {
        if (std::regex_match(station_names[i], pattern_regex)) {
          is_selected[i] = new_value;
          n++;
        }
      }
    };
    if (pattern.size() > 1 && (pattern[0] == '!' || pattern[0] == '^')) {
      SetSelected(false, PatternToRegex(pattern.substr(1)));
    } else {
      SetSelected(true, PatternToRegex(pattern));
    }
    if (n == 0) {
      aocommon::Logger::Warn << "No matching stations found for pattern "
                             << pattern << '\n';
    }
  }
  for (size_t i = 0; i < is_selected.size(); ++i) {
    if (is_selected[i]) result[i] = value;
  }
}

/**
 * This function parses a configuration string list that specifies a value per
 * station. Wildcards like "CS*" can be used to ease specifying multiple
 * stations. The station pattern and the value are separated by a colon. In @p
 * result, the value for stations for which no pattern is matching will be
 * unchanged, hence the @p result value should be initialized to sensible
 * default values before this call.
 * @param [in,out] result On input, an array equal in size to station_names with
 * the default values for the stations. On output, it will be set according to
 * @p string_list.
 * @param string_list A list of values in the form of patterns:value. The
 * 'patterns' part can be an array enclosed by square brackets, e.g. [RS_017,
 * CS*, ^CS_001]:3 would set the value to 3 for the stations RS_017, stations
 * whose name start with the text "CS", but not for CS_001.
 */
template <typename ValueType>
void ParseValuePerStation(std::span<ValueType> result,
                          std::span<const std::string> string_list,
                          std::span<const std::string> station_names) {
  assert(result.size() == station_names.size());

  for (const std::string& stations_value : string_list) {
    const size_t colon = stations_value.find(':');
    if (colon == std::string::npos)
      throw std::runtime_error(
          "Items of station-value lists must be of the form "
          "'station-pattern:value'. Item is missing a colon: '" +
          stations_value + "'");
    const std::string value_string(
        stations_value.substr(colon + 1, stations_value.size() - colon));
    std::istringstream value_stream(value_string);
    ValueType value;
    value_stream >> value;
    std::string station_pattern(stations_value.substr(0, colon));
    AssignStationValues(result, value, station_pattern, station_names);
  }
}

}  // namespace dp3::common

#endif
