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

#ifndef IDGPREDICT_PARSETATERMS_H
#define IDGPREDICT_PARSETATERMS_H

#include <EveryBeam/aterms/parsetprovider.h>

namespace DP3 {
namespace DPPP {

/**
 * @brief Parses the parameter settings (parset) related to the
 * aterm settings in an EveryBeam-acceptable format.
 *
 */
class ParsetATerms : public everybeam::aterms::ParsetProvider {
 public:
  /// Constructor.
  /// Note that this class stores references to the provided parset and prefix.
  /// If the parset or prefix object change, this class will use the new values.
  /// If the parset or prefix object become invalid, this class should not be
  /// used anymore.
  explicit ParsetATerms(const ParameterSet& parset,
                        const std::string& prefix = "")
      : parset_(parset), prefix_(prefix) {}

  /// Extracts string for given key.
  std::string GetString(const std::string& key) const final override {
    return parset_.getString(prefix_ + key);
  }

  /// Extracts string for given \p key. Defaults to \p or_value.
  std::string GetStringOr(const std::string& key,
                          const std::string& or_value) const final override {
    return parset_.getString(prefix_ + key, or_value);
  }

  /// Extracts string list for given \p key.
  std::vector<std::string> GetStringList(
      const std::string& key) const final override {
    return parset_.getStringVector(prefix_ + key);
  }

  /// Extracts double for given \p key. Defaults to \p or_value.
  double GetDoubleOr(const std::string& key,
                     double or_value) const final override {
    return parset_.getDouble(prefix_ + key, or_value);
  };

  /// Extracts bool for given \p key. Defaults to \p or_value.
  bool GetBool(const std::string& key) const final override {
    return parset_.getBool(prefix_ + key);
  };

  /// Extracts bool for given \p key. Defaults to \p or_value.
  bool GetBoolOr(const std::string& key, bool or_value) const final override {
    return parset_.getBool(prefix_ + key, or_value);
  };

 private:
  const ParameterSet& parset_;
  const std::string& prefix_;
};
}  // namespace DPPP
}  // namespace DP3
#endif