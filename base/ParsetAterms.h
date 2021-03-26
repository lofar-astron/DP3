// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef IDGPREDICT_PARSETATERMS_H
#define IDGPREDICT_PARSETATERMS_H

#include <EveryBeam/aterms/parsetprovider.h>

namespace dp3 {
namespace base {

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
  explicit ParsetATerms(const common::ParameterSet& parset,
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
  const common::ParameterSet& parset_;
  const std::string& prefix_;
};
}  // namespace base
}  // namespace dp3
#endif
