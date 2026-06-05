#ifndef DP3_DDECAL_INITIAL_SOLUTIONS_H_
#define DP3_DDECAL_INITIAL_SOLUTIONS_H_

#include <memory>
#include <string>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include <schaapcommon/h5parm/jonesparameters.h>
#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/soltab.h>

#include "Settings.h"

#include "common/ParameterSet.h"

namespace dp3::ddecal {

class InitialSolutions {
 public:
  InitialSolutions(const common::ParameterSet& parset,
                   const std::string& prefix);

  void Initialize(size_t n_directions,
                  size_t n_polarization_parameters_per_solution);

  void Read(const Settings& settings,
            std::vector<std::vector<std::complex<double>>>& solutions_data,
            const std::vector<std::string>& antenna_names,
            const std::vector<double>& frequencies,
            const std::vector<base::Direction>& source_directions,
            double average_time);

  bool Empty() const { return initial_solutions_h5_parm_name_.empty(); }

  void Show(std::ostream& os) const;

 private:
  xt::xtensor<std::complex<float>, 3> ReadJonesMatrix(
      const base::Direction& direction, double timestamp,
      schaapcommon::h5parm::GainType gain_type,
      schaapcommon::h5parm::SolTab* first_soltab,
      schaapcommon::h5parm::SolTab* second_soltab,
      const std::vector<double>& frequencies,
      const std::vector<std::string>& antenna_names);

  /// Stores the H5Parm file and loads all solutions into memory when the user
  /// requests the solver to use initial solutions.
  std::unique_ptr<schaapcommon::h5parm::H5Parm> solutions_;
  std::string initial_solutions_h5_parm_name_;
  std::vector<std::string> initial_solutions_table_;
  std::vector<schaapcommon::h5parm::SolTab> solution_tables_;
  bool initial_solutions_are_full_jones_ = false;
  schaapcommon::h5parm::JonesParameters::InterpolationType interpolation_type_;
  schaapcommon::h5parm::JonesParameters::MissingAntennaBehavior
      missing_antenna_behavior_;
  std::vector<schaapcommon::h5parm::GainType> gain_types_;

  size_t n_directions_ = 0;
  size_t n_polarization_parameters_per_solution_ = 0;
};

}  // namespace dp3::ddecal

#endif
