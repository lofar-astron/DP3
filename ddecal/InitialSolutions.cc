#include "InitialSolutions.h"

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/views/xview.hpp>

using schaapcommon::h5parm::GainType;
using schaapcommon::h5parm::JonesParameters;

namespace dp3::ddecal {

void InitialSolutions::Initialize(
    size_t n_directions, size_t n_polarization_parameters_per_solution) {
  n_directions_ = n_directions;
  n_polarization_parameters_per_solution_ =
      n_polarization_parameters_per_solution;
}

InitialSolutions::InitialSolutions(const common::ParameterSet& parset,
                                   const std::string& prefix)
    : initial_solutions_h5_parm_name_(
          parset.getString(prefix + "initialsolutions.h5parm", "")),
      missing_antenna_behavior_(
          JonesParameters::StringToMissingAntennaBehavior(parset.getString(
              prefix + "initialsolutions.missingantennabehavior", "error"))) {
  if (initial_solutions_h5_parm_name_.empty()) {
    return;
  }

  const std::vector<std::string> default_solution_tables = {"amplitude000",
                                                            "phase000"};
  initial_solutions_table_ = parset.getStringVector(
      prefix + "initialsolutions.soltab", default_solution_tables);
  solutions_ = std::make_unique<schaapcommon::h5parm::H5Parm>(
      initial_solutions_h5_parm_name_, initial_solutions_table_);

  for (const std::string& soltab_name : initial_solutions_table_) {
    solution_tables_.push_back(solutions_->GetSolTab(soltab_name));
  }

  // Check if H5Parm stores full-Jones solutions.
  initial_solutions_are_full_jones_ = false;
  if (solution_tables_[0].HasAxis("pol") &&
      initial_solutions_table_.size() == 2 &&
      initial_solutions_table_[0].find("amplitude") != std::string::npos &&
      initial_solutions_table_[1].find("phase") != std::string::npos) {
    initial_solutions_are_full_jones_ =
        solution_tables_[0].GetAxis("pol").size == 4;
  }

  if (initial_solutions_are_full_jones_) {
    gain_types_.resize(1);
    gain_types_[0] = GainType::kFullJones;
  } else {
    gain_types_.reserve(initial_solutions_table_.size());
    for (const std::string& soltab_name : initial_solutions_table_) {
      gain_types_.push_back(JonesParameters::H5ParmTypeStringToGainType(
          solutions_->GetSolTab(soltab_name).GetType()));
    }
  }

  const std::string interpolation_type =
      parset.getString(prefix + "initialsolutions.interpolation", "nearest");
  if (interpolation_type == "nearest") {
    interpolation_type_ = JonesParameters::InterpolationType::NEAREST;
  } else if (interpolation_type == "linear") {
    interpolation_type_ = JonesParameters::InterpolationType::LINEAR;
  } else {
    throw std::runtime_error("Unsupported interpolation method: " +
                             interpolation_type);
  }
}

xt::xtensor<std::complex<float>, 3> InitialSolutions::ReadJonesMatrix(
    const base::Direction& direction, double timestamp,
    schaapcommon::h5parm::GainType gain_type,
    schaapcommon::h5parm::SolTab* first_soltab,
    schaapcommon::h5parm::SolTab* second_soltab,
    const std::vector<double>& frequencies,
    const std::vector<std::string>& antenna_names) {
  // Retrieve initial solutions for the patch that's closest to the
  // directions in the sourcedb, since the patch names and the number of
  // patches in the provided skymodel might not be the same as the
  // directions in the H5Parm.
  const std::vector<double> timestamps = {timestamp};
  const std::string closest_patch =
      solutions_->GetNearestSource(direction.ra, direction.dec);
  const hsize_t soltab_direction_index =
      first_soltab->GetDirIndex(closest_patch);
  size_t n_parm_values = 0;
  if (gain_type == GainType::kTec || gain_type == GainType::kDelay) {
    n_parm_values =
        first_soltab->HasAxis("pol") ? first_soltab->GetAxis("pol").size : 1;
  }
  auto jones_parameters =
      std::make_unique<schaapcommon::h5parm::JonesParameters>(
          frequencies, timestamps, antenna_names, gain_type,
          interpolation_type_, soltab_direction_index, first_soltab,
          second_soltab, false, n_parm_values, missing_antenna_behavior_);

  // Write casacore Cube with Jones parameters to xtensor.
  const casacore::Cube<std::complex<float>>& jones_matrix =
      jones_parameters->GetParms();
  // Ensure that the xtensor has shape: frequency x antenna x polarization.
  const std::array<size_t, 3> xt_shape = {
      static_cast<size_t>(jones_matrix.shape()[2]),
      static_cast<size_t>(jones_matrix.shape()[1]),
      static_cast<size_t>(jones_matrix.shape()[0])};
  xt::xtensor<std::complex<float>, 3> jones_tensor =
      xt::adapt(jones_matrix.tovector(), xt_shape);
  return jones_tensor;
}

void InitialSolutions::Read(
    const Settings& settings,
    std::vector<std::vector<std::complex<double>>>& solutions_data,
    const std::vector<std::string>& antenna_names,
    const std::vector<double>& frequencies,
    const std::vector<base::Direction>& source_directions,
    double average_time) {
  const size_t n_subsolutions = settings.GetNSolutions();
  const size_t n_values_per_channel_block =
      n_subsolutions * antenna_names.size() *
      n_polarization_parameters_per_solution_;
  const size_t n_frequencies = frequencies.size();
  const size_t n_antennas = antenna_names.size();

  const std::array<size_t, 4> xt_shape = {
      n_directions_, n_frequencies, n_antennas,
      n_polarization_parameters_per_solution_};
  xt::xtensor<std::complex<float>, 4> jones_parameters_per_direction(xt_shape);
  if (initial_solutions_are_full_jones_) {
    for (size_t direction_index = 0; direction_index < n_directions_;
         ++direction_index) {
      xt::view(jones_parameters_per_direction, direction_index, xt::all(),
               xt::all(), xt::all()) =
          ReadJonesMatrix(source_directions[direction_index], average_time,
                          gain_types_[0], &solution_tables_[0],
                          &solution_tables_[1], frequencies, antenna_names);
    }
  } else {
    // Load solutions per soltab and multiply to get the full Jones matrix for
    // each direction.
    const size_t n_soltabs = initial_solutions_table_.size();
    // The schaapcommon h5parm functionality always converts the data to two
    // polarizations. TODO it would be nicer if it wouldn't do this
    // unconditionally.
    const size_t kNPolarizationsInH5Parm = 2;
    const std::array<size_t, 5> xt_shape = {n_soltabs, n_directions_,
                                            n_frequencies, n_antennas,
                                            kNPolarizationsInH5Parm};
    xt::xtensor<std::complex<float>, 5> jones_parameters_per_soltab(xt_shape);
    for (size_t soltab_index = 0; soltab_index < n_soltabs; ++soltab_index) {
      for (size_t direction_index = 0; direction_index < n_directions_;
           ++direction_index) {
        xt::view(jones_parameters_per_soltab, soltab_index, direction_index,
                 xt::all(), xt::all(), xt::all()) =
            ReadJonesMatrix(source_directions[direction_index], average_time,
                            gain_types_[soltab_index],
                            &solution_tables_[soltab_index], nullptr,
                            frequencies, antenna_names);
      }
    }

    jones_parameters_per_direction = xt::prod(jones_parameters_per_soltab, {0});
  }

  for (size_t channel_block = 0; channel_block < n_frequencies;
       ++channel_block) {
    solutions_data[channel_block].resize(n_values_per_channel_block);
    for (size_t antenna_index = 0; antenna_index < n_antennas;
         ++antenna_index) {
      size_t n_assigned_subsolutions = 0;
      for (size_t direction_index = 0; direction_index < n_directions_;
           ++direction_index) {
        const xt::xtensor<std::complex<float>, 3>& jones_parameters =
            xt::view(jones_parameters_per_direction, direction_index, xt::all(),
                     xt::all(), xt::all());
        size_t n_subsolutions_per_direction =
            settings.sub_solutions_per_direction[direction_index];
        for (size_t direction_solution_index = 0;
             direction_solution_index < n_subsolutions_per_direction;
             ++direction_solution_index) {
          size_t direction_dependent_solution_index =
              n_assigned_subsolutions + direction_solution_index;
          for (size_t polarization_index = 0;
               polarization_index < n_polarization_parameters_per_solution_;
               ++polarization_index) {
            const size_t flattened_index =
                antenna_index * n_subsolutions *
                    n_polarization_parameters_per_solution_ +
                direction_dependent_solution_index *
                    n_polarization_parameters_per_solution_ +
                polarization_index;
            solutions_data[channel_block][flattened_index] = jones_parameters(
                channel_block, antenna_index, polarization_index);
          }
        }
        n_assigned_subsolutions += n_subsolutions_per_direction;
      }
    }
  }
}

void InitialSolutions::Show(std::ostream& os) const {
  if (!initial_solutions_h5_parm_name_.empty()) {
    os << "  initial sols H5Parm: " << initial_solutions_h5_parm_name_ << '\n'
       << "               soltab: ";
    for (size_t i = 0; i < initial_solutions_table_.size(); ++i) {
      os << initial_solutions_table_[i]
         << (i == initial_solutions_table_.size() - 1 ? " " : ", ");
    }
    os << '\n' << "                 type: ";
    for (size_t i = 0; i < gain_types_.size(); ++i) {
      os << JonesParameters::GainTypeToHumanReadableString(gain_types_[i])
         << (i == gain_types_.size() - 1 ? " " : ", ");
    }
    os << '\n'
       << "        interp method: "
       << (interpolation_type_ == JonesParameters::InterpolationType::NEAREST
               ? "nearest"
               : "linear")
       << '\n'
       << "          missing ant: "
       << JonesParameters::MissingAntennaBehaviorToString(
              missing_antenna_behavior_)
       << '\n';
  }
}

}  // namespace dp3::ddecal
