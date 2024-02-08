// ApplyBeam.h: DP3 step class to ApplyBeam visibilities from a source model
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to apply the beam model (optionally inverted)
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_APPLYBEAM_H_
#define DP3_STEPS_APPLYBEAM_H_

#include "InputStep.h"

#include <dp3/base/DPBuffer.h>

#include <EveryBeam/telescope/telescope.h>

#include <aocommon/matrix2x2.h>
#include <aocommon/barrier.h>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/measures/Measures/MDirection.h>

#include "../common/ParameterSet.h"

namespace dp3 {
namespace steps {

/// \brief DP3 step class to ApplyBeam visibilities from a source model

/// This class is a Step class to apply the beam model, optionally inverted.
/// The input MeasurementSet it operates on, must have the LOFAR subtables
/// defining the station layout and tiles/dipoles used.

class ApplyBeam : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  ApplyBeam(const common::ParameterSet&, const string& prefix,
            bool substep = false);

  common::Fields getRequiredFields() const override {
    common::Fields fields = kDataField;
    if (itsUpdateWeights) fields |= kWeightsField;
    return fields;
  }

  common::Fields getProvidedFields() const override {
    common::Fields fields = kDataField;
    if (itsUpdateWeights) fields |= kWeightsField;
    return fields;
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override {
    return processMultithreaded(std::move(buffer), 0);
  }

  /// If apply beam is called from multiple threads, it needs the thread index
  /// to determine what scratch space to use etc.
  bool processMultithreaded(std::unique_ptr<base::DPBuffer> buffer,
                            size_t thread);

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  bool invert() { return itsInvert; }

  template <typename T>
  static void applyBeam(const base::DPInfo& info, double time, T* data0,
                        float* weight0, const everybeam::vector3r_t& srcdir,
                        const everybeam::telescope::Telescope* telescope,
                        std::vector<aocommon::MC2x2>& beamValues, bool invert,
                        everybeam::CorrectionMode mode,
                        bool doUpdateWeights = false,
                        std::mutex* mutex = nullptr);
  // This method applies the beam for processing when parallelizing over
  // baselines. Because the beam is a per-antenna effect, this requires
  // synchronisation, which is performed with the provided barrier.
  static void applyBeam(const base::DPInfo& info, double time,
                        std::complex<double>* data0, float* weight0,
                        const everybeam::vector3r_t& srcdir,
                        const everybeam::telescope::Telescope* telescope,
                        std::vector<aocommon::MC2x2>& beam_values,
                        const std::pair<size_t, size_t>& baseline_range,
                        const std::pair<size_t, size_t>& station_range,
                        aocommon::Barrier& barrier, bool invert,
                        everybeam::CorrectionMode mode,
                        bool do_update_weights = false,
                        std::mutex* mutex = nullptr);

  template <typename T>
  static void applyBeamStokesIArrayFactor(
      const base::DPInfo& info, double time, T* data0,
      const everybeam::vector3r_t& srcdir,
      const everybeam::telescope::Telescope* telescope,
      std::vector<everybeam::complex_t>& beamValues, bool invert,
      everybeam::CorrectionMode mode, std::mutex* mutex = nullptr);

  // This method applies the Stokes I Array factor for processing when
  // parallelizing over baselines and uses a barrier for synchronisation,
  // similar to the corresponding @ref applyBeam() overload.
  static void applyBeamStokesIArrayFactor(
      const base::DPInfo& info, double time, std::complex<double>* data0,
      const everybeam::vector3r_t& srcdir,
      const everybeam::telescope::Telescope* telescope,
      std::vector<everybeam::complex_t>& beam_values,
      const std::pair<size_t, size_t>& baseline_range,
      const std::pair<size_t, size_t>& station_range,
      aocommon::Barrier& barrier, bool invert, everybeam::CorrectionMode mode,
      std::mutex* mutex = nullptr);

 private:
  everybeam::vector3r_t dir2Itrf(const casacore::MDirection& dir,
                                 casacore::MDirection::Convert& measConverter);

  string itsName;
  bool itsInvert;
  bool itsUpdateWeights;
  std::vector<string> itsDirectionStr;
  casacore::MDirection itsDirection;
  bool itsUseChannelFreq;
  everybeam::CorrectionMode itsMode;
  everybeam::ElementResponseModel itsElementResponseModel;

  /// If a beam had already been applied before running this step, that beam
  /// needs to undone; hence we register that beam info here:
  ///@{
  casacore::MDirection itsDirectionAtStart;
  everybeam::CorrectionMode itsModeAtStart;
  ///@}

  unsigned int itsDebugLevel;

  /// The info needed to calculate the station beams.
  ///@{
  // TODO: a copy initialization (?) somewhere is preventing a unique_ptr to
  // work
  std::vector<std::shared_ptr<everybeam::telescope::Telescope>> telescopes_;
  std::vector<size_t> ant_to_msindex_;
  std::vector<casacore::MeasFrame> itsMeasFrames;
  std::vector<casacore::MDirection::Convert> itsMeasConverters;
  std::vector<std::vector<aocommon::MC2x2>> itsBeamValues;
  ///@}

  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
