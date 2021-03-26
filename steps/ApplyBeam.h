// ApplyBeam.h: DPPP step class to ApplyBeam visibilities from a source model
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to apply the beam model (optionally inverted)
/// @author Tammo Jan Dijkema

#ifndef DPPP_APPLYBEAM_H
#define DPPP_APPLYBEAM_H

#include "InputStep.h"

#include "../base/DPBuffer.h"
#include "../base/Position.h"

#include <EveryBeam/station.h>
#include <EveryBeam/common/types.h>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/measures/Measures/MDirection.h>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// \brief DPPP step class to ApplyBeam visibilities from a source model

/// This class is a Step class to apply the beam model, optionally inverted.
/// The input MeasurementSet it operates on, must have the LOFAR subtables
/// defining the station layout and tiles/dipoles used.

class ApplyBeam : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  ApplyBeam(InputStep*, const common::ParameterSet&, const string& prefix,
            bool substep = false);

  ApplyBeam();

  virtual ~ApplyBeam();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer& buffer) {
    return processMultithreaded(buffer, 0);
  }

  /// If apply beam is called from multiple threads, it needs the thread index
  /// to determine what scratch space to use etc.
  bool processMultithreaded(const base::DPBuffer&, size_t thread);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  bool invert() { return itsInvert; }

  template <typename T>
  static void applyBeam(
      const base::DPInfo& info, double time, T* data0, float* weight0,
      const everybeam::vector3r_t& srcdir, const everybeam::vector3r_t& refdir,
      const everybeam::vector3r_t& tiledir,
      const std::vector<std::shared_ptr<everybeam::Station>>& antBeamInfo,
      std::vector<everybeam::matrix22c_t>& beamValues, bool useChannelFreq,
      bool invert, base::BeamCorrectionMode mode,
      // everybeam::ElementResponseModel element_reponse_model,
      bool doUpdateWeights = false);

  template <typename T>
  static void applyBeamStokesIArrayFactor(
      const base::DPInfo& info, double time, T* data0, float* weight0,
      const everybeam::vector3r_t& srcdir, const everybeam::vector3r_t& refdir,
      const everybeam::vector3r_t& tiledir,
      const std::vector<std::shared_ptr<everybeam::Station>>& antBeamInfo,
      std::vector<everybeam::complex_t>& beamValues, bool useChannelFreq,
      bool invert, base::BeamCorrectionMode mode,
      // everybeam::ElementResponseModel element_reponse_model,
      bool doUpdateWeights = false);

 private:
  everybeam::vector3r_t dir2Itrf(const casacore::MDirection& dir,
                                 casacore::MDirection::Convert& measConverter);

  InputStep* itsInput;
  string itsName;
  base::DPBuffer itsBuffer;
  bool itsInvert;
  bool itsUpdateWeights;
  std::vector<string> itsDirectionStr;
  casacore::MDirection itsDirection;
  bool itsUseChannelFreq;
  base::BeamCorrectionMode itsMode;
  everybeam::ElementResponseModel itsElementResponseModel;

  /// If a beam had already been applied before running this step, that beam
  /// needs to undone; hence we register that beam info here:
  ///@{
  casacore::MDirection itsDirectionAtStart;
  base::BeamCorrectionMode itsModeAtStart;
  ///@}

  unsigned int itsDebugLevel;

  /// The info needed to calculate the station beams.
  ///@{
  std::vector<std::vector<std::shared_ptr<everybeam::Station>>> itsAntBeamInfo;
  std::vector<casacore::MeasFrame> itsMeasFrames;
  std::vector<casacore::MDirection::Convert> itsMeasConverters;
  std::vector<std::vector<everybeam::matrix22c_t>> itsBeamValues;
  ///@}

  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif  // HAVE_LOFAR_BEAM
