// InputStep.h: Abstract base class for a Step generating input
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Abstract base class for a Step generating input
/// @author Ger van Diepen

#ifndef DP3_INPUTSTEP_H
#define DP3_INPUTSTEP_H

#include "Step.h"
#include "../base/DPBuffer.h"
#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"

#include <EveryBeam/telescope/phasedarray.h>
#include <EveryBeam/elementresponse.h>

#include <casacore/casa/Arrays/Vector.h>

#include <memory>

namespace casacore {
class MeasurementSet;
class RefRows;
}  // namespace casacore

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// @brief Abstract base class for a Step generating input

/// This class is the abstract base class for a Step object that
/// handles the input. A concrete example is MSReader that reads the
/// data from a MeasurementSet. However, it is also possible to have
/// input steps generating data on the fly as done in test programs
/// like tAverager.cc.
///
/// A particular task of the class is to fetch the input for various
/// data items like weight, uvw, etc.. This is done by testing if the
/// item's data array is in the DPBuffer. If so, it will be returned.
/// Otherwise the appropriate 'get' function will be called to read the
/// data array from the input.
/// The derived classes should implement those 'get' functions, unless
/// they are sure the data arrays are always put in the buffer.

class InputStep : public Step {
 public:
  virtual ~InputStep();

  /// Read the UVW at the given row numbers into the buffer.
  /// The default implementation throws an exception.
  virtual void getUVW(const casacore::RefRows& rowNrs, double time,
                      base::DPBuffer&);

  /// Read the weights at the given row numbers into the buffer.
  /// The default implementation throws an exception.
  virtual void getWeights(const casacore::RefRows& rowNrs, base::DPBuffer&);

  /// Read the fullRes flags (LOFAR_FULL_RES_FLAG) at the given row numbers
  /// into the buffer.
  /// If undefined, false is returned.
  /// The default implementation throws an exception.
  virtual bool getFullResFlags(const casacore::RefRows& rowNrs,
                               base::DPBuffer&);

  /// Get the MS name.
  /// The default implementation returns an empty string.
  virtual std::string msName() const;

  /// Retrieve the everybeam telescope from the input source (MS).
  /// Default implementation throws an exception.
  virtual std::unique_ptr<everybeam::telescope::Telescope> GetTelescope(
      const everybeam::ElementResponseModel element_response_model,
      bool use_channel_frequency) const;

  /// Select station indices that match the input vector of antenna(aka station)
  /// names from a telescope pointer.
  static std::vector<size_t> SelectStationIndices(
      const everybeam::telescope::Telescope* telescope,
      const casacore::Vector<casacore::String>& station_names);

  /// Tell if the visibility data are to be read. If set to true once,
  /// this will stay true.
  virtual void setReadVisData(bool);

  /// Get the main MS table.
  const virtual casacore::Table& table() const;

  /// Get the time information.
  virtual double firstTime() const;
  virtual double lastTime() const;

  /// Get the selected spectral window.
  virtual unsigned int spectralWindow() const;

  /// Get the nr of averaged full resolution channels.
  virtual unsigned int nchanAvgFullRes() const;
  /// Get the nr of averaged full resolution time slots.
  virtual unsigned int ntimeAvgFullRes() const;

  /// Fetch the FullRes flags.
  /// If defined in the buffer, they are taken from there.
  /// Otherwise there are read from the input.
  /// If not defined in the input, they are filled using the flags in the
  /// buffer assuming that no averaging has been done so far.
  /// If desired, they can be merged with the buffer's FLAG which means
  /// that if an averaged channel is flagged, the corresponding FullRes
  /// flags are set.
  /// It does a stop/start of the timer when actually reading the data.
  const casacore::Cube<bool>& fetchFullResFlags(const base::DPBuffer& bufin,
                                                base::DPBuffer& bufout,
                                                common::NSTimer& timer,
                                                bool merge = false);

  /// Fetch the weights.
  /// If defined in the buffer, they are taken from there.
  /// Otherwise there are read from the input.
  /// If they have to be read and if autoweighting is in effect, the buffer
  /// must contain DATA to calculate the weights.
  /// <br>It does a stop/start of the timer when actually reading the data.
  const casacore::Cube<float>& fetchWeights(const base::DPBuffer& bufin,
                                            base::DPBuffer& bufout,
                                            common::NSTimer& timer);

  /// Fetch the UVW.
  /// If defined in the buffer, they are taken from there.
  /// Otherwise there are read from the input.
  /// <br>It does a stop/start of the timer when actually reading the data.
  const casacore::Matrix<double>& fetchUVW(const base::DPBuffer& bufin,
                                           base::DPBuffer& bufout,
                                           common::NSTimer& timer);

  /// Check if a measurement set contains Baseline Dependent Averaged data.
  /// @param ms A casacore measurement set.
  /// @return true if the measurement set has BDA data, false if it is regular.
  static bool HasBda(const casacore::MeasurementSet& ms);

  /// Creates an MS reader.
  /// Based on the MS it will create either a BDAMSReader or a regular
  static std::unique_ptr<InputStep> CreateReader(const common::ParameterSet&);
};

}  // namespace steps
}  // namespace dp3

#endif
