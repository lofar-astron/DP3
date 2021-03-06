// InputStep.h: Abstract base class for a Step generating input
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Abstract base class for a Step generating input
/// @author Ger van Diepen

#ifndef DPPP_DPINPUT_H
#define DPPP_DPINPUT_H

#include "Step.h"
#include "../base/DPBuffer.h"
#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"

#include <EveryBeam/station.h>
#include <EveryBeam/elementresponse.h>

#include <casacore/tables/Tables/TableIter.h>
#include <casacore/tables/Tables/RefRows.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>

#include <memory>

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
  typedef std::shared_ptr<InputStep> ShPtr;

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

  /// Read the model data at the given row numbers into the array.
  /// The default implementation throws an exception.
  virtual void getModelData(const casacore::RefRows& rowNrs,
                            casacore::Cube<casacore::Complex>&);

  /// Get the MS name.
  /// The default implementation returns an empty string.
  virtual std::string msName() const;

  /// Fill the vector with station beam info from the input source (MS).
  /// Only fill it for the given station names.
  /// The default implementation throws an exception.
  virtual void fillBeamInfo(
      std::vector<std::shared_ptr<everybeam::Station>>&,
      const casacore::Vector<casacore::String>& antNames,
      const everybeam::ElementResponseModel element_reponse_model) const;

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

  /// Creates an MS reader.
  /// Based on the MS it will create either a BDAMSReader or a regular
  static std::unique_ptr<InputStep> CreateReader(const common::ParameterSet&,
                                                 const std::string&);
};

}  // namespace steps
}  // namespace dp3

#endif
