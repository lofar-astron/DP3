// DPBuffer.h: Buffer holding the data of a timeslot/band
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Buffer holding the data of a timeslot/band
/// @author Ger van Diepen

#ifndef DPPP_DPBUFFER_H
#define DPPP_DPBUFFER_H

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Cube.h>

#include "../common/Types.h"

namespace dp3 {
namespace base {

/// @brief Buffer holding the data of a timeslot/band

/// This class holds the data for one time slot in Array variables.
/// It makes heavy use of reference semantics to avoid data copying
/// when data are pushed from one step to another.
/// This means that a data array can be shared between Step objects.
/// So if a Step object changes data in a buffer, it has to be sure
/// it can do it. If needed, Array::unique should be called to ensure
/// the array is not shared.
///
/// The following data can be kept in a DPBuffer object.
/// <table>
///  <tr>
///   <td>TIME</td>
///   <td>The time slot center of the current data (in MJD seconds).</td>
///  </tr>
///  <tr>
///   <td>ROWNRS</td>
///   <td>The row numbers of the current time slot. It can be empty
///       when e.g. a time slot is inserted or if data are averaged.</td>
///  </tr>
///  <tr>
///   <td>DATA</td>
///   <td>The visibility data as [ncorr,nchan,nbaseline].</td>
///  </tr>
///  <tr>
///   <td>FLAG</td>
///   <td>The data flags as [ncorr,nchan,nbaseline] (True is bad).
///       Note that the ncorr axis is redundant because NDPPP will always
///       have the same flag for all correlations. The reason all
///       correlations are there is because the MS expects them.</td>
///  </tr>
///  <tr>
///   <td>WEIGHT</td>
///   <td>The data weights as [ncorr,nchan,nbaseline].
///       Similarly to FLAG the ncorr axis is redundant because the
///       same weight is used for all correlations.</td>
///  </tr>
///  <tr>
///   <td>UVW</td>
///   <td>The UVW coordinates in meters as [3,nbaseline].</td>
///  </tr>
///  <tr>
///   <td>FULLRESFLAG</td>
///   <td>The flags before any averaging was done. In the MS they are kept
///       in column LOFAR_FULL_RES_FLAG. They are used to deal in BBS
///       in a smart way with bandwidth and time smearing.
///       The shape of the array is [nchan, ntimeavg, nbaseline], where
///       ntimeavg is the number of time slots averaged to a single one.
///       The number of channels averaged to a single one can be determined
///       by dividing nchan by the number of channels in the data (or flags).
///   </td>
///  </tr>
/// </table>
/// The FLAG data member should always be filled in, so the first Step
/// (MSReader) will do that. The DATA data member is filled in if any
/// Step needs DATA. Other data members are filled on demand.
/// The InputStep::fetch functions should be used to get data for those
/// members. They take care that the input buffer's data are used if
/// available, otherwise they get it from the InputStep object.
/// In that way as little memory as needed is used. Note that e.g. the
/// AOFlagger can use a lot of memory if a large time window is used.
///
/// Until early 2015 NDPPP used the strategy of shallow data copies.
/// I.e., a step increased the data reference counter and did not make
/// an actual copy. Only when data were changed, a new data array was made.
/// Thus MSReader allocated a new array when it read the data.
/// However, it appeared this strategy lead to memory fragmentation and
/// to sudden jumps in memory usage on Linux systems.
/// <br>Therefore the strategy was changed to having each step preallocate
/// its buffers and making deep copies when moving data from one step to
/// the next one. It appeared that it not only improved memory usage,
/// but also improved performance, possible due to far less mallocs.
///
/// The buffer/step guidelines are as follows:
/// 1. If a step keeps a buffer for later processing (e.g. AORFlagger),
///    it must make a copy of the buffer because the input data arrays
///    might have changed before that step processes the data.
/// 2. A shallow copy of a data member can be used if a step processes
///    the data immediately (e.g. Averager).
/// The InputStep::fetch functions come in those 2 flavours.
class DPBuffer {
 public:
  using Complex = std::complex<float>;

  /// Construct object with empty arrays.
  explicit DPBuffer(double time = 0.0, double exposure = 0.0);

  /// The copy constructor uses reference copies.
  DPBuffer(const DPBuffer&);

  /// The move constructor moves all data without using reference semantics.
  DPBuffer(DPBuffer&&);

  /// Copy assignment uses reference copies.
  DPBuffer& operator=(const DPBuffer&);

  /// Move assignment moves all data without using reference semantics.
  DPBuffer& operator=(DPBuffer&&);

  /// Make a deep copy of all arrays in that to this.
  void copy(const DPBuffer& that);

  /// Reference only the arrays that are filled in that.
  void referenceFilled(const DPBuffer& that);

  /// Set or get the visibility data per corr,chan,baseline.
  void setData(const casacore::Cube<Complex>& data) { itsData.reference(data); }
  const casacore::Cube<Complex>& getData() const { return itsData; }
  casacore::Cube<Complex>& getData() { return itsData; }

  /// Set or get the flags per corr,chan,baseline.
  void setFlags(const casacore::Cube<bool>& flags) {
    itsFlags.reference(flags);
  }
  const casacore::Cube<bool>& getFlags() const { return itsFlags; }
  casacore::Cube<bool>& getFlags() { return itsFlags; }

  /// Set or get the weights per corr,chan,baseline.
  void setWeights(const casacore::Cube<float>& weights) {
    itsWeights.reference(weights);
  }
  const casacore::Cube<float>& getWeights() const { return itsWeights; }
  casacore::Cube<float>& getWeights() { return itsWeights; }

  /// Set or get the flags at the full resolution per chan,timeavg,baseline.
  void setFullResFlags(const casacore::Cube<bool>& flags) {
    itsFullResFlags.reference(flags);
  }
  const casacore::Cube<bool>& getFullResFlags() const {
    return itsFullResFlags;
  }
  casacore::Cube<bool>& getFullResFlags() { return itsFullResFlags; }

  /// Get or set the time.
  void setTime(double time) { itsTime = time; }
  double getTime() const { return itsTime; }

  /// Get or set the exposure.
  void setExposure(double exposure) { itsExposure = exposure; }
  double getExposure() const { return itsExposure; }

  /// Get or set the row numbers used by the InputStep class.
  /// It can be empty (e.g. when MSReader inserted a dummy time slot).
  void setRowNrs(const casacore::Vector<common::rownr_t>& rownrs) {
    itsRowNrs.reference(rownrs);
  }
  const casacore::Vector<common::rownr_t>& getRowNrs() const {
    return itsRowNrs;
  }

  /// Get or set the UVW coordinates per baseline.
  void setUVW(const casacore::Matrix<double>& uvw) { itsUVW.reference(uvw); }
  const casacore::Matrix<double>& getUVW() const { return itsUVW; }
  casacore::Matrix<double>& getUVW() { return itsUVW; }

  /// Merge the flags into the pre-average flags.
  /// For each flagged point, the corresponding pre-average flags are set.
  static void mergeFullResFlags(casacore::Cube<bool>& fullResFlags,
                                const casacore::Cube<bool>& flags);

 private:
  double itsTime;
  double itsExposure;
  casacore::Vector<common::rownr_t> itsRowNrs;
  casacore::Cube<Complex> itsData;       ///< ncorr,nchan,nbasel
  casacore::Cube<bool> itsFlags;         ///< ncorr,nchan,nbasel
  casacore::Matrix<double> itsUVW;       ///< 3,nbasel
  casacore::Cube<float> itsWeights;      ///< ncorr,nchan,nbasel
  casacore::Cube<bool> itsFullResFlags;  ///< fullres_nchan,ntimeavg,nbl
};

}  // namespace base
}  // namespace dp3

#endif
