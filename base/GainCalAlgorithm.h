// GainCalAlgorithm.h: Perform algorithm for gain calibration
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to apply a calibration correction to the data
/// @author Tammo Jan Dijkema

#ifndef DPPP_GAINCALALGORITHM_H
#define DPPP_GAINCALALGORITHM_H

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>

namespace dp3 {

namespace base {
/// @brief DPPP step class to apply a calibration correction to the data

class GainCalAlgorithm {
 public:
  using DComplex = std::complex<double>;
  enum Status { CONVERGED = 1, NOTCONVERGED = 2, STALLED = 3, FAILED = 4 };

  enum Mode { DEFAULT, PHASEONLY, AMPLITUDEONLY, FULLJONES };

  /// mode can be "diagonal", "fulljones", "phaseonly", "scalarphase"
  GainCalAlgorithm(unsigned int solInt, unsigned int nChan, Mode mode,
                   bool scalar, double tolerance, unsigned int maxAntennas,
                   bool detectStalling, unsigned int debugLevel);

  /// Sets visibility matrices to zero
  void resetVis();

  /// Initializes a new run of gaincal, resizes all internal vectors
  /// If initSolutions is false, you are responsible for setting them
  /// before running the solver. You could set the solutions to those
  /// of the previous time step.
  void init(bool initSolutions);

  /// Perform an iteration of gaincal. Returns CONVERGED, NOTCONVERGED
  /// or STALLED
  Status doStep(unsigned int iter);

  /// Returns the solution. The return matrix has a length of maxAntennas,
  /// which is zero for antennas for which no solution was computed.
  /// The mapping is stored in the antenna map
  casacore::Matrix<std::complex<double>> getSolution(bool setNaNs);

  double getWeight() { return _totalWeight; }

  /// Increments the weight (only relevant for TEC-fitting)
  void incrementWeight(float weight);

  /// Returns a reference to the visibility matrix
  casacore::Array<std::complex<double>>& getVis() { return _vis; }

  /// Returns a reference to the model visibility matrix
  casacore::Array<std::complex<double>>& getMVis() { return _mvis; }

  casacore::Vector<bool>& getStationFlagged() { return _stationFlagged; }

  /// Number of correlations in the solution (1,2 or 4)
  unsigned int numCorrelations() { return _savedNCr; }

  /// Number of correlations (1 or 4)
  unsigned int nCr() { return _nCr; }

  /// Clear antFlagged
  void clearStationFlagged();

 private:
  /// Perform relaxation
  Status relax(unsigned int iter);
  static bool isFinite(const DComplex& val) {
    return std::isfinite(val.real()) && std::isfinite(val.imag());
  }

  void doStep_polarized();
  void doStep_unpolarized();

  double getAverageUnflaggedSolution();

  unsigned int _savedNCr;
  casacore::Vector<bool>
      _stationFlagged;  ///< Contains true for totally flagged stations
  casacore::Array<std::complex<double>> _vis;   ///< Visibility matrix
  casacore::Array<std::complex<double>> _mvis;  ///< Model visibility matrix
  casacore::Matrix<std::complex<double>>
      _g;  ///< Solution, indexed by station, correlation
  casacore::Matrix<std::complex<double>> _gx;  ///< Previous solution
  casacore::Matrix<std::complex<double>>
      _gxx;  ///< Solution before previous solution
  casacore::Matrix<std::complex<double>> _gold;  ///< Previous solution
  casacore::Matrix<std::complex<double>>
      _h;  ///< Hermitian transpose of previous solution
  casacore::Matrix<std::complex<double>> _z;  ///< Internal algorithm vector

  unsigned int _nSt;       ///< number of stations in the current solution
  unsigned int _nUn;       ///< number of unknowns
  unsigned int _nCr;       ///< number of correlations (1 or 4)
  unsigned int _nSp;       ///< number that is two for scalarphase, one else
  unsigned int _badIters;  ///< number of bad iterations, for stalling detection
  unsigned int
      _veryBadIters;     ///< number of iterations where solution got worse
  unsigned int _solInt;  ///< solution interval
  unsigned int _nChan;   ///< number of channels
  Mode _mode;            ///< diagonal, scalarphase, fulljones or phaseonly
  bool _scalar;          ///< false if each polarization has a separate solution
  double _tolerance;
  double _totalWeight;
  bool _detectStalling;
  unsigned int _debugLevel;

  double _dg, _dgx;          ///< previous convergence
  std::vector<double> _dgs;  ///< convergence history
};

}  // namespace base
}  // namespace dp3

#endif
