#ifndef MULTI_DIR_SOLVER_H
#define MULTI_DIR_SOLVER_H

#ifdef AOPROJECT
#include "phasefitter.h"
#include "Constraint.h"
#define UPTR std::unique_ptr
#else
#include <DPPP/phasefitter.h>
#include <DPPP_DDECal/Constraint.h>
#define UPTR std::auto_ptr
#endif

#include <armadillo>

#include <complex>
#include <vector>
#include <memory>

class MultiDirSolver
{
public:
  typedef std::complex<double> DComplex;
  typedef std::complex<float> Complex;
  
  struct SolveResult {
    size_t iterations;
    std::vector<std::vector<Constraint::Result> > _results;
  };
  
  enum CalibrationMode { CalibrateComplexGain, 
                         CalibrateTEC1, 
                         CalibrateTEC2, 
                         CalibratePhase };
  
  MultiDirSolver(size_t maxIterations, double accuracy, double stepSize, bool phaseOnly = true);
  
  void init(size_t nAntennas, size_t nDirections, size_t nChannels, 
            const std::vector<int>& ant1, const std::vector<int>& ant2);
  
  // data[i] is een pointer naar de data voor tijdstap i, vanaf die pointer staat het in volgorde als in MS (bl, chan, pol)
  // mdata[i] is een pointer voor tijdstap i naar arrays van ndir model data pointers (elk van die data pointers staat in zelfde volgorde als data)
  // solutions[ch] is een pointer voor channelblock ch naar antenna x directions oplossingen.
  SolveResult process(std::vector<Complex*>& data, std::vector<std::vector<Complex* > >& modelData,
    std::vector<std::vector<DComplex> >& solutions, double time) const;
  
  void set_mode(CalibrationMode mode) { _mode = mode; }
  
  void set_channel_blocks(size_t nChannelBlocks) { _nChannelBlocks = nChannelBlocks; }
  
  void set_max_iterations(size_t maxIterations) { _maxIterations = maxIterations; }
  
  void set_accuracy(double accuracy) { _accuracy = accuracy; }
  
  void set_step_size(double stepSize) { _stepSize = stepSize; }
  
  void add_constraint(Constraint* constraint) { _constraints.push_back(constraint); }
  
private:
  void performSolveIteration(size_t channelBlockIndex,
                             std::vector<arma::cx_mat>& gTimesCs,
                             std::vector<arma::cx_vec>& vs,
                             const std::vector<DComplex>& solutions,
                             std::vector<DComplex>& nextSolutions,
                             const std::vector<Complex *>& data,
                             const std::vector<std::vector<Complex *> >& modelData) const;
  
  size_t _nAntennas, _nDirections, _nChannels, _nChannelBlocks;
  std::vector<int> _ant1, _ant2;
  
  // Calibration setup
  enum CalibrationMode _mode;
  size_t _maxIterations;
  double _accuracy;
  double _stepSize;
  bool _phaseOnly;
  std::vector<Constraint*> _constraints;
};

#endif
