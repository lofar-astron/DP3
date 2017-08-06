#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#ifdef AOPROJECT
#include "PhaseFitter.h"
#define UPTR std::unique_ptr
#else
#include <DPPP/PhaseFitter.h>
#define UPTR std::auto_ptr
#endif

#include <complex>
#include <memory>
#include <set>
#include <vector>

/**
 * This class is the base class for classes that implement a constraint on
 * calibration solutions. Constraints are used to increase
 * the converge of calibration by applying these inside the solving step.
 * 
 * The MultiDirSolver class uses this class for constrained calibration.
 */
class Constraint
{
public:
  typedef std::complex<double> dcomplex;
  struct Result
  {
  public:
    std::vector<double> vals;
    std::string axes; // Comma-separated string with axes names, fastest varying last
    std::vector<size_t> dims;
    std::string name;
  };

  virtual ~Constraint() { }
   
  /**
   * TODO To be removed: constraints should have their own init specific initialization methods.
   */
  virtual void init(size_t nAntennas, size_t nDirections, 
                    size_t nChannelBlocks, const double* frequencies) = 0;

  /**
   * This method applies the constraints to the solutions.
   * @param solutions is an array of array, such that:
   * - solutions[ch] is a pointer for channelblock ch to antenna x directions solutions.
   * - directions is the dimension with the fastest changing index.
   * @param time Central time of interval.
   */
  virtual std::vector<Result> Apply(
    std::vector<std::vector<dcomplex> >& solutions,
    double time) = 0;
};

/**
 * This class constraints the amplitudes of the solution to be unity, but
 * keeps the phase.
 */
class PhaseOnlyConstraint : public Constraint
{
public:
  PhaseOnlyConstraint() {};

  virtual void init(size_t, size_t, size_t, const double*) {};

  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time);
};

/**
 * This class constraints the phases of the solution to be zero, but
 * keeps the amplitude information.
 */
class AmplitudeOnlyConstraint : public Constraint
{
public:
  AmplitudeOnlyConstraint() {};

  virtual void init(size_t, size_t, size_t, const double*) {};

  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time);
};

class TECConstraint : public Constraint
{
public:
  enum Mode {
    /** Solve for both a (differential) TEC and an XX/YY-common scalar */
    TECAndCommonScalarMode,
    /** Solve only for a (differential) TEC value */
    TECOnlyMode
  };
  
  TECConstraint(Mode mode, size_t nAntennas, size_t nDirections, 
                size_t nChannelBlocks, const double* frequencies);
  TECConstraint(Mode mode);

  virtual void init(size_t nAntennas, size_t nDirections, 
                    size_t nChannelBlocks, const double* frequencies);
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                       double time);
  
private:
  Mode _mode;
  size_t _nAntennas, _nDirections, _nChannelBlocks;
  std::vector<PhaseFitter> _phaseFitters;
};

/**
 * This constraint limits averages the solutions of a given list of antennas,
 * so that they have equal solutions.
 * 
 * The DDE solver uses this constraint to average the solutions of the core
 * antennas. Core antennas are there defined by a given maximum distance from
 * a reference antenna. In that case, the reference antenna is by default the first
 * antenna. This constraint is meant to force all core stations to
 * have the same solution, thereby decreasing the noise in their solutions.
 */
class CoreConstraint : public Constraint
{
public:
  CoreConstraint() { }

  virtual void init(size_t nAntennas, size_t nDirections, 
                    size_t nChannelBlocks, const double*)
  {
    _nAntennas = nAntennas;
    _nDirections = nDirections;
    _nChannelBlocks = nChannelBlocks;
  }
                     
  void setCoreAntennas(const std::set<size_t>& coreAntennas)
  {
    _coreAntennas = coreAntennas;
  }
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time);
  
private:
  size_t _nAntennas, _nDirections, _nChannelBlocks;
  std::set<size_t> _coreAntennas;
};

#endif
