#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <complex>
#include <memory>
#include <set>
#include <vector>
#include <ostream>

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
    std::vector<double> weights;
    std::string axes; // Comma-separated string with axis names, fastest varying last
    std::vector<size_t> dims;
    std::string name;
  };

  virtual ~Constraint() { }

  /**
   * Function that initializes the constraint for the next calibration iteration.
   * It should be called each time all antenna solutions have been calculated,
   * but before the constraint has been applied to all those antenna solutions.
   * 
   * Unlike Apply(), this method is not thread safe.
   * 
   * @param bool This can be used to specify whether the previous solution "step" is
   * smaller than the requested precision, i.e. calibration with the constrained
   * has converged. This allows a constraint to apply
   * its constraint in steps: apply a better-converging constraint as long as the
   * solutions are far from the correct answer, then switch to a different constraint
   * when hasReachedPrecision=true.
   */
  virtual void PrepareIteration(bool /*hasReachedPrecision*/, size_t /*iteration*/, bool /*finalIter*/) { }
  
  /**
   * Whether the constraint has been satisfied. The calibration process will continue
   * at least as long as Satisfied()=false, and performs at least one more iteration
   * after Satisfied()=true. Together with SetPrecisionReached(), this
   * can make the algorithm change the constraining method based on amount of
   * convergence.
   */
  virtual bool Satisfied() const { return true; }
   
  /**
   * This method applies the constraints to the solutions.
   * @param solutions is an array of array, such that:
   * - solutions[ch] is a pointer for channelblock ch to antenna x directions x pol solutions.
   * - pol is the dimension with the fastest changing index.
   * @param time Central time of interval.
   */
  virtual std::vector<Result> Apply(
    std::vector<std::vector<dcomplex> >& solutions,
    double time, std::ostream* statStream) = 0;

  /**
  * Initialize the dimensions for the constraint. Should be overridden when
  * something more than assigning dimensions is needed (e.g. resizing vectors).
  * Weights are initialized to 1. here.
  */
  virtual void InitializeDimensions(size_t nAntennas,
                                    size_t nDirections,
                                    size_t nChannelBlocks)
  {
    _nAntennas = nAntennas;
    _nDirections = nDirections;
    _nChannelBlocks = nChannelBlocks;
  }

  /**
   * Set weights. The vector should contain an array of size nAntennas * nChannelBlocks,
   * where the channel index varies fastest.
   */
  virtual void SetWeights(const std::vector<double> &) {}

  virtual void showTimings (std::ostream&, double) const {}

protected:
  size_t _nAntennas, _nDirections, _nChannelBlocks;
};

/**
 * This class constraints the amplitudes of the solution to be unity, but
 * keeps the phase.
 */
class PhaseOnlyConstraint : public Constraint
{
public:
  PhaseOnlyConstraint() {};

  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time,
                    std::ostream* statStream);
};

/**
 * This class constraints the phases of the solution to be zero, but
 * keeps the amplitude information.
 */
class AmplitudeOnlyConstraint : public Constraint
{
public:
  AmplitudeOnlyConstraint() {};

  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time,
                    std::ostream* statStream);
};

class DiagonalConstraint : public Constraint
{
public:
  DiagonalConstraint(size_t polsPerSolution) : _polsPerSolution(polsPerSolution) {};
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time,
                    std::ostream* statStream);
private:
  const size_t _polsPerSolution;
};

/**
 * This constraint averages the solutions of a given list of antennas,
 * so that they have equal solutions.
 * 
 * The DDE solver uses this constraint to average the solutions of the core
 * antennas. Core antennas are determined by a given maximum distance from
 * a reference antenna. The reference antenna is by default the first
 * antenna. This constraint is meant to force all core stations to
 * have the same solution, thereby decreasing the noise in their solutions.
 */
class CoreConstraint : public Constraint
{
public:
  CoreConstraint() { }

  void initialize(const std::set<size_t>& coreAntennas)
  {
    _coreAntennas = coreAntennas;
  }
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time,
                    std::ostream* statStream);
  
private:
  std::set<size_t> _coreAntennas;
};

#endif
