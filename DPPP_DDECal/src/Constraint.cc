#ifdef AOPROJECT
#include "Constraint.h"
#else
#include <DPPP_DDECal/Constraint.h>
#endif

std::vector<Constraint::Result> PhaseOnlyConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double)
{
  for (size_t ch=0; ch<solutions.size(); ++ch) {
    for (size_t solIndex=0; solIndex<solutions[ch].size(); ++solIndex) {
      solutions[ch][solIndex] /= std::abs(solutions[ch][solIndex]);
    }
  }

  return std::vector<Constraint::Result>();
}

std::vector<Constraint::Result> AmplitudeOnlyConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double)
{
  for (size_t ch=0; ch<solutions.size(); ++ch) {
    for (size_t solIndex=0; solIndex<solutions[ch].size(); ++solIndex) {
      solutions[ch][solIndex] = std::abs(solutions[ch][solIndex]);
    }
  }

  return std::vector<Constraint::Result>();
}

std::vector<Constraint::Result> DiagonalConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double)
{
  if(_polsPerSolution == 4)
  {
    for (size_t ch=0; ch<solutions.size(); ++ch) {
      for (size_t solIndex=0; solIndex<solutions[ch].size(); solIndex += 4) {
        solutions[ch][solIndex+1] = 0.0;
        solutions[ch][solIndex+2] = 0.0;
      }
    }
  }

  return std::vector<Constraint::Result>();
}

std::vector<Constraint::Result> CoreConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double)
{
  for (uint ch=0; ch<solutions.size(); ++ch) {
    std::vector<dcomplex> coreSolutions(_nDirections, 0.0);
    // Calculate the sum of solutions over the core stations
    for(std::set<size_t>::const_iterator antennaIter = _coreAntennas.begin(); antennaIter!=_coreAntennas.end(); ++antennaIter)
    {
      size_t startIndex = (*antennaIter)*_nDirections;
      for(size_t direction = 0; direction != _nDirections; ++direction)
        coreSolutions[direction] += solutions[ch][startIndex + direction];
    }
    
    // Divide by nr of core stations to get the mean solution
    for(std::vector<dcomplex>::iterator solutionIter = coreSolutions.begin(); solutionIter!=coreSolutions.end(); ++solutionIter)
      (*solutionIter) /= _coreAntennas.size();
    
    // Assign all core stations to the mean solution
    for(std::set<size_t>::const_iterator antennaIter = _coreAntennas.begin(); antennaIter!=_coreAntennas.end(); ++antennaIter)
    {
      size_t startIndex = (*antennaIter)*_nDirections;
      for(size_t direction = 0; direction != _nDirections; ++direction)
        solutions[ch][startIndex + direction] = coreSolutions[direction];
    }
  }
  return std::vector<Constraint::Result>();
}
