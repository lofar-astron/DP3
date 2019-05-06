#include "Constraint.h"

std::vector<Constraint::Result> PhaseOnlyConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double,
    std::ostream* statStream)
{
  for (size_t ch=0; ch<solutions.size(); ++ch) {
    for (size_t solIndex=0; solIndex<solutions[ch].size(); ++solIndex) {
      solutions[ch][solIndex] /= std::abs(solutions[ch][solIndex]);
    }
  }

  return std::vector<Constraint::Result>();
}

std::vector<Constraint::Result> AmplitudeOnlyConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double,
    std::ostream* statStream)
{
  for (size_t ch=0; ch<solutions.size(); ++ch) {
    for (size_t solIndex=0; solIndex<solutions[ch].size(); ++solIndex) {
      solutions[ch][solIndex] = std::abs(solutions[ch][solIndex]);
    }
  }

  return std::vector<Constraint::Result>();
}

std::vector<Constraint::Result> DiagonalConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double,
    std::ostream* statStream)
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

std::vector<Constraint::Result> AntennaConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double,
    std::ostream* statStream)
{
  std::vector<dcomplex> setSolutions(_nDirections);
  for (unsigned int ch=0; ch<solutions.size(); ++ch) {
    for(const std::set<size_t>& antennaSet : _antennaSets)
    {
      setSolutions.assign(_nDirections, 0.0);
      // Calculate the sum of solutions over the core stations
      for(size_t antennaIndex : antennaSet)
      {
        size_t startIndex = antennaIndex * _nDirections;
        for(size_t direction = 0; direction != _nDirections; ++direction)
          setSolutions[direction] += solutions[ch][startIndex + direction];
      }
      
      // Divide by nr of core stations to get the mean solution
      for(dcomplex& solution : setSolutions)
        solution /= antennaSet.size();
      
      // Assign all core stations to the mean solution
      for(size_t antennaIndex : antennaSet)
      {
        size_t startIndex = antennaIndex * _nDirections;
        for(size_t direction = 0; direction != _nDirections; ++direction)
          solutions[ch][startIndex + direction] = setSolutions[direction];
      }
    }
  }
  return std::vector<Constraint::Result>();
}
