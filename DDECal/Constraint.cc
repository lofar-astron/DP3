// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

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
  // nSols is nPol x nDirections (i.e., nr of sols per antenna)
  size_t nSols = solutions.front().size() / _nAntennas;
  std::vector<dcomplex> setSolutions(nSols);
  std::vector<size_t> setSolutionCounts(nSols);
  for (unsigned int ch=0; ch<solutions.size(); ++ch) {
    for(const std::set<size_t>& antennaSet : _antennaSets)
    {
      setSolutions.assign(nSols, 0.0);
      setSolutionCounts.assign(nSols, 0);
      // Calculate the sum of solutions over the set of stations
      for(size_t antennaIndex : antennaSet)
      {
        size_t startIndex = antennaIndex * nSols;
        for(size_t solIndex = 0; solIndex != nSols; ++solIndex)
        {
          dcomplex value = solutions[ch][startIndex + solIndex];
          if(isfinite(value))
          {
            setSolutions[solIndex] += value;
            ++setSolutionCounts[solIndex];
          }
        }
      }
      
      // Divide by nr of core stations to get the mean solution
      for(size_t solIndex = 0; solIndex != nSols; ++solIndex)
        setSolutions[solIndex] /= setSolutionCounts[solIndex];
      
      // Assign all core stations to the mean solution
      for(size_t antennaIndex : antennaSet)
      {
        size_t startIndex = antennaIndex * nSols;
        for(size_t solIndex = 0; solIndex != nSols; ++solIndex)
          solutions[ch][startIndex + solIndex] = setSolutions[solIndex];
      }
    }
  }
  return std::vector<Constraint::Result>();
}
