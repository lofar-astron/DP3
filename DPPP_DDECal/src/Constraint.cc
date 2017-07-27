#ifdef AOPROJECT
#include "Constraint.h"
#include <omp.h> // for tec constraints
#else
#include <DPPP_DDECal/Constraint.h>
#include <Common/OpenMP.h>
#endif


std::vector<Constraint::Result> PhaseOnlyConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double)
{
  for (uint ch=0; ch<solutions.size(); ++ch) {
    for (uint solIndex=0; solIndex<solutions[ch].size(); ++solIndex) {
      solutions[ch][solIndex] /= std::abs(solutions[ch][solIndex]);
    }
  }

  return std::vector<Constraint::Result>();
}

std::vector<Constraint::Result> AmplitudeOnlyConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double)
{
  for (uint ch=0; ch<solutions.size(); ++ch) {
    for (uint solIndex=0; solIndex<solutions[ch].size(); ++solIndex) {
      solutions[ch][solIndex] = std::abs(solutions[ch][solIndex]);
    }
  }

  return std::vector<Constraint::Result>();
}

TECConstraint::TECConstraint(Mode mode) :
  _mode(mode),
  _nAntennas(0)
  ,_nDirections(0),
  _nChannelBlocks(0),
  _phaseFitters()
{
}

TECConstraint::TECConstraint(Mode mode, size_t nAntennas, size_t nDirections, 
                             size_t nChannelBlocks, const double* frequencies) :
  _mode(mode)
{
  init(nAntennas, nDirections, nChannelBlocks, frequencies);
}

void TECConstraint::init(size_t nAntennas, size_t nDirections, 
                         size_t nChannelBlocks, const double* frequencies) {
  _nAntennas = nAntennas;
  _nDirections = nDirections;
  _nChannelBlocks = nChannelBlocks;
  _phaseFitters.resize(
#ifdef AOPROJECT
      omp_get_max_threads()
#else
      LOFAR::OpenMP::maxThreads()
#endif
   );

  for(size_t i=0; i!=_phaseFitters.size(); ++i)
  {
    _phaseFitters[i].SetChannelCount(_nChannelBlocks);
    std::memcpy(_phaseFitters[i].FrequencyData(), frequencies, sizeof(double) * _nChannelBlocks);
  }
}

std::vector<Constraint::Result> TECConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double)
{
  std::vector<Constraint::Result> res(2);

  res[0].vals.resize(_nAntennas*_nDirections);
  res[0].axes="ant,dir,freq";
  res[0].name="tec";
  res[0].dims.resize(3);
  res[0].dims[0]=_nAntennas;
  res[0].dims[1]=_nDirections;
  res[0].dims[2]=1;
  res[1]=res[0];
  res[1].name="scalarphase";

  // TODO chose this more cleverly?
  size_t refAntenna = 0;

  // Divide out the reference antenna
  for(size_t ch=0; ch!=_nChannelBlocks; ++ch)
  {
    for(size_t antennaIndex=0; antennaIndex!=_nAntennas; ++antennaIndex)
    {
      for(size_t d=0; d!=_nDirections; ++d)
      {
        size_t solutionIndex = antennaIndex*_nDirections + d;
        size_t refAntennaIndex = d + refAntenna*_nDirections;
        if(antennaIndex != refAntenna)
        {
          solutions[ch][solutionIndex] = solutions[ch][solutionIndex] / solutions[ch][refAntennaIndex];
        }
      }
    }
    for(size_t d=0; d!=_nDirections; ++d)
      solutions[ch][refAntenna*_nDirections + d] = 1.0;
  }
  
#pragma omp parallel for
  for(size_t solutionIndex = 0; solutionIndex<_nAntennas*_nDirections; ++solutionIndex)
  {
    size_t thread =
#ifdef AOPROJECT
        omp_get_thread_num();
#else
        LOFAR::OpenMP::threadNum();
#endif

    for(size_t ch=0; ch!=_nChannelBlocks; ++ch) {
      _phaseFitters[thread].PhaseData()[ch] = std::arg(solutions[ch][solutionIndex]);
    }
    /*if(solutionIndex == _nDirections*60)
    {
      std::cout << "BEFORE ";
      for(size_t ch=0; ch!=_nChannelBlocks; ++ch) {
        std::cout << _phaseFitters[thread].PhaseData()[ch] << ' ';
      }
      std::cout << '\n';
    }*/
    
    double alpha, beta=0.0;
    if(_mode == TECOnlyMode) {
      _phaseFitters[thread].FitDataToTEC1Model(alpha);
    } else {
      _phaseFitters[thread].FitDataToTEC2Model(alpha, beta);
    }

    res[0].vals[solutionIndex] = alpha / -8.44797245e9;
    res[1].vals[solutionIndex] = beta;
    
    for(size_t ch=0; ch!=_nChannelBlocks; ++ch) 
    {
     solutions[ch][solutionIndex] = std::polar<double>(1.0, _phaseFitters[thread].PhaseData()[ch]);
    }
    /*if(solutionIndex == _nDirections*60)
    {
      std::cout << "AFTER ";
      for(size_t ch=0; ch!=_nChannelBlocks; ++ch) {
        std::cout << _phaseFitters[thread].PhaseData()[ch] << ' ';
      }
      std::cout << " (a=" << (alpha / -8.44797245e9) << ", b=" << beta << ")\n";
    }*/
  }

  return res;
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
