#ifdef AOPROJECT
#include "KernelSmoother.h"
#include "SmoothnessConstraint.h"
#include "omptools.h"
#else
#include <DPPP_DDECal/KernelSmoother.h>
#include <DPPP_DDECal/SmoothnessConstraint.h>
#include <Common/OpenMP.h>
#endif

SmoothnessConstraint::SmoothnessConstraint(double bandwidthHz) :
  _kernelType(Smoother::GaussianKernel),
  _bandwidth(bandwidthHz)
{ }

void SmoothnessConstraint::Initialize(const double* frequencies, size_t n)
{
  _frequencies.assign(frequencies, frequencies+n);
}

void SmoothnessConstraint::InitializeDimensions(size_t nAntennas,
  size_t nDirections,
  size_t nChannelBlocks)
{
  Constraint::InitializeDimensions(nAntennas, nDirections, nChannelBlocks);
  size_t nthreads =
 #ifdef AOPROJECT
    omp_get_max_threads();
#else
    LOFAR::OpenMP::maxThreads();
#endif
  for(size_t i=0; i!=nthreads; ++i)
    _fitData.emplace_back(_frequencies.data(), _frequencies.size(), _kernelType, _bandwidth);
}

std::vector<Constraint::Result> SmoothnessConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double, std::ostream*)
{
  std::vector<dcomplex> data(solutions.size());
  size_t thread =
#ifdef AOPROJECT
      omp_get_thread_num();
#else
      LOFAR::OpenMP::threadNum();
#endif
  
#pragma omp parallel for
  for(size_t solutionIndex = 0; solutionIndex<_nAntennas*_nDirections; ++solutionIndex)
  {
    size_t antIndex = solutionIndex / _nDirections;
    for(size_t ch=0; ch!=_nChannelBlocks; ++ch)
    {
      // Flag channels where calibration yielded inf or nan
      if(std::isfinite(solutions[ch][solutionIndex].real()) &&
        std::isfinite(solutions[ch][solutionIndex].imag()))
      {
        _fitData[thread].data[ch] = solutions[ch][solutionIndex];
        _fitData[thread].weight[ch] = _weights[antIndex];
      }
      else {
        _fitData[thread].data[ch] = 0.0;
        _fitData[thread].weight[ch] = 0.0;
      }
    }
    
    _fitData[thread].smoother.Smooth(_fitData[thread].data.data(), _fitData[thread].weight.data());
    
    for(size_t ch=0; ch!=_nChannelBlocks; ++ch)
    {
      solutions[ch][solutionIndex] = _fitData[thread].data[ch];
    }
  }
  
  return std::vector<Constraint::Result>();
}
