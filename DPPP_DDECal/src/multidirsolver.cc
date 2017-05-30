#ifdef AOPROJECT
#include "multidirsolver.h"
#else
#include <DPPP_DDECal/multidirsolver.h>
#endif

using namespace arma;

MultiDirSolver::MultiDirSolver(size_t maxIterations, double accuracy, double stepSize) :
  _nAntennas(0),
  _nDirections(0),
  _nChannels(0),
  _nChannelBlocks(0),
  _mode(CalibrateComplexGain),
  _maxIterations(maxIterations),
  _accuracy(accuracy),
  _stepSize(stepSize)
{
}

void MultiDirSolver::init(size_t nAntennas,
                          size_t nDirections,
                          size_t nChannels,
                          const std::vector<int>& ant1,
                          const std::vector<int>& ant2)
{
  _nAntennas = nAntennas;
  _nDirections = nDirections;
  _nChannels = nChannels;
  _nChannelBlocks = nChannels;
  _ant1 = ant1;
  _ant2 = ant2;
}

MultiDirSolver::SolveResult MultiDirSolver::process(std::vector<Complex *>& data,
  std::vector<std::vector<Complex *> >& modelData,
  std::vector<std::vector<DComplex> >& solutions, double time) const
{
  const size_t nTimes = data.size();
  SolveResult result;
  
  std::vector<std::vector<DComplex> > nextSolutions(_nChannelBlocks);
  if (solutions.size() != _nChannelBlocks) {
    cout << "Error: 'solutions' parameter does not have the right shape" << endl;
    result.iterations = 0;
    return result;
  }

  result._results.resize(_constraints.size());
  
  // Model matrix ant x [N x D] and visibility matrix ant x [N x 1],
  // for each channelblock
  // The following loop allocates all structures
  std::vector<std::vector<cx_mat> > gTimesCs(_nChannelBlocks);
  std::vector<std::vector<cx_vec> > vs(_nChannelBlocks);
  for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
  {
    //solutions[chBlock].assign(_nDirections * _nAntennas, 1.0);
    nextSolutions[chBlock].resize(_nDirections * _nAntennas);
    const size_t
      channelIndexStart = chBlock * _nChannels / _nChannelBlocks,
      channelIndexEnd = (chBlock+1) * _nChannels / _nChannelBlocks,
      curChannelBlockSize = channelIndexEnd - channelIndexStart;
    gTimesCs[chBlock].resize(_nAntennas);
    vs[chBlock].resize(_nAntennas);
    
    for(size_t ant=0; ant!=_nAntennas; ++ant)
    {
      // Model matrix [N x D] and visibility matrix x [N x 1]
      // Also space for the auto correlation is reserved, but it will be set to 0.
      gTimesCs[chBlock][ant] = cx_mat(_nAntennas * nTimes * curChannelBlockSize,
                                      _nDirections, fill::zeros);
      vs[chBlock][ant] = cx_vec(_nAntennas * nTimes * curChannelBlockSize, fill::zeros);
    }
  }
  
  // TODO the data and model data needs to be preweighted.
  // Maybe we can get a non-const pointer from DPPP, that saves copying/allocating
  
  ///
  /// Start iterating
  ///
  size_t iteration = 0;
  double normSum = 0.0, sum = 0.0;
  do {
#pragma omp parallel for
    for(size_t chBlock=0; chBlock<_nChannelBlocks; ++chBlock)
    {
      performSolveIteration(chBlock, gTimesCs[chBlock], vs[chBlock],
                            solutions[chBlock], nextSolutions[chBlock],
                            data, modelData);
    }
      
    // Move the solutions towards nextSolutions
    // (the moved solutions are stored in 'nextSolutions')
    for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
    {
      for(size_t i=0; i!=_nAntennas*_nDirections; ++i)
      {
        nextSolutions[chBlock][i] = solutions[chBlock][i]*(1.0-_stepSize) +
          nextSolutions[chBlock][i] * _stepSize;
      }
    }
    for(size_t i=0; i!=_constraints.size(); ++i)
    {
      result._results[i] = _constraints[i]->Apply(nextSolutions, time);
    }
    
    //  Calculate the norm of the difference between the old and new solutions
    for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
    {
      for(size_t i=0; i!=_nAntennas*_nDirections; ++i)
      {
        double e = std::norm(nextSolutions[chBlock][i] - solutions[chBlock][i]);
        normSum += e;
        sum += std::abs(solutions[chBlock][i]);
        
        solutions[chBlock][i] = nextSolutions[chBlock][i];
        
        // For debug: output the solutions of the first antenna
        if(i<_nDirections && false)
        {
          std::cout << " |s_" << i << "|=|" << solutions[chBlock][i] << "|="
                    << std::abs(solutions[chBlock][i]);
        }
      }
    }
    normSum /= _nChannelBlocks*_nAntennas*_nDirections;
    sum /= _nChannelBlocks*_nAntennas*_nDirections;
    iteration++;
    
  } while(iteration < _maxIterations && normSum/sum > _accuracy);
  
  if(normSum/sum <= _accuracy)
    result.iterations = iteration;
  else
    result.iterations = _maxIterations+1;
  return result;
}

void MultiDirSolver::performSolveIteration(size_t channelBlockIndex,
                       std::vector<arma::cx_mat>& gTimesCs,
                       std::vector<arma::cx_vec>& vs,
                       const std::vector<DComplex>& solutions,
                       std::vector<DComplex>& nextSolutions,
                       const std::vector<Complex *>& data,
                       const std::vector<std::vector<Complex *> >& modelData) const
{
  const size_t
    channelIndexStart = channelBlockIndex * _nChannels / _nChannelBlocks,
    channelIndexEnd = (channelBlockIndex+1) * _nChannels / _nChannelBlocks,
    curChannelBlockSize = channelIndexEnd - channelIndexStart,
    nTimes = data.size();
  
  // The following loop fills the matrices for all antennas
  for(size_t timeIndex=0; timeIndex!=nTimes; ++timeIndex)
  {
    std::vector<const Complex*> modelPtrs(_nDirections);
    for(size_t baseline=0; baseline!=_ant1.size(); ++baseline)
    {
      size_t antenna1 = _ant1[baseline];
      size_t antenna2 = _ant2[baseline];
      if(antenna1 != antenna2)
      {
        cx_mat& gTimesC1 = gTimesCs[antenna1];
        cx_vec& v1 = vs[antenna1];
        cx_mat& gTimesC2 = gTimesCs[antenna2];
        cx_vec& v2 = vs[antenna2];
        for(size_t d=0; d!=_nDirections; ++d)
          modelPtrs[d] = modelData[timeIndex][d] + (channelIndexStart + baseline * _nChannels) * 4;
        const Complex* dataPtr = data[timeIndex] + (channelIndexStart + baseline * _nChannels) * 4;
        for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
        {
          size_t dataIndex1 = ch-channelIndexStart + (timeIndex + antenna1 * nTimes) * curChannelBlockSize;
          size_t dataIndex2 = ch-channelIndexStart + (timeIndex + antenna2 * nTimes) * curChannelBlockSize;
          //std::cout << "timeindex" << timeIndex << ' ';
          for(size_t d=0; d!=_nDirections; ++d)
          {
            std::complex<double> predicted = modelPtrs[d][0] + modelPtrs[d][3];
            //std::cout << predicted << ' ';
            
            size_t solIndex1 = antenna1*_nDirections + d;
            size_t solIndex2 = antenna2*_nDirections + d;
            gTimesC1(dataIndex2, d) = std::conj(solutions[solIndex2]) * std::conj(predicted);
	    gTimesC2(dataIndex1, d) = std::conj(solutions[solIndex1]) * predicted;
            
            modelPtrs[d] += 4; // Goto the next 2x2 matrix.
          }
          v2(dataIndex1) = dataPtr[0] + dataPtr[3]; // Solve using Stokes I
          v1(dataIndex2) = std::conj(v2(dataIndex1));
          dataPtr += 4; // Goto the next 2x2 matrix.
        }
      }
    }
  }
  
  // The matrices have been filled; compute the linear solution
  // for each antenna.
  for(size_t ant=0; ant!=_nAntennas; ++ant)
  {
    //std::cout << ant << '\n';
    cx_mat& gTimesC = gTimesCs[ant];
    cx_vec& v = vs[ant];
    // solve [g* C] x  = v
    cx_vec x = solve(gTimesC, v);
    for(size_t d=0; d!=_nDirections; ++d)
      nextSolutions[ant*_nDirections + d] = x(d);
  }
}
