
#include "MultiDirSolver.h"
#include "Matrix2x2.h"
#include "QRSolver.h"

#ifdef AOPROJECT
#include "ParallelFor.h"
#else
#include "../Common/ParallelFor.h"
#endif

#include <iomanip>
#include <iostream>

MultiDirSolver::MultiDirSolver() :
  _nAntennas(0),
  _nDirections(0),
  _nChannels(0),
  _nChannelBlocks(0),
  _maxIterations(100),
  _accuracy(1e-5),
  _constraintAccuracy(1e-4),
  _stepSize(0.2),
  _detectStalling(true),
  _phaseOnly(false)
{ }

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
  _buffer.SetDimensions(nDirections, nChannels, ant1.size());
}

void MultiDirSolver::makeStep(const std::vector<std::vector<DComplex> >& solutions,
                              std::vector<std::vector<DComplex> >& nextSolutions) const
{
  // Move the solutions towards nextSolutions
  // (the moved solutions are stored in 'nextSolutions')
  DP3::ParallelFor<size_t> loop(_nThreads);
  loop.Run(0, _nChannelBlocks, [&](size_t chBlock, size_t /*thread*/)
  {
    for(size_t i=0; i!=nextSolutions[chBlock].size(); ++i)
    {
      if(_phaseOnly)
      {
        // In phase only mode, a step is made along the complex circle,
        // towards the shortest direction.
        double phaseFrom = std::arg(solutions[chBlock][i]);
        double distance = std::arg(nextSolutions[chBlock][i]) - phaseFrom;
        if(distance > M_PI) distance = distance - 2.0*M_PI;
        else if(distance < -M_PI) distance = distance + 2.0*M_PI;
        nextSolutions[chBlock][i] = std::polar(1.0, phaseFrom + _stepSize * distance);
      }
      else {
        nextSolutions[chBlock][i] = solutions[chBlock][i]*(1.0-_stepSize) +
          nextSolutions[chBlock][i] * _stepSize;
      }
    }
  });
}

void MultiDirSolver::makeSolutionsFinite(std::vector<std::vector<DComplex> >& solutions, size_t perPol) const
{
  for(std::vector<DComplex>& solVector : solutions)
  {
    size_t n = solVector.size() / perPol;
    std::vector<DComplex>::iterator iter = solVector.begin();
    for(size_t i=0; i!=n; ++i)
    {
      bool hasNonFinite = false;
      for(size_t p=0; p!=perPol; ++p)
      {
        hasNonFinite = hasNonFinite || !std::isfinite(iter->real()) || !std::isfinite(iter->imag());
      }
      if(hasNonFinite)
      {
        if(perPol == 4)
        {
          iter[0] = DComplex(1.0, 0.0);
          iter[1] = DComplex(0.0, 0.0);
          iter[2] = DComplex(0.0, 0.0);
          iter[3] = DComplex(1.0, 0.0);
        }
        else {
          for(size_t p=0; p!=perPol; ++p)
          {
            iter[p] = DComplex(1.0, 0.0);
          }
        }
      }
      iter += perPol;
    }
  }
}

template<size_t NPol>
bool MultiDirSolver::assignSolutions(std::vector<std::vector<DComplex> >& solutions,
  std::vector<std::vector<DComplex> >& nextSolutions, bool useConstraintAccuracy,
  double& avgAbsDiff, std::vector<double>& stepMagnitudes) const
{
  avgAbsDiff = 0.0;
  //  Calculate the norm of the difference between the old and new solutions
  size_t n = 0;
  for(size_t chBlock=0; chBlock<_nChannelBlocks; ++chBlock)
  {
    for(size_t i=0; i!=solutions[chBlock].size(); i += NPol)
    {
      // A normalized squared difference is calculated between the solutions of this
      // and the previous step:
      //   sqrt { 1/n sum over | (t1 - t0) t0^(-1) |_2 }
      // This criterion is scale independent: all solutions can be scaled without
      // affecting the number of iterations. Also, when the polarized version is given
      // scalar matrices, it will use the same number of iterations as the scalar
      // version.
      if(NPol == 1)
      {
        if(solutions[chBlock][i] != 0.0)
        {
          double a = std::abs((nextSolutions[chBlock][i] - solutions[chBlock][i]) / solutions[chBlock][i]);
          if(std::isfinite(a))
          {
            avgAbsDiff += a;
            ++n;
          }
        }
        solutions[chBlock][i] = nextSolutions[chBlock][i];
      }
      else {
        MC2x2 s(&solutions[chBlock][i]), sInv(s);
        if(sInv.Invert())
        {
          MC2x2 ns(&nextSolutions[chBlock][i]);
          ns -= s;
          ns *= sInv;
          double sumabs = 0.0;
          for(size_t p=0; p!=NPol; ++p)
          {
            sumabs += std::abs(ns[p]);
          }
          if(std::isfinite(sumabs))
          {
            avgAbsDiff += sumabs;
            n += 4;
          }
        }
        for(size_t p=0; p!=NPol; ++p)
        {
          solutions[chBlock][i+p] = nextSolutions[chBlock][i+p];
        }
      }
    }
  }
  // The polarized version needs a factor of two normalization to make it work
  // like the scalar version would and when given only scalar matrices.
  //if(NPol == 4)
  //  avgSquaredDiff = sqrt(avgSquaredDiff*0.5/n) ;
  //else
  //  avgSquaredDiff = sqrt(avgSquaredDiff/n);

  // The stepsize is taken out, so that a small stepsize won't cause
  // a premature stopping criterion.
  double stepMagnitude = (n==0 ? 0 : avgAbsDiff/_stepSize/n);
  stepMagnitudes.emplace_back(stepMagnitude);

  if(useConstraintAccuracy)
    return stepMagnitude <= _constraintAccuracy;
  else {
    return stepMagnitude <= _accuracy;
  }
}

MultiDirSolver::SolveResult MultiDirSolver::processScalar(
  const std::vector<Complex *>& dataNoW,
  const std::vector<float*>& weights,
  const std::vector<std::vector<Complex *> >& modelDataNoW,
  std::vector<std::vector<DComplex> >& solutions, double time,
  std::ostream* statStream)
{
  const size_t nTimes = dataNoW.size();
  
  _buffer.CopyAndWeight(dataNoW, weights, modelDataNoW);
  
  for(size_t i=0; i!=_constraints.size(); ++i)
    _constraints[i]->PrepareIteration(false, 0, false);
  
  std::vector<std::vector<DComplex> > nextSolutions(_nChannelBlocks);

  SolveResult result;
#ifndef NDEBUG
  if (solutions.size() != _nChannelBlocks) {
    std::cout << "Error: 'solutions' parameter does not have the right shape\n";
    result.iterations = 0;
    return result;
  }
#endif

  result._results.resize(_constraints.size());
  
  // Model matrix ant x [N x D] and visibility vector ant x [N x 1],
  // for each channelblock
  // The following loop allocates all structures
  std::vector<std::vector<Matrix> > gTimesCs(_nChannelBlocks);
  std::vector<std::vector<Matrix> > vs(_nChannelBlocks);
  for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
  {
    nextSolutions[chBlock].resize(_nDirections * _nAntennas);
    const size_t
      channelIndexStart = chBlock * _nChannels / _nChannelBlocks,
      channelIndexEnd = (chBlock+1) * _nChannels / _nChannelBlocks,
      curChannelBlockSize = channelIndexEnd - channelIndexStart;
    gTimesCs[chBlock].resize(_nAntennas);
    vs[chBlock].resize(_nAntennas);
    
    for(size_t ant=0; ant!=_nAntennas; ++ant)
    {
      // Model matrix [N x D] and visibility vector [N x 1]
      // Also space for the auto correlation is reserved, but they will be set to 0.
      size_t
        m = _nAntennas * nTimes * curChannelBlockSize * 4,
        n = _nDirections, nrhs = 1;
      gTimesCs[chBlock][ant] = Matrix(m, n);
      vs[chBlock][ant] = Matrix(std::max(m, n), nrhs);
    }
  }
  
  ///
  /// Start iterating
  ///
  size_t iteration = 0, constrainedIterations = 0;
  bool
    hasConverged = false,
    hasPreviouslyConverged = false,
    constraintsSatisfied = false,
    hasStalled = false;

  std::vector<double> stepMagnitudes;
  stepMagnitudes.reserve(_maxIterations);

  do {
    makeSolutionsFinite(solutions, 1);
    
    DP3::ParallelFor<size_t> loop(_nThreads);
    loop.Run(0, _nChannelBlocks, [&](size_t chBlock, size_t /*thread*/)
    {
      performScalarIteration(chBlock, gTimesCs[chBlock], vs[chBlock],
                            solutions[chBlock], nextSolutions[chBlock]);
    });
      
    makeStep(solutions, nextSolutions);
    
    constraintsSatisfied = true;

    if(statStream)
    {
      (*statStream) << iteration << '\t';
    }

    for(size_t i=0; i!=_constraints.size(); ++i)
    {
      // PrepareIteration() might change Satisfied(), and since we always want to
      // iterate at least once more when a constraint is not yet satisfied, we
      // evaluate Satisfied() before preparing.
      constraintsSatisfied = _constraints[i]->Satisfied() && constraintsSatisfied;
      _constraints[i]->PrepareIteration(hasPreviouslyConverged, iteration, iteration+1 >= _maxIterations);
      result._results[i] = _constraints[i]->Apply(nextSolutions, time, statStream);
    }
    
    if(!constraintsSatisfied)
      constrainedIterations = iteration+1;
    
    double avgSquaredDiff;
    hasConverged = assignSolutions<1>(solutions, nextSolutions, !constraintsSatisfied, avgSquaredDiff, stepMagnitudes);
    if(statStream)
    {
      (*statStream) << stepMagnitudes.back() << '\t' << avgSquaredDiff << '\n';
    }
    iteration++;
    
    hasPreviouslyConverged = hasConverged || hasPreviouslyConverged;

    if (_detectStalling && constraintsSatisfied)
      hasStalled = detectStall(iteration, stepMagnitudes);

  } while(iteration < _maxIterations && (!hasConverged || !constraintsSatisfied) && !hasStalled);
  
  // When we have not converged yet, we set the nr of iterations to the max+1, so that
  // non-converged iterations can be distinguished from converged ones.
  if((!hasConverged || !constraintsSatisfied) && !hasStalled)
    result.iterations = iteration+1;
  else
    result.iterations = iteration;
  result.constraintIterations = constrainedIterations;
  return result;
}

void MultiDirSolver::performScalarIteration(size_t channelBlockIndex,
                       std::vector<Matrix>& gTimesCs,
                       std::vector<Matrix>& vs,
                       const std::vector<DComplex>& solutions,
                       std::vector<DComplex>& nextSolutions)
{
  for(size_t ant=0; ant!=_nAntennas; ++ant)
  {
    gTimesCs[ant].zeros();
    vs[ant].zeros();
  }
  
  const size_t
    channelIndexStart = channelBlockIndex * _nChannels / _nChannelBlocks,
    channelIndexEnd = (channelBlockIndex+1) * _nChannels / _nChannelBlocks,
    curChannelBlockSize = channelIndexEnd - channelIndexStart,
    nTimes = _buffer.Data().size();
  
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
        Matrix& gTimesC1 = gTimesCs[antenna1];
        Matrix& v1 = vs[antenna1];
        Matrix& gTimesC2 = gTimesCs[antenna2];
        Matrix& v2 = vs[antenna2];
        for(size_t d=0; d!=_nDirections; ++d)
        {
          modelPtrs[d] = &_buffer.ModelData()[timeIndex][d][(channelIndexStart + baseline * _nChannels) * 4];
        }
        const Complex* dataPtr = &_buffer.Data()[timeIndex][(channelIndexStart + baseline * _nChannels) * 4];
        const size_t p1top2[4] = {0, 2, 1, 3};
        for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
        {
          const size_t
            dataIndex1 = ch-channelIndexStart + (timeIndex + antenna1 * nTimes) * curChannelBlockSize,
            dataIndex2 = ch-channelIndexStart + (timeIndex + antenna2 * nTimes) * curChannelBlockSize;
          for(size_t p1=0; p1!=4; ++p1)
          {
            size_t p2 = p1top2[p1];
            for(size_t d=0; d!=_nDirections; ++d)
            {
              std::complex<double> predicted = *modelPtrs[d];
              
              size_t solIndex1 = antenna1*_nDirections + d;
              size_t solIndex2 = antenna2*_nDirections + d;
              gTimesC2(dataIndex1*4+p1, d) = std::conj(solutions[solIndex1] * predicted); // using a* b* = (ab)*
              gTimesC1(dataIndex2*4+p2, d) = std::conj(solutions[solIndex2]) * predicted;
              
              ++modelPtrs[d]; // Goto the next polarization of this 2x2 matrix.
            }
            v1(dataIndex2*4+p2, 0) = *dataPtr;
            v2(dataIndex1*4+p1, 0) = std::conj(*dataPtr);
            ++dataPtr; // Goto the next polarization of this 2x2 matrix.
          }
        }
      }
    }
  }
  
  // The matrices have been filled; compute the linear solution
  // for each antenna.
  size_t m = _nAntennas * nTimes * curChannelBlockSize * 4;
  size_t n = _nDirections, nrhs = 1;
  QRSolver solver(m, n, nrhs);
  for(size_t ant=0; ant!=_nAntennas; ++ant) {
    // solve x^H in [g C] x^H  = v
    bool success = solver.Solve(gTimesCs[ant].data(), vs[ant].data());
    Matrix& x = vs[ant];
    if(success && x(0, 0) != 0.)
    {
      for(size_t d=0; d!=_nDirections; ++d)
        nextSolutions[ant*_nDirections + d] = x(d, 0);
    }
    else {
      for(size_t d=0; d!=_nDirections; ++d)
        nextSolutions[ant*_nDirections + d] = std::numeric_limits<double>::quiet_NaN();
    }
  }
}

MultiDirSolver::SolveResult MultiDirSolver::processFullMatrix(
  const std::vector<Complex *>& dataNoW,
  const std::vector<float*>& weights,
  const std::vector<std::vector<Complex *> >& modelDataNoW,
  std::vector<std::vector<DComplex> >& solutions, double time,
  std::ostream* statStream)
{
  // This algorithm is basically the same as the scalar algorithm,
  // but visibility values are extended to 2x2 matrices and concatenated
  // in the matrix equations as block matrices. One difference is that
  // order of operations are important because of the non-commutativity of
  // matrix multiplication, as well as that A^H B^H = [BA]^H.
  //  
  // The approach:
  // First we pre-apply the left-hand solutions to the model to make JM. Each
  // 2x2 coherence matrix Ji is matrix-multied by the lh solutions, for all
  // directions, and visibilities (times x channels).
  //   JMi = Ji Mi
  // These are stacked in matrix JM :
  //        JM0_d0 JM1_d0 ...
  //   JM = JM0_d1 JM1_d1 
  //        ...          
  // such that JM is a (2D) rows x (2N) col matrix, N=nvis, D=ndir.
  // The solved 2D x 2 solution matrix is similarly formed with the solution
  // values:
  //       ( J0 )
  //   J = ( J1 )
  //       ( .. )
  // And the 2N x 2 visibility matrix as well:
  //       ( V0 )
  //   V = ( V1 )
  //       ( .. )
  // And we solve the equation:
  //   'JM' J^H = V
  // With dimensions:
  //   [ 2N x 2D ] [ 2D x 2 ] = [ 2N x 2 ]

  const size_t nTimes = dataNoW.size();
  
  _buffer.CopyAndWeight(dataNoW, weights, modelDataNoW);
  
  for(size_t i=0; i!=_constraints.size(); ++i)
    _constraints[i]->PrepareIteration(false, 0, false);
  
  std::vector<std::vector<DComplex> > nextSolutions(_nChannelBlocks);

  SolveResult result;
#ifndef NDEBUG
  if (solutions.size() != _nChannelBlocks) {
    std::cout << "Error: 'solutions' parameter does not have the right shape\n";
    result.iterations = 0;
    return result;
  }
#endif

  result._results.resize(_constraints.size());
  
  // Dimensions for each channelblock:
  // Model matrix ant x [2N x 2D] and visibility matrix ant x [2N x 2],
  // The following loop allocates all structures
  std::vector<std::vector<Matrix> > gTimesCs(_nChannelBlocks);
  std::vector<std::vector<Matrix> > vs(_nChannelBlocks);
  for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
  {
    nextSolutions[chBlock].resize(_nDirections * _nAntennas * 4);
    const size_t
      channelIndexStart = chBlock * _nChannels / _nChannelBlocks,
      channelIndexEnd = (chBlock+1) * _nChannels / _nChannelBlocks,
      curChannelBlockSize = channelIndexEnd - channelIndexStart;
    gTimesCs[chBlock].resize(_nAntennas);
    vs[chBlock].resize(_nAntennas);
    
    for(size_t ant=0; ant!=_nAntennas; ++ant)
    {
      // Model matrix [2N x 2D] and visibility matrix [2N x 2]
      // Space for the auto correlation is also reserved, but they will be set to 0.
      size_t m = _nAntennas * nTimes * curChannelBlockSize * 2;
      size_t n = _nDirections * 2, nrhs = 2;
      gTimesCs[chBlock][ant] = Matrix(m, n);
      vs[chBlock][ant] = Matrix(std::max(n, m), nrhs);
    }
  }
  
  ///
  /// Start iterating
  ///
  size_t iteration = 0, constrainedIterations = 0;
  bool hasConverged = false,
       hasPreviouslyConverged = false,
       constraintsSatisfied = false,
       hasStalled = false;

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(_maxIterations);

  do {
    makeSolutionsFinite(solutions, 4);
    
    DP3::ParallelFor<size_t> loop(_nThreads);
    loop.Run(0, _nChannelBlocks, [&](size_t chBlock, size_t /*thread*/)
    {
      performFullMatrixIteration(chBlock, gTimesCs[chBlock], vs[chBlock],
                                solutions[chBlock], nextSolutions[chBlock]);
    });
    
    makeStep(solutions, nextSolutions);

    if(statStream)
    {
      (*statStream) << iteration << '\t';
    }
    
    constraintsSatisfied = true;
    for(size_t i=0; i!=_constraints.size(); ++i)
    {
      constraintsSatisfied = _constraints[i]->Satisfied() && constraintsSatisfied;
      _constraints[i]->PrepareIteration(hasPreviouslyConverged, iteration, iteration+1 >= _maxIterations);
      result._results[i] = _constraints[i]->Apply(nextSolutions, time, statStream);
    }
    
    if(!constraintsSatisfied)
      constrainedIterations = iteration+1;
    
    double avgSquaredDiff;
    hasConverged = assignSolutions<4>(solutions, nextSolutions, !constraintsSatisfied, avgSquaredDiff, step_magnitudes);
    if(statStream)
    {
      (*statStream) << step_magnitudes.back() << '\t' << avgSquaredDiff << '\n';
    }
    iteration++;
    
    hasPreviouslyConverged = hasConverged || hasPreviouslyConverged;

    if (_detectStalling && constraintsSatisfied)
      hasStalled = detectStall(iteration, step_magnitudes);

  } while(iteration < _maxIterations && (!hasConverged || !constraintsSatisfied) && !hasStalled);
 
  // When we have not converged yet, we set the nr of iterations to the max+1, so that
  // non-converged solves can be distinguished from converged ones.
  if((!hasConverged || !constraintsSatisfied) && !hasStalled)
    result.iterations = iteration+1;
  else
    result.iterations = iteration;
  result.constraintIterations = constrainedIterations;
  return result;
}

bool MultiDirSolver::detectStall(size_t iteration, const std::vector<double>& step_magnitudes) const
{
  if (iteration<30) {
    return false;
  } else {
    return std::abs(step_magnitudes[iteration-1]/step_magnitudes[iteration-2]-1) < 1.e-4 &&
           std::abs(step_magnitudes[iteration-2]/step_magnitudes[iteration-3]-1) < 1.e-4;
  }
}

void MultiDirSolver::performFullMatrixIteration(size_t channelBlockIndex,
                             std::vector<Matrix>& gTimesCs,
                             std::vector<Matrix>& vs,
                             const std::vector<DComplex>& solutions,
                             std::vector<DComplex>& nextSolutions)
{
  for(size_t ant=0; ant!=_nAntennas; ++ant)
  {
    gTimesCs[ant].zeros();
    vs[ant].zeros();
  }
  
  const size_t
    channelIndexStart = channelBlockIndex * _nChannels / _nChannelBlocks,
    channelIndexEnd = (channelBlockIndex+1) * _nChannels / _nChannelBlocks,
    curChannelBlockSize = channelIndexEnd - channelIndexStart,
    nTimes = _buffer.Data().size();
  
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
        // This equation is solved:
        //   J_1 M J_2^H = V
        // for visibilities of the 'normal' correlation ant1 x ant2^H.
        // Since in this equation antenna2 is solved, the solve matrices are
        // called gTimesC2 and v2. The index into these matrices is depending
        // on antenna1, hence the index is called dataIndex1.
        //
        // To use visibilities of correlation ant2 x ant1^H to solve ant1, an
        // extra Herm transpose on M and V is required. The equation is:
        //   J_2 M^H J_1^H = V^H,
        // and the relevant matrices/index are called gTimesC1, v1 and dataIndex2.
        Matrix
          &gTimesC1 = gTimesCs[antenna1],
          &v1 = vs[antenna1],
          &gTimesC2 = gTimesCs[antenna2],
          &v2 = vs[antenna2];
        for(size_t d=0; d!=_nDirections; ++d)
          modelPtrs[d] = &_buffer.ModelData()[timeIndex][d][(channelIndexStart + baseline * _nChannels) * 4];
        const Complex* dataPtr = &_buffer.Data()[timeIndex][(channelIndexStart + baseline * _nChannels) * 4];
        for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
        {
          const size_t
            dataIndex1 = 2 * (ch-channelIndexStart + (timeIndex + antenna1 * nTimes) * curChannelBlockSize),
            dataIndex2 = 2 * (ch-channelIndexStart + (timeIndex + antenna2 * nTimes) * curChannelBlockSize);
            
          for(size_t d=0; d!=_nDirections; ++d)
          {
            MC2x2
              modelMat(modelPtrs[d]),
              gTimesC1Mat, gTimesC2Mat;
            size_t solIndex1 = (antenna1*_nDirections + d) * 4;
            size_t solIndex2 = (antenna2*_nDirections + d) * 4;
            Matrix2x2::ATimesB(gTimesC2Mat.Data(), &solutions[solIndex1], modelMat.Data());
            Matrix2x2::ATimesHermB(gTimesC1Mat.Data(), &solutions[solIndex2], modelMat.Data());
            for(size_t p=0; p!=4; ++p)
            {
              gTimesC2(dataIndex1+(p/2), d*2+p%2) = gTimesC2Mat[p];
              gTimesC1(dataIndex2+(p/2), d*2+p%2) = gTimesC1Mat[p];
            }
            
            modelPtrs[d] += 4; // Goto the next 2x2 matrix.
          }
          for(size_t p=0; p!=4; ++p)
          {
            v1(dataIndex2+(p%2), p/2) = std::conj(*dataPtr);
            v2(dataIndex1+(p/2), p%2) = *dataPtr; // note that this also performs the Herm transpose
            ++dataPtr;  // Goto the next element of the 2x2 matrix.
          }
        }
      }
    }
  }
  
  // The matrices have been filled; compute the linear solution
  // for each antenna.

  size_t m = _nAntennas * nTimes * curChannelBlockSize * 2;
  size_t n = _nDirections * 2, nrhs = 2;
  QRSolver solver(m, n, nrhs);

  for(size_t ant=0; ant!=_nAntennas; ++ant) {
    // solve x^H in [g C] x^H  = v
    bool success = solver.Solve(gTimesCs[ant].data(), vs[ant].data());
    Matrix& x = vs[ant];
    if(success && x(0, 0) != 0.)
    {
      for(size_t d=0; d!=_nDirections; ++d)
      {
        for(size_t p=0; p!=4; ++p) {
          // The conj transpose is also performed at this point (note swap of % and /)
          nextSolutions[(ant*_nDirections + d)*4 + p] = std::conj(x(d*2+p%2, p/2));
        }
      }
    }
    else {
      for(size_t i=0; i!=_nDirections*4; ++i) {
        nextSolutions[ant*_nDirections*4 + i] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

void MultiDirSolver::showTimings (std::ostream& os, double duration) const {
  //os << "                " << std::fixed << std::setprecision(2) << _timerSolve.Seconds()/duration << "% spent in solve" << std::endl;
  //os << "                " << std::fixed << std::setprecision(2) << _timerFillMatrices.Seconds()/duration << "% spent in filling matrices" << std::endl;
  if (!_constraints.empty()) {
    //os << "                " << std::fixed << std::setprecision(2) << _timerConstrain.Seconds()/duration << "% spent in constraints" << std::endl;
    for (auto& constraint: _constraints) {
      constraint->showTimings(os, duration);
    }
  }
}
