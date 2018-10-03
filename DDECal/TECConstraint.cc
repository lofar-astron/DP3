#include "TECConstraint.h"

#include "../Common/ParallelFor.h"

TECConstraintBase::TECConstraintBase(Mode mode, size_t nThreads) :
  _mode(mode),
  _nThreads(nThreads),
  _phaseFitters()
{
}

void TECConstraintBase::initialize(const double* frequencies) {
  _phaseFitters.resize(_nThreads);

  for(size_t i=0; i!=_phaseFitters.size(); ++i)
  {
    _phaseFitters[i].SetChannelCount(_nChannelBlocks);
    std::memcpy(_phaseFitters[i].FrequencyData(), frequencies,
                sizeof(double) * _nChannelBlocks);
  }
  _weights.assign(_nChannelBlocks*_nAntennas, 1.0);
  initializeChild();
}

void TECConstraintBase::SetWeights(const std::vector<double>& weights) {
  _weights = weights;
}

void ApproximateTECConstraint::initializeChild()
{
  _pwFitters.resize(_nThreads);
  _threadData.resize(_pwFitters.size());
  _threadFittedData.resize(_pwFitters.size());
  _threadWeights.resize(_pwFitters.size());
  for(size_t threadId=0; threadId!=_pwFitters.size(); ++threadId)
  {
    _threadData[threadId].resize(_nChannelBlocks);
    _threadFittedData[threadId].resize(_nChannelBlocks);
    _threadWeights[threadId].resize(_nChannelBlocks);
  }
  
  if(_fittingChunkSize == 0)
  {
    size_t
      n = _phaseFitters.front().Size();
    const double
      startFreq = _phaseFitters.front().FrequencyData()[0],
      endFreq = _phaseFitters.front().FrequencyData()[n-1];
      _fittingChunkSize = PieceWisePhaseFitter::CalculateChunkSize(startFreq, endFreq, n);
  }
  for(size_t i=0; i!=_pwFitters.size(); ++i)
    _pwFitters[i].SetChunkSize(_fittingChunkSize);
}

void TECConstraintBase::applyReferenceAntenna(std::vector<std::vector<dcomplex> >& solutions) const
{
  // TODO chose this more cleverly?
  size_t refAntenna = 0;

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
}

std::vector<Constraint::Result> TECConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double,
    std::ostream* /*statStream*/)
{
  size_t nRes = 3;
  if(_mode == TECOnlyMode) {
    nRes = 2; // TEC and error
  }
  else {
    nRes = 3; // TEC, phase and error
  }

  std::vector<Constraint::Result> res(nRes);
  res[0].vals.resize(_nAntennas*_nDirections);
  res[0].weights.resize(_nAntennas*_nDirections);
  res[0].axes="ant,dir,freq";
  res[0].name="tec";
  res[0].dims.resize(3);
  res[0].dims[0]=_nAntennas;
  res[0].dims[1]=_nDirections;
  res[0].dims[2]=1;
  if(_mode == TECAndCommonScalarMode) {
    res[1]=res[0];
    res[1].name="phase";
  }
  res.back()=res[0];
  res.back().name="error";

  // Divide out the reference antenna
  applyReferenceAntenna(solutions);
  
  DP3::ParallelFor<size_t> loop(_nThreads);
  loop.Run(0, _nAntennas*_nDirections, [&](size_t solutionIndex, size_t thread)
  {
    size_t antennaIndex = solutionIndex/_nDirections;

    // Flag channels where calibration yielded inf or nan
    double weightSum = 0.0;
    for(size_t ch=0; ch!=_nChannelBlocks; ++ch) {
      if(std::isfinite(solutions[ch][solutionIndex].real()) &&
        std::isfinite(solutions[ch][solutionIndex].imag()))
      {
        _phaseFitters[thread].PhaseData()[ch] = std::arg(solutions[ch][solutionIndex]);
        _phaseFitters[thread].WeightData()[ch] = _weights[antennaIndex*_nChannelBlocks + ch];
        weightSum += _weights[antennaIndex*_nChannelBlocks + ch];
      }
      else {
        _phaseFitters[thread].PhaseData()[ch] = 0.0;
        _phaseFitters[thread].WeightData()[ch] = 0.0;
      }
    }
    
    double alpha, beta=0.0;
    if(_mode == TECOnlyMode) {
      res.back().vals[solutionIndex]=_phaseFitters[thread].FitDataToTEC1Model(alpha);
    } else {
      res.back().vals[solutionIndex]=_phaseFitters[thread].FitDataToTEC2Model(alpha, beta);
    }
    res.back().weights[solutionIndex] = weightSum;

    res[0].vals[solutionIndex] = alpha / -8.44797245e9;
    res[0].weights[solutionIndex] = weightSum;
    if(_mode == TECAndCommonScalarMode) {
      res[1].vals[solutionIndex] = beta;
      res[1].weights[solutionIndex] = weightSum;
    }
    
    for(size_t ch=0; ch!=_nChannelBlocks; ++ch) 
    {
      solutions[ch][solutionIndex] = std::polar<double>(1.0, _phaseFitters[thread].PhaseData()[ch]);
    }
  });

  return res;
}

std::vector<Constraint::Result> ApproximateTECConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double time,
    std::ostream* statStream)
{
  if(_finishedApproximateStage)
    return TECConstraint::Apply(solutions, time, statStream);
  else {
    applyReferenceAntenna(solutions);
    
    DP3::ParallelFor<size_t> loop(_nThreads);
    loop.Run(0, _nAntennas*_nDirections, [&](size_t solutionIndex, size_t thread)
    {
      size_t antennaIndex = solutionIndex/_nDirections;
      std::vector<double>& data = _threadData[thread];
      std::vector<double>& fittedData = _threadFittedData[thread];
      std::vector<double>& weights = _threadWeights[thread];
      
      // Flag channels where calibration yielded inf or nan
      for(size_t ch=0; ch!=_nChannelBlocks; ++ch) {
        if(std::isfinite(solutions[ch][solutionIndex].real()) &&
          std::isfinite(solutions[ch][solutionIndex].imag()))
        {
          data[ch] = std::arg(solutions[ch][solutionIndex]);
          weights[ch] = _weights[antennaIndex*_nChannelBlocks + ch];
        }
        else {
          data[ch] = 0.0;
          weights[ch] = 0.0;
        }
      }
      
      // TODO might be nice to make it a user option whether to break or not
      _pwFitters[thread].SlidingFitWithBreak(_phaseFitters[thread].FrequencyData(), data.data(), weights.data(), fittedData.data(), data.size());

      for(size_t ch=0; ch!=_nChannelBlocks; ++ch) 
      {
        solutions[ch][solutionIndex] = std::polar<double>(1.0, fittedData[ch]);
      }
    });

    return std::vector<Constraint::Result>();
  }
}
