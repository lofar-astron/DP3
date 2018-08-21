#ifndef KERNEL_SMOOTHER_H
#define KERNEL_SMOOTHER_H

#include <cmath>
#include <stdexcept>
#include <vector>

template<typename DataType, typename NumType>
class KernelSmoother
{
public:
  enum KernelType {
    RectangularKernel,
    TriangularKernel,
    /** Gaussian, trimmed off at 3 sigma */
    GaussianKernel,
    /** The Epanechnikov kernel is a quadratic kernel, given by 3/4 (1 - x^2) */
    EpanechnikovKernel
  };
  
  KernelSmoother(const NumType* frequencies, size_t n, KernelType kernelType, NumType kernelBandwidth) :
    _frequencies(frequencies, frequencies+n),
    _scratch(n),
    _kernelType(kernelType),
    _bandwidth(kernelBandwidth)
  {
  }
  
  NumType Kernel(NumType distance) const
  {
    NumType x = distance / _bandwidth;
    if(x < NumType(-1.0) || x > NumType(1.0))
      return NumType(0.0);
    else {
      switch(_kernelType)
      {
      case RectangularKernel:
      default:
        return NumType(0.5);
      case TriangularKernel:
        return x >= NumType(0.0) ? (NumType(1.0) - x) : (NumType(1.0) + x);
      case GaussianKernel:
        // e^(-x^2 / sigma^2), sigma = bandwidth / 3.
        return std::exp(-x*x*NumType(9.0));
      case EpanechnikovKernel:
        // 3/4 * (1-x)^2;
        x = NumType(1.0) - x;
        return (NumType(3.0) / NumType(4.0)) * x * x;
      }
    }
  }
  
  void Smooth(DataType* data, const NumType* weight)
  {
    size_t n = _frequencies.size();
     
    size_t
      bandLeft = 0,
      // find right kernel value for first element
      bandRight = std::lower_bound(_frequencies.begin(), _frequencies.end(), _frequencies[0] + _bandwidth * 0.5) - _frequencies.begin() + 1;
    
    for(size_t i=0; i!=n; ++i)
    {
      // If a boundary is further than half the bandwidth away, move boundary
      while(_frequencies[bandLeft] < _frequencies[i] - _bandwidth * 0.5)
        ++bandLeft;
      while(bandRight!=n && _frequencies[bandRight] < _frequencies[i] + _bandwidth * 0.5)
        ++bandRight;
      
      // A value of 1 is added to make sure we are not skipping a value because of rounding errors
      // (kernel will be zero past boundaries, so including an unnecessary value has no effect)
      size_t start = bandLeft > 0 ? bandLeft-1 : 0;
      size_t end = bandRight < n ? bandRight+1 : n;
      
      DataType sum(0.0);
      NumType weightSum(0.0);
      //std::cout << start << " -> " << end << " (" << _frequencies[start] << " -> " << _frequencies[end] << ")\n";
      for(size_t j=start; j!=end; ++j)
      {
        double distance = _frequencies[i] - _frequencies[j];
        double w = Kernel(distance) * weight[j];
        sum += data[j] * w;
        weightSum += w;
      }
      if(weightSum == 0.0)
        _scratch[i] = 0.0;
      else
        _scratch[i] = sum / weightSum;
    }
    std::copy(_scratch.begin(), _scratch.end(), data);
  }
  
private:
  std::vector<NumType> _frequencies;
  std::vector<DataType> _scratch;
  enum KernelType _kernelType;
  NumType _bandwidth;
};

#endif

