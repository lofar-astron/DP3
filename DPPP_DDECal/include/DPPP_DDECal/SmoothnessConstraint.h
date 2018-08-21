#include "Constraint.h"
#include "KernelSmoother.h"

#ifndef SMOOTHNESS_CONSTRAINT_H
#define SMOOTHNESS_CONSTRAINT_H

class SmoothnessConstraint : public Constraint
{
public:
  typedef std::complex<double> dcomplex;
  typedef KernelSmoother<dcomplex, double> Smoother;
  
  SmoothnessConstraint(double bandwidthHz);
  
  std::vector<Constraint::Result> Apply(
    std::vector<std::vector<dcomplex> >& solutions, double, std::ostream* statStream) final override;
  
  void SetWeights(const std::vector<double> &weights) final override {
    _weights = weights;
  }
  
  void Initialize(const double* frequencies);
  
  virtual void InitializeDimensions(size_t nAntennas,
                                    size_t nDirections,
                                    size_t nChannelBlocks) final override;
                                    
  struct FitData
  {
    FitData(const double* frequencies, size_t n, Smoother::KernelType kernelType, double kernelBandwidth)
      : smoother(frequencies, n, kernelType, kernelBandwidth),
      data(n), weight(n)
    { }
    
    Smoother smoother;
    std::vector<dcomplex> data;
    std::vector<double> weight;
  };
  std::vector<FitData> _fitData;
  std::vector<double> _frequencies, _weights;
  Smoother::KernelType _kernelType;
  double _bandwidth;
};

#endif

