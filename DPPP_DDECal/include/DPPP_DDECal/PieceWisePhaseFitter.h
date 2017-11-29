#ifndef PIECE_WISE_PHASE_FITTER_H
#define PIECE_WISE_PHASE_FITTER_H

#include <algorithm>
#include <iterator>
#include <cmath>

class PieceWisePhaseFitter
{
public:
  PieceWisePhaseFitter() : _chunkSize(0) { }
  
  /**
   * Constructor.
   * @param chunkSize Size of chunk in number of samples.
   */
  PieceWisePhaseFitter(size_t chunkSize) : _chunkSize(chunkSize) { }
  
  size_t ChunkSize() const { return _chunkSize; }
  void SetChunkSize(size_t chunkSize) { _chunkSize = chunkSize; }
  
  /**
  * Unwrap a range of phase values.
  */
  template<typename Iter>
  static void Unwrap(Iter first, Iter last)
  {
    for(Iter iter = first+1; iter!=last; ++iter)
    {
      while(*iter - (*(iter-1)) > M_PI)
        *iter -= 2.0*M_PI;
      while((*(iter-1) - *iter) > M_PI)
        *iter += 2.0*M_PI;
    }
  }
  
  static size_t CalculateChunkSize(double startFrequencyHz, double endFrequencyHz, size_t channelCount)
  {
    // it seems that 10 chunks per octave seems reasonable
    double nOctaves = (log(endFrequencyHz) - log(startFrequencyHz)) / M_LN2;
    if(nOctaves > 0.0)
      return std::min<size_t>(channelCount, ceil(channelCount / (10.0 * nOctaves)));
    else
      return channelCount;
  }
  
  /**
  * Perform a piece-wise fit of linear phase gradients.
  * 
  * The data is divided into chunks with given size, and in each chunk a line is fitted (offset & gradient).
  * The fitted line of each chunk is written into fittedData. The resulting line will have a 'jump' between each
  * chunk.
  * @param nu Frequencies of the channels (in Hz).
  * @param data The values, same size as @c nu.
  * @param fittedData A vector of at least size @c nu.size(), in which the fitted lines are stored.
  */
  void PieceWiseFit(const std::vector<double>& nu, const std::vector<double>& data, std::vector<double>& fittedData)
  {
    size_t chunkSize = std::min(data.size(), _chunkSize);
    if(chunkSize==1)
      fittedData = data;
    else
    {
      for(size_t ch=0; ch<data.size(); ch+=chunkSize)
      {
        size_t pos = ch;
        if(pos > data.size() - chunkSize)
          pos = data.size() - chunkSize;
        fitChunk(chunkSize, pos, nu, data, fittedData);
      }
    }
  }
  
  void PieceWiseFit(const std::vector<double>& nu, const std::vector<double>& data, const double* weights, std::vector<double>& fittedData)
  {
    size_t chunkSize = std::min(data.size(), _chunkSize);
    if(chunkSize==1)
      fittedData = data;
    else
    {
      for(size_t ch=0; ch<data.size(); ch+=chunkSize)
      {
        size_t pos = ch;
        if(pos > data.size() - chunkSize)
          pos = data.size() - chunkSize;
        fitChunk(chunkSize, pos, nu, data, weights, fittedData);
      }
    }
  }
  
  /**
  * Performs a sliding fit of linear phase gradients (without weighting).
  * See other SlidingFit() overload.
  */
  void SlidingFit(const double* nu, const double* data, double* fittedData, size_t n)
  {
    const size_t
      chunkSize = std::min(_chunkSize, n);
    if(chunkSize==1)
      std::copy(data, data+n, fittedData);
    else
    {
      const size_t
        leftEdge = chunkSize/2,
        rightEdge = n - chunkSize + chunkSize/2;
      
      double a, b;
      PieceWisePhaseFitter::fitSlope(&data[0], &nu[0], chunkSize, a, b);
      for(size_t ch=0; ch!=leftEdge; ++ch)
        fittedData[ch] = a + b * nu[ch];
      
      for(size_t ch=0; ch!=n - chunkSize; ++ch)
      {
        PieceWisePhaseFitter::fitSlope(&data[ch], &nu[ch], chunkSize, a, b);
        fittedData[ch + chunkSize/2] = a + b * nu[ch + chunkSize/2];
      }
      
      PieceWisePhaseFitter::fitSlope(&data[n - chunkSize], &nu[n - chunkSize], chunkSize, a, b);
      for(size_t ch=rightEdge; ch!=n; ++ch)
        fittedData[ch] = a + b * nu[ch];
    }
  }
  
  /**
  * Performs a sliding fit of linear phase gradients (with weighting).
  * 
  * A window of size @c chunkSize is slid over the data, and a line is fitted at each position. The values near the left
  * and right side of the border are evaluated with the fit to the most left/right window. This is slower than
  * @ref PieceWiseFit(), but it produces a smoother & more accurate fit.
  * The fit is not necessarily entirely smooth, because the line fit tries to guess the right wrapping of the samples.
  * Since that is not a linear process, this can still cause slight jumps.
  * 
  * @param nu Frequencies of the channels (in Hz).
  * @param data The values, same size as @c nu.
  * @param weights Inverse-variance weights.
  * @param fittedData A vector of at least size @c nu.size(), in which the fitted lines are stored.
  */
  void SlidingFit(const double* nu, const double* data, const double* weights, double* fittedData, size_t n)
  {
    const size_t
      chunkSize = std::min(_chunkSize, n);
    if(chunkSize==1)
      std::copy(data, data+n, fittedData);
    else
    {
    const size_t
      leftEdge = chunkSize/2,
      rightEdge = n - chunkSize + chunkSize/2;
    
      double a, b;
      PieceWisePhaseFitter::fitSlope(&data[0], &nu[0], &weights[0], chunkSize, a, b);
      for(size_t ch=0; ch!=leftEdge; ++ch)
        fittedData[ch] = a + b * nu[ch];
      
      for(size_t ch=0; ch!=n - chunkSize; ++ch)
      {
        PieceWisePhaseFitter::fitSlope(&data[ch], &nu[ch], &weights[ch], chunkSize, a, b);
        fittedData[ch + chunkSize/2] = a + b * nu[ch + chunkSize/2];
      }
      
      PieceWisePhaseFitter::fitSlope(&data[n - chunkSize], &nu[n - chunkSize], &weights[n - chunkSize], chunkSize, a, b);
      for(size_t ch=rightEdge; ch!=n; ++ch)
        fittedData[ch] = a + b * nu[ch];
    }
  }
  
  size_t SlidingFitWithBreak(const double* nu, const double* data, const double* weights, double* fittedData, size_t n)
  {
    SlidingFit(nu, data, weights, fittedData, n);
    size_t bp = BreakPoint(data, weights, fittedData, n);
    SlidingFit(nu, data, weights, fittedData, bp);
    SlidingFit(nu+bp, data+bp, weights+bp, fittedData+bp, n-bp);
    return bp;
  }
  
  /**
   * Calculates a weighted median. The weighted median is defined as the
   * sample x for which the sum of all values smaller than x is less than half of the
   * total weight, and idem for all values larger than x.
   * For example, given the sequence
   * 1,2,3,4,5 with weights 1,100,1,1,1; the weighted median is '2'.
   * (sum of weight=104, values < 2 a have weight of 1, values > 2 have a weight of 3).
   * This function is internally used by the piece-wise phase fitter.
   * @param values A vector of (value, weight) pairs.
   * @returns The weighted median.
   */
  static double WeightedMedian(std::vector<std::pair<double, double>>& values)
  {
    std::sort(values.begin(), values.end());
    // calculate total weight
    double sum = 0.0;
    for(std::pair<double,double>& v : values)
      sum += v.second;

    int index = 0;
    // prefixSum is the weight sum of everything after `index`
    double prefixSum = sum - values[0].second;
    while(prefixSum > sum/2)
    {
      ++index;
      prefixSum -= values[index].second;
    }
    return values[index].first;
  }
  
  size_t BreakPoint(const double* data, const double* weights, const double* fittedData, size_t n)
  {
    double largest = 0.0;
    size_t index = 0;
    for(size_t i=0; i!=n-_chunkSize; ++i)
    {
      double curSum = 0.0;
      for(size_t j=0; j!=_chunkSize; ++j)
      {
        double dist = data[i+j] - fittedData[i+j];
        curSum += fabs(dist) * weights[i+j];
      }
      if(curSum > largest)
      {
        largest = curSum;
        index = i + _chunkSize/2;
      }
    }
    return index;
  }
  
private:
  size_t _chunkSize;
  /**
   * This is used to avoid having to reallocate data every fit
   */
  std::vector<double> _tempArray;
  std::vector<std::pair<double,double>> _tempWeightedArray;
  
  /**
  * Fit a + bx
  */
  void fitSlope(const double* data, const double* frequencies, size_t n, double& a, double& b)
  {
    double
      sum_x = 0.0,
      sum_xSq = 0.0,
      mean_y = 0.0,
      ss_xy = 0.0;
    std::vector<double>& d = _tempArray;
    d.assign(data, data+n);
    size_t qCounts[4] = {0,0,0,0};
    for(double& v : d)
    {
      if(v > M_PI*1.5)
        ++qCounts[3];
      else if(v > M_PI*1.0)
        ++qCounts[2];
      else if(v > M_PI*0.5)
        ++qCounts[1];
      else
        ++qCounts[0];
    }
    size_t maxQ = 0;
    for(size_t i=1; i!=4; ++i)
    {
      if(qCounts[i] > qCounts[maxQ])
        maxQ = i;
    }
    double wrapTo = (maxQ*0.5+0.25)*M_PI;
    for(double& v : d)
    {
      if(v - wrapTo > M_PI)
        v -= 2.0*M_PI;
      else if(v - wrapTo < -M_PI)
        v += 2.0*M_PI;
    }
    std::nth_element(d.data(), d.data()+n/2, d.data()+n);
    wrapTo = d[n/2];
    for(size_t i=0; i!=n; ++i)
    {
      double x = data[i];
      if(x - wrapTo > M_PI)
        x -= 2.0*M_PI;
      else if(x - wrapTo < -M_PI)
        x += 2.0*M_PI;
      sum_x += frequencies[i];
      sum_xSq += frequencies[i] * frequencies[i];
      mean_y += x;
      ss_xy += frequencies[i] * x;
    }
    mean_y /= n;
    double mean_x = sum_x / n;
    double ss_xx = sum_xSq - mean_x * sum_x;
    ss_xy -= sum_x * mean_y;
    b = ss_xy / ss_xx;
    a = mean_y - b * mean_x;
  }
  
  template<typename Iter1, typename Iter2>
  double weightedMedian(Iter1 valBegin, Iter1 valEnd, Iter2 weightBegin)
  {
    _tempWeightedArray.resize(std::distance(valBegin, valEnd));
    auto wIter = _tempWeightedArray.begin();
    while(valBegin != valEnd)
    {
      *wIter = std::make_pair(*valBegin, *weightBegin);
      ++wIter;
      ++valBegin;
      ++weightBegin;
    }
    return WeightedMedian(_tempWeightedArray);
  }
  
  /**
  * Fit a + bx with weights
  */
  void fitSlope(const double* data, const double* frequencies, const double* weights, size_t n, double& a, double& b)
  {
    size_t qCounts[4] = {0,0,0,0};
    for(size_t i=0; i!=n; ++i)
    {
      if(data[i] > M_PI*1.5)
        qCounts[3] += weights[i];
      else if(data[i] > M_PI*1.0)
        qCounts[2] += weights[i];
      else if(data[i] > M_PI*0.5)
        qCounts[1] += weights[i];
      else
        qCounts[0] += weights[i];
    }
    size_t maxQ = 0;
    for(size_t i=1; i!=4; ++i)
    {
      if(qCounts[i] > qCounts[maxQ])
        maxQ = i;
    }
    double wrapTo = (maxQ*0.5+0.25)*M_PI;
    std::vector<double>& d = _tempArray;
    d.resize(n);
    for(size_t i=n/4; i!=n*3/4; ++i)
    {
      double v = data[i];
      if(v - wrapTo > M_PI)
        v -= 2.0*M_PI;
      else if(v - wrapTo < -M_PI)
        v += 2.0*M_PI;
      d[i] = v;
    }
    wrapTo = weightedMedian(d.data()+n/4, d.data()+n*3/4, weights+n/4);
    
    double
      sum_x = 0.0,
      sum_xSq = 0.0,
      mean_y = 0.0,
      ss_xy = 0.0,
      sum_weight = 0.0;
      
    for(size_t i=0; i!=n; ++i)
    {
      double x = data[i];
      if(x - wrapTo > M_PI)
        x -= 2.0*M_PI;
      else if(x - wrapTo < -M_PI)
        x += 2.0*M_PI;
      sum_x += frequencies[i] * weights[i];
      sum_xSq += frequencies[i] * frequencies[i] * weights[i];
      mean_y += x * weights[i];
      ss_xy += frequencies[i] * x * weights[i];
      sum_weight += weights[i];
    }
    if(sum_weight == 0.0)
    {
      a = 0.0;
      b = 0.0;
    }
    else {
      mean_y /= sum_weight;
      double mean_x = sum_x / sum_weight;
      double ss_xx = sum_xSq - mean_x * sum_x;
      ss_xy -= sum_x * mean_y;
      b = ss_xy / ss_xx;
      a = mean_y - b * mean_x;
    }
  }
    
  void fitChunk(size_t chunkSize, size_t pos, const std::vector<double>& nu, const std::vector<double>& data, std::vector<double>& fittedData)
  {
    double a, b;
    PieceWisePhaseFitter::fitSlope(&data[pos], &nu[pos], chunkSize, a, b);
    for(size_t ch=0; ch!=chunkSize; ++ch)
      fittedData[ch+pos] = a + b * nu[ch+pos];
  }
  
  void fitChunk(size_t chunkSize, size_t pos, const std::vector<double>& nu, const std::vector<double>& data, const double* weights, std::vector<double>& fittedData)
  {
    double a, b;
    PieceWisePhaseFitter::fitSlope(&data[pos], &nu[pos], &weights[pos], chunkSize, a, b);
    for(size_t ch=0; ch!=chunkSize; ++ch)
      fittedData[ch+pos] = a + b * nu[ch+pos];
  }
  
  /**
  * No longer used -- I tried removing samples that deviated far from the fitted line and refit, but this seemed not to help.
  */
  static bool removeDeviations(std::vector<double>& data, std::vector<double>& frequencies, double a, double b, double phaseThreshold)
  {
    if(data.empty())
      return false;
    double maxDev = 0.0;
    size_t index = 0;
    for(size_t i=0; i!=data.size(); ++i)
    {
      double deviation = std::fabs(data[i] - (a + b * frequencies[i]));
      if(deviation >= maxDev)
      {
        maxDev = deviation;
        index = i;
      }
    }
    if(maxDev > phaseThreshold)
    {
      data.erase(data.begin() + index);
      frequencies.erase(frequencies.begin() + index);
      return true;
    }
    else return false;
  }
};

#endif
