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
    for(size_t ch=0; ch<data.size(); ch+=chunkSize)
    {
      size_t pos = ch;
      if(pos > data.size() - chunkSize)
        pos = data.size() - chunkSize;
      fitChunk(chunkSize, pos, nu, data, fittedData);
    }
  }
  
  void PieceWiseFit(const std::vector<double>& nu, const std::vector<double>& data, const double* weights, std::vector<double>& fittedData)
  {
    size_t chunkSize = std::min(data.size(), _chunkSize);
    for(size_t ch=0; ch<data.size(); ch+=chunkSize)
    {
      size_t pos = ch;
      if(pos > data.size() - chunkSize)
        pos = data.size() - chunkSize;
      fitChunk(chunkSize, pos, nu, data, weights, fittedData);
    }
  }
  
  /**
  * Performs a sliding fit of linear phase gradients (without weighting).
  * See other SlidingFit() overload.
  */
  void SlidingFit(const double* nu, const std::vector<double>& data, std::vector<double>& fittedData)
  {
    const size_t
      chunkSize = std::min(_chunkSize, data.size()),
      leftEdge = chunkSize/2,
      rightEdge = data.size() - chunkSize + chunkSize/2;
    
    double a, b;
    PieceWisePhaseFitter::fitSlope(&data[0], &nu[0], chunkSize, a, b);
    for(size_t ch=0; ch!=leftEdge; ++ch)
      fittedData[ch] = a + b * nu[ch];
    
    for(size_t ch=0; ch!=data.size() - chunkSize; ++ch)
    {
      PieceWisePhaseFitter::fitSlope(&data[ch], &nu[ch], chunkSize, a, b);
      fittedData[ch + chunkSize/2] = a + b * nu[ch + chunkSize/2];
    }
    
    PieceWisePhaseFitter::fitSlope(&data[data.size() - chunkSize], &nu[data.size() - chunkSize], chunkSize, a, b);
    for(size_t ch=rightEdge; ch!=data.size(); ++ch)
      fittedData[ch] = a + b * nu[ch];
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
  * @param chunkSize Size of chunk in number of samples.
  * @param nu Frequencies of the channels (in Hz).
  * @param data The values, same size as @c nu.
  * @param weights Inverse-variance weights.
  * @param fittedData A vector of at least size @c nu.size(), in which the fitted lines are stored.
  */
  void SlidingFit(const double* nu, const std::vector<double>& data, const double* weights, std::vector<double>& fittedData)
  {
    const size_t
      chunkSize = std::min(_chunkSize, data.size()),
      leftEdge = chunkSize/2,
      rightEdge = data.size() - chunkSize + chunkSize/2;
    
    double a, b;
    PieceWisePhaseFitter::fitSlope(&data[0], &nu[0], &weights[0], chunkSize, a, b);
    for(size_t ch=0; ch!=leftEdge; ++ch)
      fittedData[ch] = a + b * nu[ch];
    
    for(size_t ch=0; ch!=data.size() - chunkSize; ++ch)
    {
      PieceWisePhaseFitter::fitSlope(&data[ch], &nu[ch], &weights[ch], chunkSize, a, b);
      fittedData[ch + chunkSize/2] = a + b * nu[ch + chunkSize/2];
    }
    
    PieceWisePhaseFitter::fitSlope(&data[data.size() - chunkSize], &nu[data.size() - chunkSize], &weights[data.size() - chunkSize], chunkSize, a, b);
    for(size_t ch=rightEdge; ch!=data.size(); ++ch)
      fittedData[ch] = a + b * nu[ch];
  }
  
private:
  size_t _chunkSize;
  /**
   * This is used to avoid having to reallocate data every fit
   */
  std::vector<double> _tempArray;
  
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
    std::nth_element(d.data(), d.data()+n/2, d.data()+n);
    double wrapTo = d[n/2];
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
  
  /**
  * Fit a + bx with weights
  */
  void fitSlope(const double* data, const double* frequencies, const double* weights, size_t n, double& a, double& b)
  {
    std::vector<double>& d = _tempArray;
    d.clear();
    for(size_t i=0; i!=n; ++i)
    {
      if(weights[i] > 0.0)
        d.push_back(data[i]);
    }
    std::nth_element(d.data(), d.data()+n/2, d.data()+n);
    double wrapTo = d[n/2];
    for(double& v : d)
    {
      if(v - wrapTo > M_PI)
        v -= 2.0*M_PI;
      else if(v - wrapTo < -M_PI)
        v += 2.0*M_PI;
    }
    std::nth_element(d.data(), d.data()+n/2, d.data()+n);
    wrapTo = d[n/2];
    
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
    mean_y /= sum_weight;
    double mean_x = sum_x / sum_weight;
    double ss_xx = sum_xSq - mean_x * sum_x;
    ss_xy -= sum_x * mean_y;
    b = ss_xy / ss_xx;
    a = mean_y - b * mean_x;
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
