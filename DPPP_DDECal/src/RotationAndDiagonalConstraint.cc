#include <lofar_config.h>

#include <DPPP_DDECal/RotationAndDiagonalConstraint.h>
#include <DPPP_DDECal/RotationConstraint.h>
#include <Common/OpenMP.h>
#include <cmath>
#include <assert.h>

using namespace std;

namespace LOFAR {

void RotationAndDiagonalConstraint::InitializeDimensions(size_t nAntennas,
                                                         size_t nDirections,
                                                         size_t nChannelBlocks) {
  Constraint::InitializeDimensions(nAntennas, nDirections, nChannelBlocks);

  assert(_nDirections == 1);

  _resTemplate.resize(3);
  _resTemplate[0].vals.resize(_nAntennas*_nChannelBlocks);
  _resTemplate[0].weights.resize(_nAntennas*_nChannelBlocks);
  _resTemplate[0].axes="ant,freq";
  _resTemplate[0].dims.resize(2);
  _resTemplate[0].dims[0]=_nAntennas;
  _resTemplate[0].dims[1]=_nChannelBlocks;
  _resTemplate[0].name="rotation";

  _resTemplate[1].vals.resize(_nAntennas*_nChannelBlocks*2);
  _resTemplate[1].weights.resize(_nAntennas*_nChannelBlocks*2);
  _resTemplate[1].axes="ant,freq,pol";
  _resTemplate[1].dims.resize(3);
  _resTemplate[1].dims[0]=_nAntennas;
  _resTemplate[1].dims[1]=_nChannelBlocks;
  _resTemplate[1].dims[2]=2;
  _resTemplate[1].name="amplitude";

  _resTemplate[2] = _resTemplate[1];
  _resTemplate[2].name="phase";
}

void RotationAndDiagonalConstraint::SetWeights(vector<double>& weights) {
  _resTemplate[0].weights.assign(weights.begin(), weights.end());

  // Duplicate weights for two polarizations
  _resTemplate[1].weights.resize(_nAntennas*_nChannelBlocks*2);
  size_t indexInWeights = 0;
  for (auto weight: weights) {
    _resTemplate[1].weights[indexInWeights++] = weight;
    _resTemplate[1].weights[indexInWeights++] = weight;
  }

  _resTemplate[2].weights = _resTemplate[1].weights;
}

vector<Constraint::Result> RotationAndDiagonalConstraint::Apply(
    vector<vector<dcomplex> >& solutions, double) {
  // Convert to circular
  complex<double> ll, rr;
  complex<double> i(0,1.);

  // Find angle
  for (uint ch=0; ch<_nChannelBlocks; ++ch) {
    for (uint ant=0; ant<_nAntennas; ++ant) {
      // Compute rotation
      complex<double> *data = &(solutions[ch][4*ant]);

      double angle = RotationConstraint::get_rotation(data);
      // Restrict angle between -pi/2 and pi/2
      // Add 2pi to make sure that fmod doesn't see negative numbers
      angle = fmod(angle + 3.5*M_PI, M_PI) - 0.5*M_PI;
      _resTemplate[0].vals[ant*_nChannelBlocks + ch] = angle;
 
      // Right multiply solution with inverse rotation,
      // save only the diagonal
      // Use sin(-phi) == -sin(phi)
      complex<double> a, b;
      a = data[0]*cos(angle) - data[1]*sin(angle);
      b = data[3]*cos(angle) + data[2]*sin(angle);
      _resTemplate[1].vals[ant*_nChannelBlocks*2 + 2*ch    ] = abs(a);
      _resTemplate[1].vals[ant*_nChannelBlocks*2 + 2*ch + 1] = abs(b);
      _resTemplate[2].vals[ant*_nChannelBlocks*2 + 2*ch    ] = arg(a);
      _resTemplate[2].vals[ant*_nChannelBlocks*2 + 2*ch  +1] = arg(b);

      // Do the actual constraining
      data[0] =  a * cos(angle);
      data[1] = -a * sin(angle);
      data[2] =  b * sin(angle);
      data[3] =  b * cos(angle);
    }
  }

  return _resTemplate;
}

} //namespace LOFAR
