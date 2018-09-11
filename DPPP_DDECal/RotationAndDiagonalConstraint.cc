#include "RotationAndDiagonalConstraint.h"
#include "RotationConstraint.h"

#include "../Common/OpenMP.h"

#include <cassert>
#include <cmath>

using namespace std;

namespace LOFAR {

void RotationAndDiagonalConstraint::InitializeDimensions(size_t nAntennas,
                                                         size_t nDirections,
                                                         size_t nChannelBlocks) {
  Constraint::InitializeDimensions(nAntennas, nDirections, nChannelBlocks);

  assert(_nDirections == 1);

  _res.resize(3);
  _res[0].vals.resize(_nAntennas*_nChannelBlocks);
  _res[0].weights.resize(_nAntennas*_nChannelBlocks);
  _res[0].axes="ant,freq";
  _res[0].dims.resize(2);
  _res[0].dims[0]=_nAntennas;
  _res[0].dims[1]=_nChannelBlocks;
  _res[0].name="rotation";

  _res[1].vals.resize(_nAntennas*_nChannelBlocks*2);
  _res[1].weights.resize(_nAntennas*_nChannelBlocks*2);
  _res[1].axes="ant,freq,pol";
  _res[1].dims.resize(3);
  _res[1].dims[0]=_nAntennas;
  _res[1].dims[1]=_nChannelBlocks;
  _res[1].dims[2]=2;
  _res[1].name="amplitude";

  _res[2] = _res[1];
  _res[2].name="phase";
}

void RotationAndDiagonalConstraint::SetWeights(const vector<double>& weights) {
  _res[0].weights = weights;

  // Duplicate weights for two polarizations
  _res[1].weights.resize(_nAntennas*_nChannelBlocks*2);
  size_t indexInWeights = 0;
  for (auto weight: weights) {
    _res[1].weights[indexInWeights++] = weight;
    _res[1].weights[indexInWeights++] = weight;
  }

  _res[2].weights = _res[1].weights;
}

vector<Constraint::Result> RotationAndDiagonalConstraint::Apply(
    vector<vector<dcomplex> >& solutions, double,
    std::ostream* statStream) {
  if (statStream) *statStream<<"["; // begin channel
  double angle0;
  for (uint ch=0; ch<_nChannelBlocks; ++ch) {
    if (statStream) *statStream<<"["; // begin antenna
    for (uint ant=0; ant<_nAntennas; ++ant) {
      // Compute rotation
      complex<double> *data = &(solutions[ch][4*ant]);

      double angle = RotationConstraint::get_rotation(data);
      // Restrict angle between -pi/2 and pi/2
      // Add 2pi to make sure that fmod doesn't see negative numbers
      angle = fmod(angle + 3.5*M_PI, M_PI) - 0.5*M_PI;

      // Right multiply solution with inverse rotation,
      // save only the diagonal
      // Use sin(-phi) == -sin(phi)
      complex<double> a, b;
      a = data[0]*cos(angle) - data[1]*sin(angle);
      b = data[3]*cos(angle) + data[2]*sin(angle);

      // Use station 0 as reference station (for every chanblock), to work
      // around unitary ambiguity
      if (ant==0) {
        angle0 = angle;
        angle = 0.;
      } else {
        angle -= angle0;
        angle = fmod(angle + 3.5*M_PI, M_PI) - 0.5*M_PI;
      }
      _res[0].vals[ant*_nChannelBlocks + ch] = angle;

      _res[1].vals[ant*_nChannelBlocks*2 + 2*ch    ] = abs(a);
      _res[1].vals[ant*_nChannelBlocks*2 + 2*ch + 1] = abs(b);
      _res[2].vals[ant*_nChannelBlocks*2 + 2*ch    ] = arg(a);
      _res[2].vals[ant*_nChannelBlocks*2 + 2*ch  +1] = arg(b);

      // Do the actual constraining
      data[0] =  a * cos(angle);
      data[1] = -a * sin(angle);
      data[2] =  b * sin(angle);
      data[3] =  b * cos(angle);
      if (statStream) *statStream<<"["<<a.real()<<"+"<<a.imag()<<"j,"<<b.real()<<"+"<<b.imag()<<"j,"<<angle<<"]";
      //if (pd) cout<<angle;
      if (statStream && ant<_nAntennas-1) *statStream<<",";
    }
    if (statStream) *statStream<<"]"; // end antenna
    if (statStream && ch<_nChannelBlocks-1) *statStream<<",";
  }
  if (statStream) *statStream<<"]\t"; //end channel

  return _res;
}

} //namespace LOFAR
