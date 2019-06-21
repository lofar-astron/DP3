#include "RotationConstraint.h"

#include <cmath>
#include <assert.h>

using namespace std;

namespace DP3 {

  void RotationConstraint::InitializeDimensions(size_t nAntennas,
                                                size_t nDirections,
                                                size_t nChannelBlocks) {
  Constraint::InitializeDimensions(nAntennas, nDirections, nChannelBlocks);

  if(_nDirections != 1) // TODO directions!
    throw std::runtime_error("RotationConstraint can't handle multiple directions yet");

  _res.resize(1);
  _res[0].vals.resize(_nAntennas*_nChannelBlocks);
  _res[0].axes="ant,dir,freq";
  _res[0].dims.resize(3);
  _res[0].dims[0]=_nAntennas;
  _res[0].dims[1]=_nDirections;
  _res[0].dims[2]=_nChannelBlocks;
  _res[0].name="rotation";
}

void RotationConstraint::SetWeights(const vector<double>& weights) {
  _res[0].weights = weights; // TODO directions!
}

double RotationConstraint::get_rotation(std::complex<double>* data) {
  // Convert to circular
  complex<double> i(0,1.);

  complex<double> ll = data[0] + data[3] - i*data[1] + i*data[2];
  complex<double> rr = data[0] + data[3] + i*data[1] - i*data[2];
  double angle = 0.5 * (arg(ll) - arg(rr));

  return angle;
}

vector<Constraint::Result> RotationConstraint::Apply(
    vector<vector<dcomplex> >& solutions, double,
    std::ostream* /*statStream*/) {

  for (unsigned int ch=0; ch<_nChannelBlocks; ++ch) {
    for (unsigned int ant=0; ant<_nAntennas; ++ant) {
      // Compute rotation
      complex<double> *data= &(solutions[ch][4*ant]);
      double angle = get_rotation(data);
      _res[0].vals[ant*_nChannelBlocks+ch] = angle;  // TODO directions!

      // Constrain the data
      data[0] = cos(angle);
      data[1] = -sin(angle);
      data[2] = -data[1];
      data[3] = data[0];
    }
  }

  return _res;
}

} //namespace LOFAR
