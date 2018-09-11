#include "RotationConstraint.h"

#include "../Common/OpenMP.h"

#include <cmath>
#include <assert.h>

using namespace std;

namespace LOFAR {

  void RotationConstraint::InitializeDimensions(size_t nAntennas,
                                                size_t nDirections,
                                                size_t nChannelBlocks) {
  Constraint::InitializeDimensions(nAntennas, nDirections, nChannelBlocks);

  assert(_nDirections == 1);

  _res.resize(1);
  _res[0].vals.resize(_nAntennas*_nChannelBlocks);
  _res[0].axes="ant,freq";
  _res[0].dims.resize(2);
  _res[0].dims[0]=_nAntennas;
  _res[0].dims[1]=_nChannelBlocks;
  _res[0].name="rotation";
}

void RotationConstraint::SetWeights(const vector<double>& weights) {
  _res[0].weights = weights;
}

double RotationConstraint::get_rotation(std::complex<double>* data) {
  // Convert to circular
  complex<double> ll, rr;
  complex<double> i(0,1.);

  ll = data[0] + data[3] - i*data[1] + i*data[2];
  rr = data[0] + data[3] + i*data[1] - i*data[2];
  double angle = 0.5 * (arg(ll) - arg(rr));

  return angle;
}

vector<Constraint::Result> RotationConstraint::Apply(
    vector<vector<dcomplex> >& solutions, double,
    std::ostream* statStream) {
  // Convert to circular
  complex<double> ll, rr;
  complex<double> i(0,1.);

  for (uint ch=0; ch<_nChannelBlocks; ++ch) {
    for (uint ant=0; ant<_nAntennas; ++ant) {
      // Compute rotation
      complex<double> *data= &(solutions[ch][4*ant]);
      double angle = get_rotation(data);
      _res[0].vals[ant*_nChannelBlocks+ch] = angle;

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
