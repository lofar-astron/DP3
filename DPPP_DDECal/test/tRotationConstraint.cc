#include <lofar_config.h>
#include <Common/LofarLogger.h>  // ASSERT
#include <casa/BasicMath/Math.h> // near 

#include <vector>
#include <iostream>
#include <complex>

#include <DPPP_DDECal/RotationConstraint.h>
#include <DPPP_DDECal/RotationAndDiagonalConstraint.h>

using namespace std;
using namespace casacore;
using namespace LOFAR;

void test_rotation() {
  RotationConstraint constraint;
  constraint.InitializeDimensions(1, 1, 1);

  vector<vector<complex<double> > > onesolution(1);
  onesolution[0].resize(4);
  double pi = 3.1415;
  for (double phi=-pi; phi<pi; phi+=pi/6) {
    cout<<"test phi = "<<phi<<endl;
    /* Solution is of the form ((a,0),(0,b))*rot(phi)
     with rot(phi) = ((cos(phi),-sin(phi)),(sin(phi),cos(phi)))
    */
    onesolution[0][0] =  cos(phi);
    onesolution[0][1] = -sin(phi);
    onesolution[0][2] =  sin(phi);
    onesolution[0][3] =  cos(phi);

    vector<Constraint::Result> constraint_result;
    constraint_result = constraint.Apply(onesolution, 0.);

    ASSERT( constraint_result.size() == 1 );
    ASSERT( constraint_result[0].axes == "ant,freq" );
    cout<<" got phi = "<<constraint_result[0].vals[0] <<endl;
    ASSERT( near(constraint_result[0].vals[0], phi) );
    ASSERT( constraint_result[0].name == "rotation" );
    ASSERT( constraint_result[0].dims.size() == 2 );
    ASSERT( constraint_result[0].dims[0] == 1 );
    ASSERT( constraint_result[0].dims[1] == 1 );
  }
}

void test_rotation_and_diagonal() {
  RotationAndDiagonalConstraint constraint;
  constraint.InitializeDimensions(1, 1, 1);

  vector<vector<complex<double> > > onesolution(1);
  onesolution[0].resize(4);
  double pi = 3.1415;
  double phi= pi/6;
  const complex<double> i(0,1.);
  complex<double> a=2.*exp(i*0.3), b=3.*exp(i*-0.2);

  /* Solution is of the form ((a,0),(0,b))*rot(phi)
     with rot(phi) = ((cos(phi),-sin(phi)),(sin(phi),cos(phi)))
  */
  onesolution[0][0] = a * cos(phi);
  onesolution[0][1] = a *-sin(phi);
  onesolution[0][2] = b * sin(phi);
  onesolution[0][3] = b * cos(phi);

  vector<Constraint::Result> constraint_result;
  constraint_result = constraint.Apply(onesolution, 0.);

  ASSERT( constraint_result.size() == 3 );
  ASSERT( constraint_result[0].name == "rotation" );
  ASSERT( constraint_result[0].axes == "ant,freq" );
  ASSERT( near(constraint_result[0].vals[0], phi) );
  ASSERT( constraint_result[0].dims.size() == 2 );
  ASSERT( constraint_result[0].dims[0] == 1 );
  ASSERT( constraint_result[0].dims[1] == 1 );

  ASSERT( constraint_result[1].name == "amplitude" );
  ASSERT( constraint_result[1].axes == "ant,freq,pol" );
  ASSERT( near(constraint_result[1].vals[0], abs(a)) );
  ASSERT( near(constraint_result[1].vals[1], abs(b)) );
  ASSERT( constraint_result[1].dims.size() == 3 );
  ASSERT( constraint_result[1].dims[0] == 1 );
  ASSERT( constraint_result[1].dims[1] == 1 );
  ASSERT( constraint_result[1].dims[2] == 2 );

  ASSERT( constraint_result[2].name == "phase" );
  ASSERT( constraint_result[2].axes == "ant,freq,pol" );
  ASSERT( near(constraint_result[2].vals[0], arg(a)) );
  ASSERT( near(constraint_result[2].vals[1], arg(b)) );
  ASSERT( constraint_result[2].dims.size() == 3 );
  ASSERT( constraint_result[2].dims[0] == 1 );
  ASSERT( constraint_result[2].dims[1] == 1 );
  ASSERT( constraint_result[2].dims[2] == 2 );
}

int main(int, char**) {
  test_rotation();
  test_rotation_and_diagonal();

  return 0;
}
