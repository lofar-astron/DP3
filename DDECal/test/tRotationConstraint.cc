#include <casacore/casa/BasicMath/Math.h> // near

#include <vector>
#include <iostream>
#include <complex>

#include "../RotationConstraint.h"
#include "../RotationAndDiagonalConstraint.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace DP3;
using namespace casacore;

BOOST_AUTO_TEST_SUITE(rotationconstraint)

BOOST_AUTO_TEST_CASE( test_rotation )
{
  RotationConstraint constraint;
  constraint.InitializeDimensions(1, 1, 1);

  vector<vector<complex<double> > > onesolution(1);
  onesolution[0].resize(4);
  double pi = 3.1415;
  for (double phi=-pi; phi<pi; phi+=pi/6) {
    //cout<<"test phi = "<<phi<<endl;

    /* Solution is of the form ((a,0),(0,b))*rot(phi)
     with rot(phi) = ((cos(phi),-sin(phi)),(sin(phi),cos(phi)))
    */
    onesolution[0][0] =  cos(phi);
    onesolution[0][1] = -sin(phi);
    onesolution[0][2] =  sin(phi);
    onesolution[0][3] =  cos(phi);
    vector<Constraint::Result> constraint_result;
    constraint_result = constraint.Apply(onesolution, 0., nullptr);

    BOOST_CHECK( constraint_result.size() == 1 );
    BOOST_CHECK( constraint_result[0].axes == "ant,dir,freq" );
    //cout<<" got phi = "<<constraint_result[0].vals[0] <<endl;
    BOOST_CHECK( near(constraint_result[0].vals[0], phi) );
    BOOST_CHECK( constraint_result[0].name == "rotation" );
    BOOST_CHECK( constraint_result[0].dims.size() == 3 );
    BOOST_CHECK( constraint_result[0].dims[0] == 1 );
    BOOST_CHECK( constraint_result[0].dims[1] == 1 );
    BOOST_CHECK( constraint_result[0].dims[2] == 1 );
  }
}

BOOST_AUTO_TEST_CASE( test_rotation_and_diagonal )
{
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
  constraint_result = constraint.Apply(onesolution, 0., nullptr);
  BOOST_CHECK( constraint_result.size() == 3 );
  BOOST_CHECK( constraint_result[0].name == "rotation" );
  BOOST_CHECK( constraint_result[0].axes == "ant,dir,freq" );
  BOOST_CHECK( near(constraint_result[0].vals[0], 0.) );
  BOOST_CHECK( constraint_result[0].dims.size() == 3 );
  BOOST_CHECK( constraint_result[0].dims[0] == 1 );
  BOOST_CHECK( constraint_result[0].dims[1] == 1 );
  BOOST_CHECK( constraint_result[0].dims[2] == 1 );

  BOOST_CHECK( constraint_result[1].name == "amplitude" );
  BOOST_CHECK( constraint_result[1].axes == "ant,dir,freq,pol" );
  BOOST_CHECK( near(constraint_result[1].vals[0], abs(a)) );
  BOOST_CHECK( near(constraint_result[1].vals[1], abs(b)) );
  BOOST_CHECK( constraint_result[1].dims.size() == 4 );
  BOOST_CHECK( constraint_result[1].dims[0] == 1 );
  BOOST_CHECK( constraint_result[1].dims[1] == 1 );
  BOOST_CHECK( constraint_result[1].dims[2] == 1 );
  BOOST_CHECK( constraint_result[1].dims[3] == 2 );

  BOOST_CHECK( constraint_result[2].name == "phase" );
  BOOST_CHECK( constraint_result[2].axes == "ant,dir,freq,pol" );
  BOOST_CHECK( near(constraint_result[2].vals[0], arg(a)) );
  BOOST_CHECK( near(constraint_result[2].vals[1], arg(b)) );
  BOOST_CHECK( constraint_result[2].dims.size() == 4 );
  BOOST_CHECK( constraint_result[2].dims[0] == 1 );
  BOOST_CHECK( constraint_result[2].dims[1] == 1 );
  BOOST_CHECK( constraint_result[2].dims[2] == 1 );
  BOOST_CHECK( constraint_result[2].dims[3] == 2 );
}

BOOST_AUTO_TEST_SUITE_END()
