#include "PiercePoint.h"

using namespace arma;

namespace DP3 {

  const double PiercePoint::IONOheight = 300000.;
  const double PiercePoint::EarthRadius = 6371000.;

PiercePoint::PiercePoint(double height):
  itsValue(3)
{
  casacore::MPosition ant;//ITRF
  casacore::MDirection source;//J2000 pole
  init(ant,source,height);
}

PiercePoint::PiercePoint(const casacore::MPosition &ant,const casacore::MDirection &source,const double height):
  itsValue(3)
{
  init(ant,source,height);
};

PiercePoint::PiercePoint(const casacore::MPosition &ant,const casacore::MDirection &source):
   itsValue(3)
  { 
    init(ant,source,PiercePoint::IONOheight);
};

void  PiercePoint::init(const casacore::MPosition &ant,const casacore::MDirection &source,const double height){
  
  itsPosition=casacore::MPosition::Convert(ant,casacore::MPosition::ITRF)();
  itsDirection=source;
  itsIonoHeight=height;
  const casacore::MVPosition &mPosition = itsPosition.getValue();
  itsC = mPosition(0)*mPosition(0)+mPosition(1)*mPosition(1)+mPosition(2)*mPosition(2)-
    (itsIonoHeight+PiercePoint::EarthRadius)*(itsIonoHeight+PiercePoint::EarthRadius);
  
}

void  PiercePoint::evaluate(casacore::MEpoch time){
  //Convert direction to ITRF vector
  casacore::MeasFrame myframe(itsPosition,time);
  casacore::MDirection::Ref myref(casacore::MDirection::ITRF,myframe);
  const casacore::MDirection dir = casacore::MDirection::Convert(itsDirection,myref)();
  const casacore::MVDirection &mDir = dir.getValue();
  const casacore::MVPosition &mPos = itsPosition.getValue();
  double A = mDir(0)*mDir(0)+mDir(1)*mDir(1)+mDir(2)*mDir(2);
  double B = mDir(0)*mPos(0) + mDir(1)*mPos(1) +mDir(2)*mPos(2);
  double alpha = (-B + sqrt(B*B - A*itsC))/A;
  for(uword i=0;i<3;i++)
    itsValue(i) = mPos(i) + alpha*mDir(i);
};
}
