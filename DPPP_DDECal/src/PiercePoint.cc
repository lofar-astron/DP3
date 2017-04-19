#include <DPPP_DDECal/PiercePoint.h>
using namespace arma;

const double PiercePoint::IONOheight = 300000.; //default height in meter
const double PiercePoint::EarthRadius = 6371000.; //default Earth radius in meter

PiercePoint::PiercePoint():
  itsValue(3)
{
  casa::MPosition ant;//ITRF
  casa::MDirection source;//J2000 pole
  init(ant,source,PiercePoint::IONOheight);
}

PiercePoint::PiercePoint(const casa::MPosition &ant,const casa::MDirection &source,const double height):
  itsValue(3)
{
  init(ant,source,height);
};

PiercePoint::PiercePoint(const casa::MPosition &ant,const casa::MDirection &source):
   itsValue(3)
  { 
    init(ant,source,PiercePoint::IONOheight);
};

void  PiercePoint::init(const casa::MPosition &ant,const casa::MDirection &source,const double height){
  
  itsPosition=casa::MPosition::Convert(ant,casa::MPosition::ITRF)();
  itsDirection=source;
  itsIonoHeight=height;
  const casa::MVPosition &mPosition = itsPosition.getValue();
  itsC = mPosition(0)*mPosition(0)+mPosition(1)*mPosition(1)+mPosition(2)*mPosition(2)-
    (itsIonoHeight+PiercePoint::EarthRadius)*(itsIonoHeight+PiercePoint::EarthRadius);
  
}

void  PiercePoint::evaluate(casa::MEpoch time){
  //Convert direction to ITRF vector
  casa::MeasFrame myframe(itsPosition,time);
  casa::MDirection::Ref myref(casa::MDirection::ITRF,myframe);
  const casa::MDirection dir = casa::MDirection::Convert(itsDirection,myref)();
  const casa::MVDirection &mDir = dir.getValue();
  const casa::MVPosition &mPos = itsPosition.getValue();
  double A = mDir(0)*mDir(0)+mDir(1)*mDir(1)+mDir(2)*mDir(2);
  double B = mDir(0)*mPos(0) + mDir(1)*mPos(1) +mDir(2)*mPos(2);
  double alpha = (-B + sqrt(B*B - A*itsC))/A;
  for(uword i=0;i<3;i++)
    itsValue(i) = mPos(i) + alpha*mDir(i);
};
