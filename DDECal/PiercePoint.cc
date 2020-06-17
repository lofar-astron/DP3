// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

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
