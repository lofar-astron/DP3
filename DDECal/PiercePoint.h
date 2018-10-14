#ifndef PIERCEPOINT_H
#define PIERCEPOINT_H

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MEpoch.h>

#include <armadillo>

#include <vector>

namespace DP3 {
 
class PiercePoint 
{
  //default height in meter (300000)
  static const double IONOheight;
  
  //default Earth radius in meter (6371000)
  static const double EarthRadius;
  
public:
  PiercePoint(double height=PiercePoint::IONOheight);
  PiercePoint(const casacore::MPosition &ant,const casacore::MDirection &source,const double height);
  PiercePoint(const casacore::MPosition &ant,const casacore::MDirection &source);
  
  void init(const casacore::MPosition &ant,const casacore::MDirection &source,const double height);
  
  void evaluate(casacore::MEpoch time);
  arma::Col<double>  getValue() const {return itsValue;} 
  casacore::MPosition  getPos() const {return itsPosition;}
  casacore::MDirection  getDir() const {return itsDirection;}
  
private:
  //station position
  casacore::MPosition     itsPosition;
  //source position
  casacore::MDirection     itsDirection;
  // Ionospheric layer height.
  double              itsIonoHeight;
  //  square of length antenna vector (int ITRF) minus square of vector to piercepoint. This is constant for a assumed spherical Earth
  double              itsC;
  arma::Col<double>         itsValue; //PiercePoint in ITRF coordinates
};

}

#endif
