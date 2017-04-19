#ifndef PIERCEPOINT_H
#define PIERCEPOINT_H

#include <vector>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MEpoch.h>

#include <armadillo>

using namespace arma;

class PiercePoint 
{
public:
  PiercePoint();
  PiercePoint(const casacore::MPosition &ant,const casacore::MDirection &source,const double height);
  PiercePoint(const casacore::MPosition &ant,const casacore::MDirection &source);
  void init(const casacore::MPosition &ant,const casacore::MDirection &source,const double height);
  void evaluate(casacore::MEpoch time);
  Col<double>  getValue() const {return itsValue;} 
  casacore::MPosition  getPos() const {return itsPosition;}
  casacore::MDirection  getDir() const {return itsDirection;}
private:
  static const double IONOheight; //default height in meter
  static const double EarthRadius; //default Earth radius in meter
  //station position
  casacore::MPosition     itsPosition;
  //source position
  casacore::MDirection     itsDirection;
  // Ionospheric layer height.
  double              itsIonoHeight;
  //  square of length antenna vector (int ITRF) minus square of vector to piercepoint. This is constant for a assumed spherical Earth
  double              itsC;
  Col<double>         itsValue; //PiercePoint in ITRF coordinates
};

#endif
