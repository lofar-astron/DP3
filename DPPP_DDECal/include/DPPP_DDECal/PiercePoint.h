#ifndef PIERCEPOINT_H
#define PIERCEPOINT_H

#include <vector>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MEpoch.h>

#include <armadillo>

using namespace arma;

class PiercePoint 
{
public:
  PiercePoint();
  PiercePoint(const casa::MPosition &ant,const casa::MDirection &source,const double height);
  PiercePoint(const casa::MPosition &ant,const casa::MDirection &source);
  void init(const casa::MPosition &ant,const casa::MDirection &source,const double height);
  void evaluate(casa::MEpoch time);
  Col<double>  getValue() const {return itsValue;} 
  casa::MPosition  getPos() const {return itsPosition;}
  casa::MDirection  getDir() const {return itsDirection;}
private:
  static const double IONOheight; //default height in meter
  static const double EarthRadius; //default Earth radius in meter
  //station position
  casa::MPosition     itsPosition;
  //source position
  casa::MDirection     itsDirection;
  // Ionospheric layer height.
  double              itsIonoHeight;
  //  square of length antenna vector (int ITRF) minus square of vector to piercepoint. This is constant for a assumed spherical Earth
  double              itsC;
  Col<double>         itsValue; //PiercePoint in ITRF coordinates
};

#endif
