#ifndef SCREEN_CONSTRAINT_H
#define SCREEN_CONSTRAINT_H

#include <DPPP_DDECal/multidirsolver.h>
#include <DPPP_DDECal/PiercePoint.h>
#include <DPPP_DDECal/KLFitter.h>

#include <cmath>
#include <vector>
#include <memory>

class ScreenConstraint : public Constraint
{
public:
  ScreenConstraint();
  virtual void init(size_t nAntennas, size_t nDirections, 
                    size_t nChannelBlocks, const double* frequencies);

  virtual std::vector<Constraint::Result> Apply(std::vector<std::vector<MultiDirSolver::DComplex> >& solutions,double time);
  virtual void CalculatePiercepoints();

  void setAntennaPositions(const std::vector<std::vector<double> > antenna_pos);
  void setDirections(const std::vector<std::pair<double, double> > source_pos);
  void setTime(double time);
  void initPiercePoints();

private:
  size_t _nAntennas, _nDirections, _nChannelBlocks;
  std::vector<std::vector<double> > itsAntennaPos;
  std::vector<std::vector<double> > itsSourcePos;
  std::vector<double>               itsFrequencies;
  // antenna positions
  // source positions
  // measures instance ofzo               
  std::vector<std::vector<PiercePoint> > itsPiercePoints;  //temporary hold calculated piercepoints per antenna
  std::vector<KLFitter> _screenFitters;
  double itsCurrentTime;
};

#endif
