#include <DPPP_DDECal/ScreenConstraint.h>
#include <Common/OpenMP.h>


ScreenConstraint::ScreenConstraint() :
  _nAntennas(0),
  _nDirections(0),
  _nChannelBlocks(0),
  itsCurrentTime(0)
{
}

void ScreenConstraint::init(size_t nAntennas, size_t nDirections, size_t nChannelBlocks, const double* frequencies) {
  _nAntennas = nAntennas;
  _nDirections = nDirections;
  _nChannelBlocks = nChannelBlocks;
  itsFrequencies.resize(_nChannelBlocks);
  std::memcpy( itsFrequencies.data(),frequencies, sizeof(double) * _nChannelBlocks);
  itsAntennaPos.resize(_nAntennas);
  itsSourcePos.resize(_nDirections);
  itsPiercePoints.resize(_nAntennas);
  for(uint i=0;i<itsPiercePoints.size();i++)
    itsPiercePoints[i].resize(_nDirections);
  _screenFitters.resize(_nAntennas);

  /*  for(size_t i=0; i!=_screenFitters.size(); ++i)
  {
    _screenFitters[i].SetChannelCount(_nChannelBlocks);
    std::memcpy(_screenFitters[i].FrequencyData(), frequencies, sizeof(double) * _nChannelBlocks);
    _screenFitters[i].SetPPCount(_nDirections);

    }
  */
}

void ScreenConstraint::setAntennaPositions(const std::vector<std::vector<double> > antenna_pos) {
  for(uint i=0;i<antenna_pos.size();i++){
    itsAntennaPos[i].resize(3);
    for(int j=0;j<3;j++)
      itsAntennaPos[i][j]=antenna_pos[i][j];
  }
}

void ScreenConstraint::setDirections(const std::vector<std::pair<double, double> > source_pos) {
  for(uint i=0;i<source_pos.size();i++){
    itsSourcePos[i].resize(2);
    itsSourcePos[i][0]=source_pos[i].first;
    itsSourcePos[i][1]=source_pos[i].second;
  }

}

void ScreenConstraint::initPiercePoints(){
  for(uint ipos=0;ipos<itsAntennaPos.size();ipos++){
    casacore::MPosition ant(casacore::MVPosition(itsAntennaPos[ipos][0],itsAntennaPos[ipos][1],itsAntennaPos[ipos][2]),
			casacore::MPosition::ITRF);
    for(uint isrc=0;isrc<itsSourcePos.size();isrc++){
      casacore::MDirection src(casacore::MVDirection(itsSourcePos[isrc][0],itsSourcePos[isrc][1]),
			   casacore::MDirection::J2000);
	  
      itsPiercePoints[ipos][isrc]=PiercePoint(ant,src);
    }
  }
}

void ScreenConstraint::setTime(double time){
  if (itsCurrentTime!=time){
      
      itsCurrentTime=time;
      CalculatePiercepoints();
#pragma omp parallel for
      for(uint ipos=0;ipos<itsAntennaPos.size();ipos++){
	_screenFitters[ipos].calculateCorrMatrix(itsPiercePoints[ipos]);
      }

    }
}

void ScreenConstraint::CalculatePiercepoints(){
   casacore::MEpoch time(casacore::MVEpoch(itsCurrentTime/(24.*3600.))); //convert to MJD
  for (uint i=0;i<itsPiercePoints.size();i++)
    for (uint j=0;j<itsPiercePoints[i].size();j++)
      itsPiercePoints[i][j].evaluate(time);
}


std::vector<Constraint::Result> ScreenConstraint::Apply(std::vector<std::vector<MultiDirSolver::DComplex> >& solutions,double time) {
  //check if we need to reinitialize piercepoints
  setTime(time);

  std::vector<Result> res(4);
  size_t numberofPar=std::min(_screenFitters[0].getOrder(),_nDirections);
  res[0].vals.resize(_nAntennas*numberofPar);
  res[0].axes="antenna,par";
  res[0].dims.resize(2);
  res[0].dims[0]=_nAntennas;
  res[0].dims[1]=numberofPar;
  res[0].name="screenpar";
  res[1].vals.resize(_nAntennas*_nDirections*3);
  res[1].axes="antenna,dir,xyz";
  res[1].dims.resize(3);
  res[1].dims[0]=_nAntennas;
  res[1].dims[1]=_nDirections;
  res[1].dims[2]=3;
  res[1].name="piercepoints";
  res[2].vals.resize(_nAntennas*_nDirections);
  res[2].axes="antenna,dir";
  res[2].dims.resize(2);
  res[2].dims[0]=_nAntennas;
  res[2].dims[1]=_nDirections;
  res[2].name="TECfitwhite";
  res[3].vals.resize(_nAntennas*_nDirections);
  res[3].axes="antenna,dir";
  res[3].dims.resize(2);
  res[3].dims[0]=_nAntennas;
  res[3].dims[1]=_nDirections;
  res[3].name="phases";

  //TODOEstimate Weights
  

#pragma omp parallel for
  for(size_t antIndex = 0; antIndex<_nAntennas; ++antIndex)
    {
      for(size_t dirIndex = 0; dirIndex<_nDirections; ++dirIndex){
	double avgTEC=0;
	size_t solutionIndex=antIndex*_nDirections+dirIndex;
	for(size_t ch=0;ch<_nChannelBlocks; ++ch){
	  double refphase=std::arg(solutions[ch][dirIndex]);
	  //TODO: more advance frequency averaging...
	  avgTEC += std::arg(solutions[ch][solutionIndex]*std::polar<double>(1.0,-1*refphase))*itsFrequencies[ch];
	}
 	res[3].vals[antIndex*_nDirections+dirIndex]= avgTEC/_nChannelBlocks;
	_screenFitters[antIndex].PhaseData()[dirIndex] = avgTEC/_nChannelBlocks;
      }
      
        _screenFitters[antIndex].doFit();
    
      for(size_t dirIndex = 0; dirIndex<_nDirections; ++dirIndex){
	size_t solutionIndex=antIndex*_nDirections+dirIndex;
	for(size_t ch=0;ch<_nChannelBlocks; ++ch)
	  solutions[ch][solutionIndex] = std::polar<double>(1.0, _screenFitters[antIndex].PhaseData()[dirIndex]/itsFrequencies[ch]);
	//res[3].vals[antIndex*_nDirections+dirIndex]= _screenFitters[antIndex].PhaseData()[dirIndex];
	for (size_t i=0;i<3;i++)
	  res[1].vals[antIndex*_nDirections*3+dirIndex*3+i]= _screenFitters[antIndex].PPData()[i*_nDirections+dirIndex]; 
      }
      for(size_t i=0;i<numberofPar;i++)
	res[0].vals[antIndex*numberofPar+i]= _screenFitters[antIndex].ParData()[i];
      for(size_t dirIndex = 0; dirIndex<_nDirections; ++dirIndex)
	res[2].vals[antIndex*_nDirections+dirIndex]=
	  _screenFitters[antIndex].TECFitWhiteData()[dirIndex];
	
      
    }
   
  
  return res;
}
