#include <DPPP_DDECal/ScreenConstraint.h>
#include <Common/OpenMP.h>

namespace LOFAR{

ScreenConstraint::ScreenConstraint(const ParameterSet& parset,
                      const string& prefix)
 :
  _nAntennas(0),
  _nDirections(0),
  _nChannelBlocks(0),
  itsCurrentTime(0)
{
  cout<<"=========="<<(prefix + "order")<<"========"<<endl;
  itsBeta=parset.getDouble (prefix + "beta", 5./3.);
  itsHeight=parset.getDouble (prefix + "height", 400e3);
  itsOrder=parset.getInt(prefix + "order", 3);
  itsRdiff=parset.getDouble (prefix +"rdiff",1e3);
  itsMode=toLower(parset.getString(prefix+"mode","station") );
  itsAVGMode=toLower(parset.getString(prefix+"average","tec") );
}

void ScreenConstraint::init(size_t nAntennas, size_t nDirections, size_t nChannelBlocks, const double* frequencies) {
  _nAntennas = nAntennas;
  _nDirections = nDirections;
  _nChannelBlocks = nChannelBlocks;
  itsFrequencies.resize(_nChannelBlocks);
  itsprevsol.assign(_nDirections*_nAntennas,-999.);
  std::memcpy( itsFrequencies.data(),frequencies, sizeof(double) * _nChannelBlocks);
  itsAntennaPos.resize(_nAntennas);
  itsSourcePos.resize(_nDirections);
  itsPiercePoints.resize(_nAntennas);
  for(uint i=0;i<itsPiercePoints.size();i++)
    itsPiercePoints[i].resize(_nDirections);
  
  if (itsMode=="station")
    _screenFitters.resize(_nAntennas);
  else if (itsMode=="direction")
    _screenFitters.resize(_nDirections);
  else if (itsMode=="full")
    _screenFitters.resize(1);
  else if (itsMode=="csfull")
    _screenFitters.resize(_nAntennas-_coreAntennas.size()+1);
  else
    THROW (Exception, "Unexpected tecscreen mode: " << itsMode); 
  

  for(size_t i=0; i!=_screenFitters.size(); ++i)
  {
    _screenFitters[i].setR0(itsRdiff);
    _screenFitters[i].setBeta(itsBeta);
    _screenFitters[i].setOrder(itsOrder);
    
  }
  
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
	  
      itsPiercePoints[ipos][isrc]=PiercePoint(ant,src,itsHeight);
    }
  }
}

void ScreenConstraint::setTime(double time){
  if (itsCurrentTime!=time){
      
      itsCurrentTime=time;
      CalculatePiercepoints();
      
      if (itsMode=="station")
#pragma omp parallel for
	for(uint ipos=0;ipos<_nAntennas;ipos++)
	  _screenFitters[ipos].calculateCorrMatrix(itsPiercePoints[ipos]);
      else if (itsMode=="direction")
#pragma omp parallel for
	for(uint idir=0;idir<_nDirections;idir++){
	  std::vector<PiercePoint *> tmpV(_nAntennas);
	  for(uint ipos=0;ipos<_nAntennas;ipos++)
	    tmpV[ipos]=&(itsPiercePoints[ipos][idir]);
	  _screenFitters[idir].calculateCorrMatrix(tmpV);
	}
      else if (itsMode=="full")
	{
	  std::vector<PiercePoint *> tmpV(_nAntennas*_nDirections);
	  size_t i=0;
	  for(uint idir=0;idir<_nDirections;idir++)
	    for(uint ipos=0;ipos<_nAntennas;ipos++)
	      tmpV[i++]=&(itsPiercePoints[ipos][idir]);
	  _screenFitters[0].calculateCorrMatrix(tmpV);
	}
      else if (itsMode=="csfull")
	{
	  std::vector<PiercePoint *> tmpV(_coreAntennas.size()*itsSourcePos.size());
	  for(size_t iant=0; iant<_coreAntennas.size(); iant++){
	    size_t ipos=_coreAntennas[iant];
	    for(uint idir=0;idir<_nDirections;idir++)
	      tmpV[iant*_nDirections+idir]=&(itsPiercePoints[ipos][idir]);
	  }
	  _screenFitters[0].calculateCorrMatrix(tmpV);
#pragma omp parallel for
	  for(size_t iant=0; iant<_otherAntennas.size(); iant++){
	    size_t ipos=_otherAntennas[iant];
	    _screenFitters[iant+1].calculateCorrMatrix(itsPiercePoints[ipos]);
	  }
	}
      else
	THROW (Exception, "Unexpected tecscreen mode: " << itsMode); 
      
  }
  
}


void ScreenConstraint::CalculatePiercepoints(){
   casacore::MEpoch time(casacore::MVEpoch(itsCurrentTime/(24.*3600.))); //convert to MJD
  for (uint i=0;i<itsPiercePoints.size();i++)
    for (uint j=0;j<itsPiercePoints[i].size();j++)
      itsPiercePoints[i][j].evaluate(time);
}



  void  ScreenConstraint::getPPValue(std::vector<std::vector<MultiDirSolver::DComplex> >& solutions,size_t solutionIndex,size_t dirIndex,double &avgTEC) const {
    if (itsAVGMode=="simple"){
      avgTEC=0;
      for(size_t ch=0;ch<_nChannelBlocks; ++ch){
	double refphase=std::arg(solutions[ch][dirIndex]);
	//TODO: more advance frequency averaging...
	avgTEC += std::arg(solutions[ch][solutionIndex]*std::polar<double>(1.0,-1*refphase))*itsFrequencies[ch]*phtoTEC;
      }
      avgTEC/=_nChannelBlocks;
    }
    else{
      PhaseFitter phfit(_nChannelBlocks) ;
      for(size_t ch=0;ch<_nChannelBlocks; ++ch){
	phfit.FrequencyData()[ch]=itsFrequencies.data()[ch];
	phfit.PhaseData()[ch] = std::arg(solutions[ch][solutionIndex]);
      }
      if (itsprevsol[solutionIndex]<-100)
	phfit.FitTEC1ModelParameters(avgTEC);
      else {
	avgTEC=itsprevsol[solutionIndex]*TECtoph;
	phfit.FitTEC1ModelParameters(avgTEC);
      }
	  
      avgTEC*=phtoTEC;
    }

}


std::vector<Constraint::Result> ScreenConstraint::Apply(std::vector<std::vector<MultiDirSolver::DComplex> >& solutions,double time) {
  //check if we need to reinitialize piercepoints
  setTime(time);

  std::vector<Result> res(4);
  size_t numberofPar=_screenFitters[0].getOrder();
  res[0].vals.resize(_screenFitters.size()*numberofPar);
  res[0].axes="screennr,par";
  res[0].dims.resize(2);
  res[0].dims[0]=_screenFitters.size();
  res[0].dims[1]=numberofPar;
  res[0].name="screenpar";
  res[1].vals.resize(_nAntennas*_nDirections*3);
  res[1].axes="ant,dir,xyz";
  res[1].dims.resize(3);
  res[1].dims[0]=_nAntennas;
  res[1].dims[1]=_nDirections;
  res[1].dims[2]=3;
  res[1].name="piercepoints";
  res[2].vals.resize(_nAntennas*_nDirections);
  res[2].axes="ant,dir";
  res[2].dims.resize(2);
  res[2].dims[0]=_nAntennas;
  res[2].dims[1]=_nDirections;
  res[2].name="TECfitwhite";
  res[3].vals.resize(_nAntennas*_nDirections);
  res[3].axes="ant,dir,freq";
  res[3].dims.resize(3);
  res[3].dims[0]=_nAntennas;
  res[3].dims[1]=_nDirections;
  res[3].dims[2]=1;
  res[3].name="tec";

  //TODOEstimate Weights
  

#pragma omp parallel for
  for(size_t antIndex = 0; antIndex<_nAntennas; ++antIndex)
    {

      int foundantcs=-999;
      int foundantoth=-999;
      if (itsMode=="csfull")
	for(size_t iant=0; iant<_coreAntennas.size(); iant++){
	  if (_coreAntennas[iant]==antIndex){
	    foundantcs=iant;
	    break;
	  }
	}
      if (foundantcs<0)
	for(size_t iant=0; iant<_otherAntennas.size(); iant++){
	  if(_otherAntennas[iant]==antIndex){
	    foundantoth=iant;
	    break;
	  }
	}
      for(size_t dirIndex = 0; dirIndex<_nDirections; ++dirIndex){
	double avgTEC;
	size_t solutionIndex=antIndex*_nDirections+dirIndex;
	getPPValue(solutions,solutionIndex,dirIndex,avgTEC);
	if (itsMode=="station")
	  _screenFitters[antIndex].PhaseData()[dirIndex] = avgTEC;
	else if (itsMode=="direction")
	  _screenFitters[dirIndex].PhaseData()[antIndex] = avgTEC;
	else if (itsMode=="full")
	  _screenFitters[0].PhaseData()[dirIndex*_nAntennas+antIndex] = avgTEC;
	else
	  {//csfull mode
	    if (foundantcs>=0){
	      _screenFitters[0].PhaseData()[foundantcs*_nDirections+dirIndex]= avgTEC;
	    }
	    else if (foundantoth>=0){
	      _screenFitters[foundantoth+1].PhaseData()[dirIndex]= avgTEC;
	    }
	  }
      }
    }

#pragma omp parallel for
  for(size_t isft=0;isft<_screenFitters.size();isft++)
    _screenFitters[isft].doFit();
 #pragma omp parallel for
  for(size_t antIndex = 0; antIndex<_nAntennas; ++antIndex)
    { 
      int foundantcs=-999;
      int foundantoth=-999;
      if (itsMode=="csfull")
	for(size_t iant=0; iant<_coreAntennas.size(); iant++){
	  if (_coreAntennas[iant]==antIndex){
	    foundantcs=iant;
	    break;
	  }
	}
      if (foundantcs<0)
	for(size_t iant=0; iant<_otherAntennas.size(); iant++){
	  if(_otherAntennas[iant]==antIndex){
	    foundantoth=iant;
	    break;
	  }
	}
      for(size_t dirIndex = 0; dirIndex<_nDirections; ++dirIndex){
	size_t solutionIndex=antIndex*_nDirections+dirIndex;
	double avgTEC=0;
	if (itsMode=="station")
	  avgTEC=_screenFitters[antIndex].PhaseData()[dirIndex];
	else if (itsMode=="direction")
	  avgTEC=_screenFitters[dirIndex].PhaseData()[antIndex];
	else if (itsMode=="full")
	  avgTEC=_screenFitters[0].PhaseData()[dirIndex*_nAntennas+antIndex];
	else
	  {//csfull
	    if (foundantcs>=0)
	      avgTEC=_screenFitters[0].PhaseData()[foundantcs*_nDirections+dirIndex];
	    else if (foundantoth>=0)
	      avgTEC=_screenFitters[foundantoth+1].PhaseData()[dirIndex];
	  
	  
	  }
	
      
 	

	for(size_t ch=0;ch<_nChannelBlocks; ++ch)
	  solutions[ch][solutionIndex] = std::polar<double>(1.0, avgTEC*TECtoph/itsFrequencies[ch]);
	
	res[3].vals[antIndex*_nDirections+dirIndex]= avgTEC;
	itsprevsol[antIndex*_nDirections+dirIndex]=avgTEC;
	for (size_t i=0;i<3;i++)
	  if (itsMode=="station")
	    res[1].vals[antIndex*_nDirections*3+dirIndex*3+i]= _screenFitters[antIndex].PPData()[i*_nDirections+dirIndex]; 
	  else if (itsMode=="direction")
	    res[1].vals[antIndex*_nDirections*3+dirIndex*3+i]= _screenFitters[dirIndex].PPData()[i*_nAntennas+antIndex]; 
	  else if (itsMode=="full")
	    res[1].vals[antIndex*_nDirections*3+dirIndex*3+i]= _screenFitters[0].PPData()[i*_nDirections*_nAntennas+dirIndex*_nAntennas+antIndex]; 
	    
	  else
	    {//csfull
	      if (foundantcs>=0)
		res[1].vals[antIndex*_nDirections*3+dirIndex*3+i]= _screenFitters[0].PPData()[i*_coreAntennas.size()*_nDirections+foundantcs*_nDirections+dirIndex]; 
	      else if (foundantoth>=0)
		res[1].vals[antIndex*_nDirections*3+dirIndex*3+i]= _screenFitters[foundantoth].PPData()[i*_nDirections+dirIndex]; 
	    
	  
	    }
      }
	    
	
      
      for(size_t dirIndex = 0; dirIndex<_nDirections; ++dirIndex){
	if (itsMode=="station")
	  res[2].vals[antIndex*_nDirections+dirIndex]=
	    _screenFitters[antIndex].TECFitWhiteData()[dirIndex];
	else  //not implemented yet for other modes
	  res[2].vals[antIndex*_nDirections+dirIndex]=0;
      }
    }
  for(size_t i=0;i<_screenFitters.size();i++)
    for(size_t j=0;j<numberofPar;j++)
      res[0].vals[i*numberofPar+j]= _screenFitters[i].ParData()[j];
   
  
  return res;
}
} //namespace LOFAR
