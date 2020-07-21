#include "KLFitter.h"

using namespace arma;

namespace DP3{
KLFitter::KLFitter(double r0,double beta,int order):
  itsOrder(order),
  itsR0(r0),
  itsBeta(beta)
{
  
}



void KLFitter::calculateCorrMatrix(const std::vector<PiercePoint> pp){
  itsPiercePoints.set_size(pp.size(),3);
  _phases.set_size(pp.size());
  _weights=eye<mat>(pp.size(),pp.size()); //TODO, make weights sensible
  Mat<double> Distance=zeros<mat>(pp.size(),pp.size());
  for(size_t i=0; i<pp.size();i++){
    Mat<double> A(pp[i].getValue().memptr(),1,3);
    itsPiercePoints.row(i)=A;
  }
  for(size_t n=0;n<pp.size();n++)
    for(size_t m=0;m<pp.size();m++)
      for(size_t i=0;i<3;i++)
	Distance(n,m)+=pow(itsPiercePoints.col(i)[n]-itsPiercePoints.col(i)[m],2);
  itsCorrMatrix=-pow((Distance / ( itsR0*itsR0 ) ),( itsBeta / 2.0 ))/2.0;
  itsinvC=pinv(itsCorrMatrix);
  mat V,U;
  Col<double> S;
  svd(U,S,V,itsCorrMatrix);
  itsU=U(span::all,span(0,itsOrder));
  itsinvU=inv(itsU.t()*(_weights*itsU)); 
}

void KLFitter::calculateCorrMatrix(const std::vector<PiercePoint*> pp){
  itsPiercePoints.set_size(pp.size(),3);
  _phases.set_size(pp.size());
  _weights=eye<mat>(pp.size(),pp.size()); //TODO, make weights sensible
  Mat<double> Distance=zeros<mat>(pp.size(),pp.size());
  for(size_t i=0; i<pp.size();i++){
    Mat<double> A(pp[i]->getValue().memptr(),1,3);
    itsPiercePoints.row(i)=A;
  }
  for(size_t n=0;n<pp.size();n++)
    for(size_t m=0;m<pp.size();m++)
      for(size_t i=0;i<3;i++)
	Distance(n,m)+=pow(itsPiercePoints.col(i)[n]-itsPiercePoints.col(i)[m],2);
  itsCorrMatrix=-pow((Distance / ( itsR0*itsR0 ) ),( itsBeta / 2.0 ))/2.0;
  itsinvC=pinv(itsCorrMatrix);
  mat V,U;
  Col<double> S;
  svd(U,S,V,itsCorrMatrix);
  itsU=U(span::all,span(0,itsOrder));
  itsinvU=inv(itsU.t()*(_weights*itsU)); 
}


void KLFitter::doFit(){
  Mat<double> A=itsU.t()*(_weights*_phases);
  itsPar=itsinvU* A; 
  itsTECFitWhite=(itsinvC*(itsU*itsPar));

  _phases=itsCorrMatrix*itsTECFitWhite;
}
}
