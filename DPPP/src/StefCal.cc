//# GainCal.cc: DPPP step class to do a gain calibration
//# Copyright (C) 2013
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/StefCal.h>

#include <vector>
#include <algorithm>

#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    StefCal::StefCal(uint solInt, uint nChan, string mode, uint maxAntennas)
    : _solInt (solInt),
      _nChan  (nChan),
      _mode   (mode)
    {
      _antMap.resize(maxAntennas, -1);
      init();
    }

    void StefCal::init() {
      dg=1.0e29;
      dgx=1.0e30;
      dgs.clear();

      nSt=vis.shape()(2);

      if (nSt==0) {
        converged=true;
      } else {
        converged=false;
      }
      if (_mode=="fulljones") {
        nUn=nSt;
        nCr=4;
        nSp=1;
      } else if (_mode=="scalarphase") {
        nUn=nSt;
        nCr=1;
        nSp=2;
      } else { // mode=="phaseonly", mode=="diagonal"
        nUn=nSt*2;
        nCr=1;
        nSp=1;
      }

      g.resize(nUn,nCr);
      gold.resize(nUn,nCr);
      gx.resize(nUn,nCr);
      gxx.resize(nUn,nCr);
      h.resize(nUn,nCr);
      z.resize(nUn*_nChan*_solInt*nSp,nCr);

      // Initialize all vectors
      double fronormvis=0;
      double fronormmod=0;

      DComplex* t_vis_p=vis.data();
      DComplex* t_mvis_p=mvis.data();

      uint vissize=vis.size();
      for (uint i=0;i<vissize;++i) {
        fronormvis+=norm(t_vis_p[i]);
        fronormmod+=norm(t_mvis_p[i]);
      }

      fronormvis=sqrt(fronormvis);
      fronormmod=sqrt(fronormmod);

      double ginit=1;
      if (nSt>0 && abs(fronormmod)>1.e-15) {
        ginit=sqrt(fronormvis/fronormmod);
      }
      if (nCr==4) {
        for (uint st=0;st<nUn;++st) {
            g(st,0)=ginit;
            g(st,1)=0.;
            g(st,2)=0.;
            g(st,3)=ginit;
        }
      } else {
        g=ginit;
      }

      gx = g;
    }

    void StefCal::doStep_polarized() {
      gold = g;

      for (uint st=0;st<nSt;++st) {
        h(st,0)=conj(g(st,0));
        h(st,1)=conj(g(st,1));
        h(st,2)=conj(g(st,2));
        h(st,3)=conj(g(st,3));
      }

      for (uint st1=0;st1<nSt;++st1) {
        DComplex* vis_p;
        DComplex* mvis_p;
        Vector<DComplex> w(nCr);
        Vector<DComplex> t(nCr);

        for (uint time=0;time<_solInt;++time) {
          for (uint ch=0;ch<_nChan;++ch) {
            uint zoff=nSt*ch+nSt*_nChan*time;
            mvis_p=&mvis(IPosition(6,0,0,time,ch,0,st1)); for (uint st2=0;st2<nSt;++st2) { z(st2+zoff,0)  = h(st2,0) * mvis_p[st2]; } // itsMVis(IPosition(6,st2,0,time,ch,0,st1))
            mvis_p=&mvis(IPosition(6,0,1,time,ch,0,st1)); for (uint st2=0;st2<nSt;++st2) { z(st2+zoff,0) += h(st2,2) * mvis_p[st2]; } // itsMVis(IPosition(6,st2,0,time,ch,1,st1))
            mvis_p=&mvis(IPosition(6,0,0,time,ch,1,st1)); for (uint st2=0;st2<nSt;++st2) { z(st2+zoff,1)  = h(st2,0) * mvis_p[st2]; } // itsMVis(IPosition(6,st2,1,time,ch,0,st1))
            mvis_p=&mvis(IPosition(6,0,1,time,ch,1,st1)); for (uint st2=0;st2<nSt;++st2) { z(st2+zoff,1) += h(st2,2) * mvis_p[st2]; } // itsMVis(IPosition(6,st2,1,time,ch,1,st1))
            mvis_p=&mvis(IPosition(6,0,0,time,ch,0,st1)); for (uint st2=0;st2<nSt;++st2) { z(st2+zoff,2)  = h(st2,1) * mvis_p[st2]; } // itsMVis(IPosition(6,st2,0,time,ch,0,st1))
            mvis_p=&mvis(IPosition(6,0,1,time,ch,0,st1)); for (uint st2=0;st2<nSt;++st2) { z(st2+zoff,2) += h(st2,3) * mvis_p[st2]; } // itsMVis(IPosition(6,st2,0,time,ch,1,st1))
            mvis_p=&mvis(IPosition(6,0,0,time,ch,1,st1)); for (uint st2=0;st2<nSt;++st2) { z(st2+zoff,3)  = h(st2,1) * mvis_p[st2]; } // itsMVis(IPosition(6,st2,1,time,ch,0,st1))
            mvis_p=&mvis(IPosition(6,0,1,time,ch,1,st1)); for (uint st2=0;st2<nSt;++st2) { z(st2+zoff,3) += h(st2,3) * mvis_p[st2]; } // itsMVis(IPosition(6,st2,1,time,ch,1,st1))
          }
        }

        w=0;

        for (uint time=0;time<_solInt;++time) {
          for (uint ch=0;ch<_nChan;++ch) {
            for (uint st2=0;st2<nSt;++st2) {
              uint zoff=st2+nSt*ch+nSt*_nChan*time;
              w(0) += conj(z(zoff,0))*z(zoff,0) + conj(z(zoff,2))*z(zoff,2);
              w(1) += conj(z(zoff,0))*z(zoff,1) + conj(z(zoff,2))*z(zoff,3);
              w(3) += conj(z(zoff,1))*z(zoff,1) + conj(z(zoff,3))*z(zoff,3);
            }
          }
        }
        w(2)=conj(w(1));

        t=0;

        for (uint time=0;time<_solInt;++time) {
          for (uint ch=0;ch<_nChan;++ch) {
            vis_p=&vis(IPosition(6,0,0,time,ch,0,st1)); for (uint st2=0;st2<nSt;++st2) { t(0) += conj(z(st2+nSt*ch+nSt*_nChan*time,0)) * vis_p[st2]; }// itsVis(IPosition(6,st2,0,time,ch,0,st1))
            vis_p=&vis(IPosition(6,0,1,time,ch,0,st1)); for (uint st2=0;st2<nSt;++st2) { t(0) += conj(z(st2+nSt*ch+nSt*_nChan*time,2)) * vis_p[st2]; }// itsVis(IPosition(6,st2,0,time,ch,1,st1))
            vis_p=&vis(IPosition(6,0,0,time,ch,1,st1)); for (uint st2=0;st2<nSt;++st2) { t(1) += conj(z(st2+nSt*ch+nSt*_nChan*time,0)) * vis_p[st2]; }// itsVis(IPosition(6,st2,1,time,ch,0,st1))
            vis_p=&vis(IPosition(6,0,1,time,ch,1,st1)); for (uint st2=0;st2<nSt;++st2) { t(1) += conj(z(st2+nSt*ch+nSt*_nChan*time,2)) * vis_p[st2]; }// itsVis(IPosition(6,st2,1,time,ch,1,st1))
            vis_p=&vis(IPosition(6,0,0,time,ch,0,st1)); for (uint st2=0;st2<nSt;++st2) { t(2) += conj(z(st2+nSt*ch+nSt*_nChan*time,1)) * vis_p[st2]; }// itsVis(IPosition(6,st2,0,time,ch,0,st1))
            vis_p=&vis(IPosition(6,0,1,time,ch,0,st1)); for (uint st2=0;st2<nSt;++st2) { t(2) += conj(z(st2+nSt*ch+nSt*_nChan*time,3)) * vis_p[st2]; }// itsVis(IPosition(6,st2,0,time,ch,1,st1))
            vis_p=&vis(IPosition(6,0,0,time,ch,1,st1)); for (uint st2=0;st2<nSt;++st2) { t(3) += conj(z(st2+nSt*ch+nSt*_nChan*time,1)) * vis_p[st2]; }// itsVis(IPosition(6,st2,1,time,ch,0,st1))
            vis_p=&vis(IPosition(6,0,1,time,ch,1,st1)); for (uint st2=0;st2<nSt;++st2) { t(3) += conj(z(st2+nSt*ch+nSt*_nChan*time,3)) * vis_p[st2]; }// itsVis(IPosition(6,st2,1,time,ch,1,st1))
          }
        }
        DComplex invdet= 1./(w(0) * w (3) - w(1)*w(2));
        g(st1,0) = invdet * ( w(3) * t(0) - w(1) * t(2) );
        g(st1,1) = invdet * ( w(3) * t(1) - w(1) * t(3) );
        g(st1,2) = invdet * ( w(0) * t(2) - w(2) * t(0) );
        g(st1,3) = invdet * ( w(0) * t(3) - w(2) * t(1) );
      }
    }

    void StefCal::doStep_unpolarized(bool phaseOnly) {
      gold=g;

      for (uint st=0;st<nUn;++st) {
        h(st,0)=conj(g(st,0));
      }

      for (uint st1=0;st1<nUn;++st1) {
        DComplex* vis_p;
        DComplex* mvis_p;
        double ww=0; // Same as w, but specifically for pol==false
        DComplex tt=0; // Same as t, but specifically for pol==false

        DComplex* z_p=z.data();
        mvis_p=&mvis(IPosition(6,0,0,0,0,st1/nSt,st1%nSt));
        vis_p = &vis(IPosition(6,0,0,0,0,st1/nSt,st1%nSt));
        for (uint st1pol=0;st1pol<nSp;++st1pol) {
          for (uint ch=0;ch<_nChan;++ch) {
            for (uint time=0;time<_solInt;++time) {
              DComplex* h_p=h.data();
              for (uint st2=0;st2<nUn;++st2) {
                *z_p = h_p[st2] * *mvis_p; //itsMVis(IPosition(6,st2%nSt,st2/nSt,time,ch,st1/nSt,st1%nSt));
                ww+=norm(*z_p);
                tt+=conj(*z_p) * *vis_p; //itsVis(IPosition(6,st2%nSt,st2/nSt,time,ch,st1/nSt,st1%nSt));
                mvis_p++;
                vis_p++;
                z_p++;
              }
              //cout<<"iS.z bij ch="<<ch<<"="<<iS.z<<endl<<"----"<<endl;
            }
          }
        }
        //cout<<"st1="<<st1%nSt<<(st1>=nSt?"y":"x")<<", t="<<tt<<"       ";
        //cout<<", w="<<ww<<"       ";
        g(st1,0)=tt/ww;
        //cout<<", g="<<iS.g(st1,0)<<endl;
        if (phaseOnly) {
          g(st1,0)/=abs(g(st1,0));
        }
      }
    }

    casa::Matrix<casa::DComplex> StefCal::getSolution() {
      casa::Matrix<casa::DComplex> sol;
      sol.resize(_antMap.size(), nCr);

      uint sSt=0; // Index in stefcal numbering
      for (uint st=0; st<_antMap.size(); ++st) {
        for (uint cr=0; cr<nCr; ++cr) {
          if (_antMap[st]==-1) {
            sol(st,cr)=std::numeric_limits<double>::quiet_NaN();
          } else
            sol(st,cr)=g(sSt,cr);
            if (cr==nCr-1) {
              sSt++;
            }
        }
      }

      ASSERT(sSt==g.size()-1);
      return sol;
    }

    StefCal::Status StefCal::relax(uint iter) {
      double f2 = -1.0;
      double f3 = -0.5;
      double f1 = 1 - f2 - f3;
      double f2q = -0.5;
      double f1q = 1 - f2q;
      double omega = 0.5;
      uint nomega = 24;
      double c1 = 0.5;
      double c2 = 1.2;
      double dgxx;
      bool threestep = false;
      int badIters=0;
      int maxBadIters=5;
      uint itsDebugLevel=0;
      bool itsDetectStalling=true;
      double itsTolerance = 1.0e-12;

      int sstep=0;

      if (itsDetectStalling && iter > 20 && dgx-dg <= 5.0e-3*dg) {
      // This iteration did not improve much upon the previous
      // Stalling detection only after 20 iterations, to account for
      // ''startup problems''
        if (itsDebugLevel>3) {
          cout<<"**"<<endl;
        }
        badIters++;
      } else {
        badIters=0;
      }

      if (badIters>=maxBadIters && itsDetectStalling) {
        if (itsDebugLevel>3) {
          cout<<"Detected stall"<<endl;
        }
        return Status::STALLED;
      }

      dgxx = dgx;
      dgx  = dg;

      double fronormdiff=0;
      double fronormg=0;
      for (uint ant=0;ant<nUn;++ant) {
        for (uint cr=0;cr<nCr;++cr) {
          DComplex diff=g(ant,cr)-gold(ant,cr);
          fronormdiff+=abs(diff*diff);
          fronormg+=abs(g(ant,cr)*g(ant,cr));
        }
      }
      fronormdiff=sqrt(fronormdiff);
      fronormg=sqrt(fronormg);

      dg = fronormdiff/fronormg;
      if (itsDebugLevel>1) {
        dgs.push_back(dg);
      }

      if (dg <= itsTolerance) {
        converged = true;
        return Status::CONVERGED;
      }

      if (itsDebugLevel>7) {
        cout<<"Averaged"<<endl;
      }

      for (uint ant=0;ant<nUn;++ant) {
        for (uint cr=0;cr<nCr;++cr) {
          g(ant,cr) = (1-omega) * g(ant,cr) +
                      omega     * gold(ant,cr);
        }
      }

      if (!threestep) {
        threestep = (iter+1 >= nomega) ||
            ( max(dg,max(dgx,dgxx)) <= 1.0e-3 && dg<dgx && dgx<dgxx);
        if (itsDebugLevel>7) {
          cout<<"Threestep="<<boolalpha<<threestep<<endl;
        }
      }

      if (threestep) {
        if (sstep <= 0) {
          if (dg <= c1 * dgx) {
            if (itsDebugLevel>7) {
              cout<<"dg<=c1*dgx"<<endl;
            }
            for (uint ant=0;ant<nUn;++ant) {
              for (uint cr=0;cr<nCr;++cr) {
                g(ant,cr) = f1q * g(ant,cr) +
                            f2q * gx(ant,cr);
              }
            }
          } else if (dg <= dgx) {
            if (itsDebugLevel>7) {
              cout<<"dg<=dgx"<<endl;
            }
            for (uint ant=0;ant<nUn;++ant) {
              for (uint cr=0;cr<nCr;++cr) {
                g(ant,cr) = f1 * g(ant,cr) +
                            f2 * gx(ant,cr) +
                            f3 * gxx(ant,cr);
              }
            }
          } else if (dg <= c2 *dgx) {
            if (itsDebugLevel>7) {
              cout<<"dg<=c2*dgx"<<endl;
            }
            g = gx;
            sstep = 1;
          } else {
            //cout<<"else"<<endl;
            g = gxx;
            sstep = 2;
          }
        } else {
          if (itsDebugLevel>7) {
            cout<<"no sstep"<<endl;
          }
          sstep = sstep - 1;
        }
      }
      gxx = gx;
      gx = g;

      return Status::NOTCONVERGED;
    }
  } //# end namespace
}
