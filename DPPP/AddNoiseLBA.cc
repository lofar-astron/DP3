// AddNoiseLBA.h: DPPP step class to add LBA random noise to data
// Copyright (C) 2013
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
//
// $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//
// @author Claudio Gheller, Henrik Edler

#include "AddNoiseLBA.h"

#include <iostream>

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <random>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    AddNoiseLBA::AddNoiseLBA (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix)
    : itsInput(input)
    {
      mode = parset.getInt (prefix+"mode",0);
    }

    AddNoiseLBA::~AddNoiseLBA()
    {}

    void AddNoiseLBA::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();
    }

    void AddNoiseLBA::show (std::ostream& os) const
    {
      os << "AddNoiseLBA " << itsName << endl;
    }

    void AddNoiseLBA::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " AddNoiseLBA " << itsName << endl;
    }

    bool AddNoiseLBA::process (const DPBuffer& buf)
    {
      double nu;    
      double stddev;
      double sefd;
      itsTimer.start();


      // Name of the column to add the noise (at the moment not used, just a placeholder)
      string column = "DATA";
      DPBuffer itsBuf;
      itsBuf.copy(buf);
      //////itsBuf.getData().assign (buf.getData());
      //itsBuf.setData(buf.getData());

      //Array<Complex>::const_contiter indIter = itsBuf.getData().cbegin();
      Array<Complex>::contiter indIter = itsBuf.getData().cbegin();

      // Set the exposure
      double exposure = buf.getExposure();

      // Load the Antenna columns
      Vector<Int> antenna1 = info().getAnt1();
      Vector<Int> antenna2 = info().getAnt2();
      //cout << endl;
      //cout << "ANTENNA1 " << antenna1.size() << endl;
      //cout << "ANTENNA2 " << antenna2.size() << endl;

      // Set Number of baselines
      int n_baselines = antenna1.size();
      //cout << "Number of baselines = " << n_baselines << endl;

      // Set the number of correlations
      int n_corr = info().ncorr();
      //cout << "Number of correlations = " << n_corr << endl;

      // Set the LOFAR_ANTENNA_SET
      string antennaSet1 = "LBA_OUTER";
      string antennaSet2 = "LBA_INNER";
      string antennaSet= getInfo().antennaSet();;
      
      // Set the Polynomial coefficients
      double * coeff;
      if (antennaSet.compare(antennaSet1)) coeff = coeffs_outer;
      if (antennaSet.compare(antennaSet2)) coeff = coeffs_inner;

      //int nant = getInfo().antennaNames().size();
      //cout << "N antennas " << nant << endl;

      // Set the number of frequency channels
      int n_freq = getInfo().chanFreqs().size();
      //cout << "Number of frequencies = " << n_freq << endl;
      Vector<Double> chan_freq = getInfo().chanFreqs();
      Vector<Double> chan_width = getInfo().chanWidths();

      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator (seed);

      //casacore::Cube<Complex> noise_cube(n_corr,n_freq,n_baselines);
      //casacore::Cube<Complex> c_noise_cube(n_corr,n_freq,n_baselines);

      DPBuffer runBuf1;
      runBuf1.copy(buf);
      Array<Complex>::contiter run1Iter = runBuf1.getData().cbegin();
      DPBuffer runBuf2;
      runBuf2.copy(buf);
      Array<Complex>::contiter run2Iter = runBuf2.getData().cbegin();
      

      // Add noise

      int ibegin = 0;
      int iend = n_baselines;
      long icount = 0;
      /////////for (int icorr=0; icorr<n_corr; icorr++)
      for (int ibase=ibegin; ibase<iend; ibase++)
      {
        for (int ifreq=0; ifreq<n_freq; ifreq++) 
        {
          nu = chan_freq(ifreq);
	  sefd = coeff[0]+
		 coeff[1]*nu+
		 coeff[2]*pow(nu,2.0)+
		 coeff[3]*pow(nu,3.0)+
		 coeff[4]*pow(nu,4.0)+
		 coeff[5]*pow(nu,5.0);
          stddev = eta * sefd;
	  stddev = stddev / sqrt(2.0*exposure*chan_width[ifreq]);
	  //cout << "STDDEV ======== " << stddev << " " << nu << endl;
          std::normal_distribution<double> distribution(0.0,stddev);

	  //////////for (int ibase=ibegin; ibase<iend; ibase++)
          for (int icorr=0; icorr<n_corr; icorr++)
	  { 
	      double noise_real = distribution(generator);
	      double noise_img  = distribution(generator);
              std::complex<float> c_noise((float)noise_real, (float)noise_img);
	      // update this with direct bracketing
	      *run1Iter = c_noise;
	      *run2Iter = *indIter + c_noise;
	      run1Iter++;
	      run2Iter++;
	      //noise_cube.at(icorr,ifreq,ibase) = c_noise;
	      //c_noise_cube.at(icorr,ifreq,ibase) = *indIter;// + c_noise;
	      indIter++;
	      icount++;
          }
	}  

      }

      if (mode == 0)
      {
	      //itsBuf.setData(noise_cube);
	      itsBuf.setData(runBuf1.getData());
	      /*
              Array<Complex>::contiter outdIter = itsBuf.getData().cbegin();
              Array<Complex>::const_contiter indIter = buf.getData().cbegin();
	      for(int i=0; i<10; i++)
	      {
		      cout << "--- " << *outdIter << "  " << *indIter << endl;
	              outdIter++;
	              indIter++;
	      }
	      */

      } else if (mode == 1)
      {
	      //itsBuf.setData(c_noise_cube);
	      itsBuf.setData(runBuf2.getData());
	      /*
              Array<Complex>::contiter outdIter = itsBuf.getData().cbegin();
              Array<Complex>::const_contiter indIter = buf.getData().cbegin();
              for(int i=0; i<10; i++)
              {
		      cout << "--- " << *outdIter << "  " << *indIter << endl;
	              outdIter++;
	              indIter++;
              }
	      */
      } else
      {
	      cout << "Mode not supported" << endl;
              exit(100);
      }	      

      itsTimer.stop();
      getNextStep()->process(itsBuf);
      return false;
    }


    void AddNoiseLBA::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } // end namespace
}
