//# GainCal.cc: DPPP step class to SageCal visibilities
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
#include <DPPP_SageCal/SageCal.h>
#include <sagecal/data.h>
#include <sagecal/Common.h>

#include <iostream>
#include <Common/ParameterSet.h>
#include <Common/Timer.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using namespace casacore;

namespace LOFAR {
  namespace DPPP {

    SageCal::SageCal (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix)
    : _input(input)
    {
      _iodata.tilesz = parset.getInt(prefix+"tilesize", Data::TileSize);
      string mode = parset.getString("mode", "replace");
      if (mode=="replace") {
        Data::DoSim = SIMUL_ONLY;
      } else if (mode=="add") {
        Data::DoSim = SIMUL_ADD;
      } else if (mode=="subtract") {
        Data::DoSim = SIMUL_SUB;
      } else {
        THROW(Exception, "Unknown mode: " << mode);
      }

      _skymodelfile = parset.getString(prefix+"skymodelfile");
      _clusterfile = parset.getString(prefix+"clusterfile");

      openblas_set_num_threads(1);
    }

    SageCal::~SageCal()
    {}

    DPStep::ShPtr SageCal::makeStep (DPInput* input,
                                    const ParameterSet& parset,
                                    const std::string& prefix)
    {
      return DPStep::ShPtr(new SageCal(input, parset, prefix));
    }

    void SageCal::readAuxData() {
      /* Do the same as SageCal's readAuxData */
      _iodata.N = info().nantenna();
      _iodata.Nbase = _iodata.N * (_iodata.N - 1) / 2;

      _iodata.deltat = info().timeInterval();
      _iodata.totalt = info().ntime();

      MDirection ref_dir = info().phaseCenter();
      _iodata.ra0 = ref_dir.getValue().get()[0];
      _iodata.dec0 = ref_dir.getValue().get()[1];

      _iodata.Nms = 1;
      try {
        _iodata.u=new double[_iodata.Nbase*_iodata.tilesz];
        _iodata.v=new double[_iodata.Nbase*_iodata.tilesz];
        _iodata.w=new double[_iodata.Nbase*_iodata.tilesz];
        _iodata.x=new double[8*_iodata.Nbase*_iodata.tilesz];
        _iodata.xo=new double[8*_iodata.Nbase*_iodata.tilesz*_iodata.Nchan];
        _iodata.freqs=new double[_iodata.Nchan];
        _iodata.flag=new double[_iodata.Nbase*_iodata.tilesz];
        _iodata.NchanMS=new int[_iodata.Nms];
      } catch (const std::bad_alloc& e) {
        THROW(Exception, "Allocating memory for _iodata failed. Quitting."<< e.what());
      }

      _iodata.NchanMS[0] = info().nchan();
      _iodata.freq0 = 0.0;
      for (uint ci=0; ci<info().nchan(); ++ci) {
        _iodata.freqs[ci] = info().chanFreqs()[ci];
        _iodata.freq0 += _iodata.freqs[ci];
      }
      _iodata.freq0 /= (double)_iodata.Nchan;

      _iodata.deltaf = info().chanWidths()[0];
    }

    void SageCal::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();

      readAuxData();

      // Initialize baselines
      clus_source_t *carr;
      baseline_t *barr;
      if ((barr = (baseline_t*)calloc((size_t)_iodata.Nbase*_iodata.tilesz,
                                      sizeof(baseline_t)))==0) {
        THROW(Exception, "No memory could be allocated for baselines");
       }
       generate_baselines(_iodata.Nbase, _iodata.tilesz, _iodata.N, barr,
                          Data::Nt);

      // Read skymodel
      read_sky_cluster(_skymodelfile.c_str(), _clusterfile.c_str(), &carr,
                       &_num_clusters, _iodata.freq0, _iodata.ra0,
                       _iodata.dec0, Data::format);
      ASSERTSTR(_num_clusters > 0, "No clusters found");

      /* update cluster array with correct pointers to parameters */
      size_t cj=0;

      readAuxData();
      for (int ci=0; ci<_num_clusters; ci++) {
        if ((carr[ci].p=(int*)calloc((size_t)carr[ci].nchunk,sizeof(int)))==0) {
          THROW(Exception, "No memory could be allocated for clusters");
        }
        for (int ck=0; ck<carr[ci].nchunk; ck++) {
          carr[ci].p[ck]=cj*8*_iodata.N;
          cj++;
        }
      }
    }

    void SageCal::show (std::ostream& os) const
    {
      os << "SageCal " << _name << endl;
      os << "  skymodelfile  : " << _skymodelfile << endl;
      os << "  clusterfile   : " << _clusterfile << endl;
      os << "  num_clusters  : " << _num_clusters << endl;
    }

    void SageCal::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, _timer.getElapsed(), duration);
      os << " SageCal " << _name << endl;
    }

    bool SageCal::process (const DPBuffer& bufin)
    {
      _timer.start();
      _buffer.copy (bufin);
      _input->fetchUVW(bufin, _buffer, _timer);
      _input->fetchWeights(bufin, _buffer, _timer);

      // Fill U, V, W
      // loadData(bufin, _iodata, &iodata.fratio)

      precalculate_coherencies(_iodata.u, _iodata.v, _iodata.w, coh,
                               _iodata.N, _iodata.Nbase * _iodata.tilesz,
                               barr, carr, 
                               _num_clusters, _iodata.freq0, _iodata.deltaf,
                               _iodata.deltat,
                               _iodata.dec0, Data::min_uvcut, Data::max_uvcut,
                               Data::Nt);

      _timer.stop();
      getNextStep()->process(_buffer);
      return false;
    }


    void SageCal::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } //# end namespace
}
