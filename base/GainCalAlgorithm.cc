// GainCalAlgorithm.cc: Perform algorithm for gain calibration
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "GainCalAlgorithm.h"

#include <vector>
#include <algorithm>
#include <limits>

#include <iostream>

using casacore::IPosition;
using casacore::Vector;

namespace dp3 {
namespace base {

GainCalAlgorithm::GainCalAlgorithm(unsigned int solInt, unsigned int nChan,
                                   Mode mode, bool scalar, double tolerance,
                                   unsigned int maxAntennas,
                                   bool detectStalling, unsigned int debugLevel)
    : _nSt(maxAntennas),
      _badIters(0),
      _veryBadIters(0),
      _solInt(solInt),
      _nChan(nChan),
      _mode(mode),
      _scalar(scalar),
      _tolerance(tolerance),
      _totalWeight(0.),
      _detectStalling(detectStalling),
      _debugLevel(debugLevel) {
  resetVis();

  _nSt = maxAntennas;
  if (_mode == FULLJONES) {
    if (_scalar)
      throw std::invalid_argument(
          "FULLJONES mode can't be set together with scalar=true");
    _nCr = 4;
    _nSp = 1;
    _savedNCr = 4;
  } else if (_scalar) {
    _nCr = 1;
    _nSp = 2;
    _savedNCr = 1;
  } else {
    _nCr = 1;
    _nSp = 1;
    _savedNCr = 2;
  }

  _vis.resize(IPosition(6, _nSt, 2, _solInt, _nChan, 2, _nSt));
  _mvis.resize(IPosition(6, _nSt, 2, _solInt, _nChan, 2, _nSt));

  if (_scalar || _mode == FULLJONES) {
    _nUn = _nSt;
  } else {
    _nUn = 2 * _nSt;
  }

  _g.resize(_nUn, _nCr);
  _gold.resize(_nUn, _nCr);
  _gx.resize(_nUn, _nCr);
  _gxx.resize(_nUn, _nCr);
  _h.resize(_nUn, _nCr);
  _z.resize(_nUn * _nChan * _solInt * _nSp, _nCr);

  _stationFlagged.resize(_nSt, false);

  init(true);
}

void GainCalAlgorithm::resetVis() {
  _vis = 0;
  _mvis = 0;
  _totalWeight = 0.;
}

void GainCalAlgorithm::clearStationFlagged() { _stationFlagged = false; }

void GainCalAlgorithm::init(bool initSolutions) {
  _dg = 1.0e29;
  _dgx = 1.0e30;
  _dgs.clear();

  _badIters = 0;
  _veryBadIters = 0;

  if (initSolutions) {
    double ginit = 1.0;
    if (_mode != PHASEONLY) {
      // Initialize solution with sensible amplitudes
      double fronormvis = 0;
      double fronormmod = 0;

      DComplex* t_vis_p = _vis.data();
      DComplex* t_mvis_p = _mvis.data();

      unsigned int vissize = _vis.size();
      for (unsigned int i = 0; i < vissize; ++i) {
        fronormvis += norm(t_vis_p[i]);
        fronormmod += norm(t_mvis_p[i]);
      }

      fronormvis = sqrt(fronormvis);
      fronormmod = sqrt(fronormmod);
      if (std::abs(fronormmod) > 1.e-15) {
        ginit = sqrt(fronormvis / fronormmod);
      } else {
        ginit = 1.0;
      }
    }

    if (_nCr == 4) {
      for (unsigned int ant = 0; ant < _nUn; ++ant) {
        _g(ant, 0) = ginit;
        _g(ant, 1) = 0.;
        _g(ant, 2) = 0.;
        _g(ant, 3) = ginit;
      }
    } else {
      _g = ginit;
    }
  } else {  // Take care of NaNs in solution
    double ginit = 0.;
    bool ginitcomputed = false;
    for (unsigned int ant = 0; ant < _nUn; ++ant) {
      if (!std::isfinite(_g(ant, 0).real())) {
        if (!ginitcomputed && !_stationFlagged[ant % _nSt]) {
          // Avoid calling getAverageUnflaggedSolution for stations that are
          // always flagged
          ginit = getAverageUnflaggedSolution();
          ginitcomputed = true;
          if (ginit == 0) {
            init(true);
            return;
          }
        }
        if (_nCr == 4) {
          _g(ant, 0) = ginit;
          _g(ant, 1) = 0.;
          _g(ant, 2) = 0.;
          _g(ant, 3) = ginit;
        } else {
          _g(ant, 0) = ginit;
        }
      }
    }
  }
}

double GainCalAlgorithm::getAverageUnflaggedSolution() {
  // Get average solution of unflagged antennas only once
  // Unflagged means unflagged in previous time slot, so
  // look at NaNs, don't look at stationFlagged (that's for
  // the current timeslot).
  double total = 0.;
  unsigned int unflaggedstations = 0;
  for (unsigned int ant2 = 0; ant2 < _nUn; ++ant2) {
    if (std::isfinite(_g(ant2, 0).real())) {
      total += std::abs(_g(ant2, 0));
      unflaggedstations++;
      if (_nCr == 4) {
        total += std::abs(_g(ant2, 3));
        unflaggedstations++;
      }
    }
  }
  if (unflaggedstations == 0) {
    return 0.;
  } else {
    return total / unflaggedstations;
  }
}

GainCalAlgorithm::Status GainCalAlgorithm::doStep(unsigned int iter) {
  _gxx = _gx;
  _gx = _g;

  bool allFlagged = true;
  for (unsigned int st1 = 0; st1 < _nSt; ++st1) {
    if (!_stationFlagged[st1]) {
      allFlagged = false;
      break;
    }
  }
  if (allFlagged) {
    return CONVERGED;
  }

  if (_mode == FULLJONES) {
    doStep_polarized();
    doStep_polarized();
    return relax(2 * iter);
  } else {
    doStep_unpolarized();
    doStep_unpolarized();
    return relax(2 * iter);
  }
}

void GainCalAlgorithm::doStep_polarized() {
  _gold = _g;

  for (unsigned int st = 0; st < _nSt; ++st) {
    _h(st, 0) = conj(_g(st, 0));
    _h(st, 1) = conj(_g(st, 1));
    _h(st, 2) = conj(_g(st, 2));
    _h(st, 3) = conj(_g(st, 3));
  }

  for (unsigned int st1 = 0; st1 < _nSt; ++st1) {
    if (_stationFlagged[st1]) {
      continue;
    }

    DComplex* vis_p;
    DComplex* mvis_p;
    Vector<DComplex> w(_nCr);
    Vector<DComplex> t(_nCr);

    for (unsigned int time = 0; time < _solInt; ++time) {
      for (unsigned int ch = 0; ch < _nChan; ++ch) {
        unsigned int zoff = _nSt * ch + _nSt * _nChan * time;
        mvis_p = &_mvis(IPosition(6, 0, 0, time, ch, 0, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          _z(st2 + zoff, 0) = _h(st2, 0) * mvis_p[st2];
        }  // itsMVis(IPosition(6,st2,0,time,ch,0,st1))
        mvis_p = &_mvis(IPosition(6, 0, 1, time, ch, 0, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          _z(st2 + zoff, 0) += _h(st2, 2) * mvis_p[st2];
        }  // itsMVis(IPosition(6,st2,0,time,ch,1,st1))
        mvis_p = &_mvis(IPosition(6, 0, 0, time, ch, 1, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          _z(st2 + zoff, 1) = _h(st2, 0) * mvis_p[st2];
        }  // itsMVis(IPosition(6,st2,1,time,ch,0,st1))
        mvis_p = &_mvis(IPosition(6, 0, 1, time, ch, 1, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          _z(st2 + zoff, 1) += _h(st2, 2) * mvis_p[st2];
        }  // itsMVis(IPosition(6,st2,1,time,ch,1,st1))
        mvis_p = &_mvis(IPosition(6, 0, 0, time, ch, 0, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          _z(st2 + zoff, 2) = _h(st2, 1) * mvis_p[st2];
        }  // itsMVis(IPosition(6,st2,0,time,ch,0,st1))
        mvis_p = &_mvis(IPosition(6, 0, 1, time, ch, 0, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          _z(st2 + zoff, 2) += _h(st2, 3) * mvis_p[st2];
        }  // itsMVis(IPosition(6,st2,0,time,ch,1,st1))
        mvis_p = &_mvis(IPosition(6, 0, 0, time, ch, 1, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          _z(st2 + zoff, 3) = _h(st2, 1) * mvis_p[st2];
        }  // itsMVis(IPosition(6,st2,1,time,ch,0,st1))
        mvis_p = &_mvis(IPosition(6, 0, 1, time, ch, 1, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          _z(st2 + zoff, 3) += _h(st2, 3) * mvis_p[st2];
        }  // itsMVis(IPosition(6,st2,1,time,ch,1,st1))
      }
    }

    w = 0;

    for (unsigned int time = 0; time < _solInt; ++time) {
      for (unsigned int ch = 0; ch < _nChan; ++ch) {
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          unsigned int zoff = st2 + _nSt * ch + _nSt * _nChan * time;
          w(0) +=
              conj(_z(zoff, 0)) * _z(zoff, 0) + conj(_z(zoff, 2)) * _z(zoff, 2);
          w(1) +=
              conj(_z(zoff, 0)) * _z(zoff, 1) + conj(_z(zoff, 2)) * _z(zoff, 3);
          w(3) +=
              conj(_z(zoff, 1)) * _z(zoff, 1) + conj(_z(zoff, 3)) * _z(zoff, 3);
        }
      }
    }
    w(2) = conj(w(1));

    t = 0;

    for (unsigned int time = 0; time < _solInt; ++time) {
      for (unsigned int ch = 0; ch < _nChan; ++ch) {
        vis_p = &_vis(IPosition(6, 0, 0, time, ch, 0, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          t(0) +=
              conj(_z(st2 + _nSt * ch + _nSt * _nChan * time, 0)) * vis_p[st2];
        }  // itsVis(IPosition(6,st2,0,time,ch,0,st1))
        vis_p = &_vis(IPosition(6, 0, 1, time, ch, 0, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          t(0) +=
              conj(_z(st2 + _nSt * ch + _nSt * _nChan * time, 2)) * vis_p[st2];
        }  // itsVis(IPosition(6,st2,0,time,ch,1,st1))
        vis_p = &_vis(IPosition(6, 0, 0, time, ch, 1, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          t(1) +=
              conj(_z(st2 + _nSt * ch + _nSt * _nChan * time, 0)) * vis_p[st2];
        }  // itsVis(IPosition(6,st2,1,time,ch,0,st1))
        vis_p = &_vis(IPosition(6, 0, 1, time, ch, 1, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          t(1) +=
              conj(_z(st2 + _nSt * ch + _nSt * _nChan * time, 2)) * vis_p[st2];
        }  // itsVis(IPosition(6,st2,1,time,ch,1,st1))
        vis_p = &_vis(IPosition(6, 0, 0, time, ch, 0, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          t(2) +=
              conj(_z(st2 + _nSt * ch + _nSt * _nChan * time, 1)) * vis_p[st2];
        }  // itsVis(IPosition(6,st2,0,time,ch,0,st1))
        vis_p = &_vis(IPosition(6, 0, 1, time, ch, 0, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          t(2) +=
              conj(_z(st2 + _nSt * ch + _nSt * _nChan * time, 3)) * vis_p[st2];
        }  // itsVis(IPosition(6,st2,0,time,ch,1,st1))
        vis_p = &_vis(IPosition(6, 0, 0, time, ch, 1, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          t(3) +=
              conj(_z(st2 + _nSt * ch + _nSt * _nChan * time, 1)) * vis_p[st2];
        }  // itsVis(IPosition(6,st2,1,time,ch,0,st1))
        vis_p = &_vis(IPosition(6, 0, 1, time, ch, 1, st1));
        for (unsigned int st2 = 0; st2 < _nSt; ++st2) {
          t(3) +=
              conj(_z(st2 + _nSt * ch + _nSt * _nChan * time, 3)) * vis_p[st2];
        }  // itsVis(IPosition(6,st2,1,time,ch,1,st1))
      }
    }
    DComplex invdet = w(0) * w(3) - w(1) * w(2);
    if (std::abs(invdet) == 0) {
      _stationFlagged[st1] = true;
      _g(st1, 0) = 0;
      continue;
    }
    invdet = 1. / invdet;
    _g(st1, 0) = invdet * (w(3) * t(0) - w(1) * t(2));
    _g(st1, 1) = invdet * (w(3) * t(1) - w(1) * t(3));
    _g(st1, 2) = invdet * (w(0) * t(2) - w(2) * t(0));
    _g(st1, 3) = invdet * (w(0) * t(3) - w(2) * t(1));
  }
}

void GainCalAlgorithm::doStep_unpolarized() {
  _gold = _g;

  for (unsigned int ant = 0; ant < _nUn; ++ant) {
    _h(ant, 0) = conj(_g(ant, 0));
  }

  for (unsigned int st1 = 0; st1 < _nUn; ++st1) {
    if (_stationFlagged[st1 % _nSt]) {
      continue;
    }
    DComplex* vis_p;
    DComplex* mvis_p;
    double ww = 0;    // Same as w, but specifically for pol==false
    DComplex tt = 0;  // Same as t, but specifically for pol==false

    mvis_p = &_mvis(IPosition(6, 0, 0, 0, 0, st1 / _nSt, st1 % _nSt));
    vis_p = &_vis(IPosition(6, 0, 0, 0, 0, st1 / _nSt, st1 % _nSt));
    for (unsigned int st1pol = 0; st1pol < _nSp; ++st1pol) {
      for (unsigned int ch = 0; ch < _nChan; ++ch) {
        for (unsigned int time = 0; time < _solInt; ++time) {
          for (unsigned int st2pol = 0; st2pol < _nSp; ++st2pol) {
            DComplex* h_p = _h.data();
            for (unsigned int st2 = 0; st2 < _nUn; ++st2) {
              DComplex z(
                  h_p[st2] *
                  *mvis_p);  // itsMVis(IPosition(6,st2%nSt,st2/nSt,time,ch,st1/nSt,st1%nSt));
              ww += norm(z);
              tt +=
                  conj(z) *
                  *vis_p;  // itsVis(IPosition(6,st2%nSt,st2/nSt,time,ch,st1/nSt,st1%nSt));
              mvis_p++;
              vis_p++;
            }
          }
        }
      }
    }

    // Flag a station if all baselines are flagged or all data is zero
    if (ww == 0 || std::abs(tt) == 0) {
      _stationFlagged[st1 % _nSt] = true;
      _g(st1, 0) = 0;
      continue;
    }
    _g(st1, 0) = tt / ww;

    // Constrain solutions
    if (_mode == PHASEONLY) {
      if (std::abs(_g(st1, 0)) == 0.0)
        throw std::runtime_error(
            "One of the gains solved to zero in gaincal algorithm");
      _g(st1, 0) /= std::abs(_g(st1, 0));
      if (!isFinite(_g(st1, 0)))
        throw std::runtime_error(
            "One of the gains was found to converge to a non-finite value in "
            "gaincal algorithm");
    } else if (_mode == AMPLITUDEONLY) {
      _g(st1, 0) = std::abs(_g(st1, 0));
    }
  }
}

void GainCalAlgorithm::incrementWeight(float weight) { _totalWeight += weight; }

casacore::Matrix<GainCalAlgorithm::DComplex> GainCalAlgorithm::getSolution(
    bool setNaNs) {
  if (setNaNs) {
    for (unsigned int ant = 0; ant < _nUn; ++ant) {
      if (_stationFlagged[ant % _nSt]) {
        for (unsigned int cr = 0; cr < _nCr; ++cr) {
          _g(ant, cr) = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }
  }

  return _g;
}

GainCalAlgorithm::Status GainCalAlgorithm::relax(unsigned int iter) {
  if (_nSt == 0) {
    return CONVERGED;
  }

  double f2 = -1.0;
  double f3 = -0.5;
  double f1 = 1 - f2 - f3;
  double f2q = -0.5;
  double f1q = 1 - f2q;
  double omega = 0.5;
  unsigned int nomega = 24;
  double c1 = 0.5;
  double c2 = 1.2;
  double dgxx;
  bool threestep = false;
  unsigned int maxBadIters = 3;

  int sstep = 0;

  if (_detectStalling && iter > 3) {
    double improvement = _dgx - _dg;

    if (std::abs(improvement) < 5.0e-2 * _dg) {
      // This iteration did not improve much upon the previous
      // Stalling detection only after 4 iterations, to account for
      // ''startup problems'' (not for tec, where stalling happens very soon)
      if (_debugLevel > 3) {
        std::cout << "**\n";
      }
      _badIters++;
    } else if (improvement < 0) {
      _veryBadIters++;
    } else {
      _badIters = 0;
    }

    if (_badIters >= maxBadIters) {
      if (_debugLevel > 3) {
        std::cout << "Detected stall\n";
      }
      return STALLED;
    } else if (_veryBadIters > maxBadIters) {
      if (_debugLevel > 3) {
        std::cout << "Detected fail\n";
      }
      return STALLED;
    }
  }

  dgxx = _dgx;
  _dgx = _dg;

  double fronormdiff = 0;
  double fronormg = 0;
  for (unsigned int ant = 0; ant < _nUn; ++ant) {
    for (unsigned int cr = 0; cr < _nCr; ++cr) {
      DComplex diff = _g(ant, cr) - _gold(ant, cr);
      fronormdiff += std::abs(diff * diff);
      fronormg += std::abs(_g(ant, cr) * _g(ant, cr));
    }
  }
  fronormdiff = sqrt(fronormdiff);
  fronormg = sqrt(fronormg);

  _dg = fronormdiff / fronormg;
  if (_debugLevel > 2) {
    _dgs.push_back(_dg);
  }

  if (_dg <= _tolerance) {
    return CONVERGED;
  }

  if (_debugLevel > 7) {
    std::cout << "Averaged\n";
  }

  for (unsigned int ant = 0; ant < _nUn; ++ant) {
    for (unsigned int cr = 0; cr < _nCr; ++cr) {
      _g(ant, cr) = (1 - omega) * _g(ant, cr) + omega * _gold(ant, cr);
    }
  }

  if (!threestep) {
    threestep = (iter + 1 >= nomega) ||
                (std::max(_dg, std::max(_dgx, dgxx)) <= 1.0e-3 && _dg < _dgx &&
                 _dgx < dgxx);
    if (_debugLevel > 7) {
      std::cout << "Threestep=" << std::boolalpha << threestep << '\n';
    }
  }

  if (threestep) {
    if (sstep <= 0) {
      if (_dg <= c1 * _dgx) {
        if (_debugLevel > 7) {
          std::cout << "dg<=c1*dgx\n";
        }
        for (unsigned int ant = 0; ant < _nUn; ++ant) {
          for (unsigned int cr = 0; cr < _nCr; ++cr) {
            _g(ant, cr) = f1q * _g(ant, cr) + f2q * _gx(ant, cr);
          }
        }
      } else if (_dg <= _dgx) {
        if (_debugLevel > 7) {
          std::cout << "dg<=dgx\n";
        }
        for (unsigned int ant = 0; ant < _nUn; ++ant) {
          for (unsigned int cr = 0; cr < _nCr; ++cr) {
            _g(ant, cr) =
                f1 * _g(ant, cr) + f2 * _gx(ant, cr) + f3 * _gxx(ant, cr);
          }
        }
      } else if (_dg <= c2 * _dgx) {
        if (_debugLevel > 7) {
          std::cout << "dg<=c2*dgx\n";
        }
        _g = _gx;
        sstep = 1;
      } else {
        // cout<<"else"<<endl;
        _g = _gxx;
        sstep = 2;
      }
    } else {
      if (_debugLevel > 7) {
        std::cout << "no sstep\n";
      }
      sstep = sstep - 1;
    }
  }
  return NOTCONVERGED;
}
}  // namespace base
}  // namespace dp3
