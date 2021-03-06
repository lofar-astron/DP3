// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ScreenConstraint.h"

#include <aocommon/parallelfor.h>

#include <boost/algorithm/string/case_conv.hpp>

#include <iostream>

namespace dp3 {
namespace ddecal {

const double ScreenConstraint::phtoTEC = 1. / 8.4479745e9;
const double ScreenConstraint::TECtoph = 8.4479745e9;
const size_t ScreenConstraint::maxIter = 30;

ScreenConstraint::ScreenConstraint(const common::ParameterSet& parset,
                                   const string& prefix)
    : itsCurrentTime(0), itsIter(0) {
  std::cout << "==========" << (prefix + "order") << "========\n";
  itsBeta = parset.getDouble(prefix + "beta", 5. / 3.);
  itsHeight = parset.getDouble(prefix + "height", 400e3);
  itsOrder = parset.getInt(prefix + "order", 3);
  itsRdiff = parset.getDouble(prefix + "rdiff", 1e3);
  itsMode = boost::to_lower_copy(parset.getString(prefix + "mode", "station"));
  itsAVGMode =
      boost::to_lower_copy(parset.getString(prefix + "average", "tec"));
  itsDebugMode = parset.getInt(prefix + "debug", 0);
}

void ScreenConstraint::Initialize(size_t nAntennas, size_t nDirections,
                                  const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, nDirections, frequencies);
  itsFrequencies = frequencies;
  itsprevsol.assign(NDirections() * NAntennas(), -999.);
  itsAntennaPos.resize(NAntennas());
  itsSourcePos.resize(NDirections());
  itsPiercePoints.resize(NAntennas());
  for (size_t i = 0; i < itsPiercePoints.size(); i++)
    itsPiercePoints[i].resize(NDirections());

  if (itsMode == "station")
    _screenFitters.resize(NAntennas());
  else if (itsMode == "direction")
    _screenFitters.resize(NDirections());
  else if (itsMode == "full")
    _screenFitters.resize(1);
  else if (itsMode == "csfull")
    _screenFitters.resize(NAntennas() - _coreAntennas.size() + 1);
  else
    throw std::runtime_error("Unexpected tecscreen mode: " + itsMode);

  for (size_t i = 0; i != _screenFitters.size(); ++i) {
    _screenFitters[i].setR0(itsRdiff);
    _screenFitters[i].setBeta(itsBeta);
    _screenFitters[i].setOrder(itsOrder);
  }
  if (itsDebugMode > 0) {
    _iterphases.resize(NAntennas() * NDirections() * NChannelBlocks() *
                       maxIter);
  }
}

void ScreenConstraint::setAntennaPositions(
    const std::vector<std::array<double, 3> > antenna_pos) {
  itsAntennaPos = std::move(antenna_pos);
}

void ScreenConstraint::setDirections(
    const std::vector<std::pair<double, double> > source_pos) {
  for (unsigned int i = 0; i < source_pos.size(); i++) {
    itsSourcePos[i].resize(2);
    itsSourcePos[i][0] = source_pos[i].first;
    itsSourcePos[i][1] = source_pos[i].second;
  }
}

void ScreenConstraint::initPiercePoints() {
  for (unsigned int ipos = 0; ipos < itsAntennaPos.size(); ipos++) {
    casacore::MPosition ant(
        casacore::MVPosition(itsAntennaPos[ipos][0], itsAntennaPos[ipos][1],
                             itsAntennaPos[ipos][2]),
        casacore::MPosition::ITRF);
    for (unsigned int isrc = 0; isrc < itsSourcePos.size(); isrc++) {
      casacore::MDirection src(
          casacore::MVDirection(itsSourcePos[isrc][0], itsSourcePos[isrc][1]),
          casacore::MDirection::J2000);

      itsPiercePoints[ipos][isrc] = PiercePoint(ant, src, itsHeight);
    }
  }
}

void ScreenConstraint::setTime(double time) {
  if (itsCurrentTime != time) {
    itsCurrentTime = time;
    itsIter = 0;

    CalculatePiercepoints();

    aocommon::ParallelFor<size_t> loop(NThreads());
    if (itsMode == "station") {
      loop.Run(0, NAntennas(), [&](size_t ipos, size_t /*thread*/) {
        _screenFitters[ipos].calculateCorrMatrix(itsPiercePoints[ipos]);
      });
    } else if (itsMode == "direction") {
      loop.Run(0, NDirections(), [&](size_t idir, size_t /*thread*/) {
        std::vector<PiercePoint*> tmpV(NAntennas());
        for (size_t ipos = 0; ipos < NAntennas(); ipos++)
          tmpV[ipos] = &(itsPiercePoints[ipos][idir]);
        _screenFitters[idir].calculateCorrMatrix(tmpV);
      });
    } else if (itsMode == "full") {
      std::vector<PiercePoint*> tmpV(NAntennas() * NDirections());
      size_t i = 0;
      for (size_t idir = 0; idir < NDirections(); idir++) {
        for (size_t ipos = 0; ipos < NAntennas(); ipos++)
          tmpV[i++] = &(itsPiercePoints[ipos][idir]);
      }
      _screenFitters[0].calculateCorrMatrix(tmpV);
    } else if (itsMode == "csfull") {
      std::vector<PiercePoint*> tmpV(_coreAntennas.size() *
                                     itsSourcePos.size());
      for (size_t iant = 0; iant < _coreAntennas.size(); iant++) {
        size_t ipos = _coreAntennas[iant];
        for (size_t idir = 0; idir < NDirections(); idir++)
          tmpV[iant * NDirections() + idir] = &(itsPiercePoints[ipos][idir]);
      }
      _screenFitters[0].calculateCorrMatrix(tmpV);
      loop.Run(0, _otherAntennas.size(), [&](size_t iant, size_t /*thread*/) {
        size_t ipos = _otherAntennas[iant];
        _screenFitters[iant + 1].calculateCorrMatrix(itsPiercePoints[ipos]);
      });
    } else
      throw std::runtime_error("Unexpected tecscreen mode: " + itsMode);

  } else
    itsIter += 1;
}

void ScreenConstraint::CalculatePiercepoints() {
  casacore::MEpoch time(
      casacore::MVEpoch(itsCurrentTime / (24. * 3600.)));  // convert to MJD
  for (unsigned int i = 0; i < itsPiercePoints.size(); i++)
    for (unsigned int j = 0; j < itsPiercePoints[i].size(); j++)
      itsPiercePoints[i][j].evaluate(time);
}

void ScreenConstraint::getPPValue(
    std::vector<std::vector<DComplex> >& solutions, size_t solutionIndex,
    size_t dirIndex, double& avgTEC, double& error) const {
  if (itsAVGMode == "simple") {
    avgTEC = 0;
    error = 1.0;
    size_t nrch(0);
    for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
      if (isfinite(solutions[ch][dirIndex])) {
        double refphase = std::arg(solutions[ch][dirIndex]);
        // TODO: more advance frequency averaging...
        if (isfinite(solutions[ch][solutionIndex])) {
          avgTEC += std::arg(solutions[ch][solutionIndex] *
                             std::polar<double>(1.0, -1 * refphase)) *
                    itsFrequencies[ch] * phtoTEC;
          nrch++;
        }
      }
    }
    if (nrch > 0) avgTEC /= nrch;
    /*
    double mydelay=0;
    for(size_t ch=0;ch<NChannelBlocks(); ++ch) {
      double refphase=std::arg(solutions[ch][dirIndex]);
      double wavelength=freqtolambda/itsFrequencies[ch]
      //TODO: more advance frequency averaging...

      mydelay +=
    std::arg(solutions[ch][solutionIndex]*std::polar<double>(1.0,-1*refphase))*itsFrequencies[ch]*phtoTEC;
    */
  } else {
    PhaseFitter phfit;
    phfit.Initialize(itsFrequencies);
    double offset = 0.0;
    for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
      if (isfinite(solutions[ch][solutionIndex])) {
        phfit.PhaseData()[ch] = std::arg(solutions[ch][solutionIndex]);
      } else
        phfit.WeightData()[ch] = 0.0;
    }
    if (itsprevsol[solutionIndex] < -100) {
      avgTEC = 0.0;
      phfit.FitTEC2ModelParameters(avgTEC, offset);
      error = phfit.TEC2ModelCost(avgTEC, offset);
    } else {
      avgTEC = itsprevsol[solutionIndex] * TECtoph;
      phfit.FitTEC2ModelParameters(avgTEC, offset);
      error = phfit.TEC2ModelCost(avgTEC, offset);
    }

    avgTEC *= phtoTEC;
  }
  error = 1.0;
}

std::vector<Constraint::Result> ScreenConstraint::Apply(
    std::vector<std::vector<DComplex> >& solutions, double time,
    std::ostream* statStream) {
  // check if we need to reinitialize piercepoints
  setTime(time);
  size_t nrresults = 4;
  if (itsDebugMode > 0) nrresults = 5;
  std::vector<Result> res(nrresults);
  size_t numberofPar = _screenFitters[0].getOrder();
  res[0].vals.resize(_screenFitters.size() * numberofPar);
  res[0].axes = "screennr,par";
  res[0].dims.resize(2);
  res[0].dims[0] = _screenFitters.size();
  res[0].dims[1] = numberofPar;
  res[0].name = "screenpar";
  res[1].vals.resize(NAntennas() * NDirections() * 3);
  res[1].axes = "ant,dir,xyz";
  res[1].dims.resize(3);
  res[1].dims[0] = NAntennas();
  res[1].dims[1] = NDirections();
  res[1].dims[2] = 3;
  res[1].name = "piercepoints";
  res[2].vals.resize(NAntennas() * NDirections());
  res[2].axes = "ant,dir";
  res[2].dims.resize(2);
  res[2].dims[0] = NAntennas();
  res[2].dims[1] = NDirections();
  res[2].name = "TECfitwhite";
  res[3].vals.resize(NAntennas() * NDirections());
  res[3].axes = "ant,dir,freq";
  res[3].dims.resize(3);
  res[3].dims[0] = NAntennas();
  res[3].dims[1] = NDirections();
  res[3].dims[2] = 1;
  res[3].name = "tec";
  if (itsDebugMode > 0) {
    res[4].vals.resize(NAntennas() * NDirections() * NChannelBlocks() *
                       maxIter);
    res[4].axes = "ant,dir,freq,iter";
    res[4].dims.resize(4);
    res[4].dims[0] = NAntennas();
    res[4].dims[1] = NDirections();
    res[4].dims[2] = NChannelBlocks();
    res[4].dims[3] = maxIter;
    res[4].name = "phases";
  }

  // TODOEstimate Weights

  aocommon::ParallelFor<size_t> loop(NThreads());
  loop.Run(0, NAntennas(), [&](size_t antIndex, size_t /*thread*/) {
    int foundantcs = -999;
    int foundantoth = -999;
    if (itsMode == "csfull") {
      for (size_t iant = 0; iant < _coreAntennas.size(); iant++) {
        if (_coreAntennas[iant] == antIndex) {
          foundantcs = iant;
          break;
        }
      }
      if (foundantcs < 0) {
        for (size_t iant = 0; iant < _otherAntennas.size(); iant++) {
          if (_otherAntennas[iant] == antIndex) {
            foundantoth = iant;
            break;
          }
        }
      }
    }
    for (size_t dirIndex = 0; dirIndex < NDirections(); ++dirIndex) {
      double avgTEC = 0.0, error = 1.0;
      size_t solutionIndex = antIndex * NDirections() + dirIndex;
      if (itsDebugMode > 0 and itsIter < maxIter) {
        for (size_t ch = 0; ch < NChannelBlocks(); ch++) {
          // cout<<"writing
          // "<<antIndex<<":"<<dirIndex<<":"<<ch<<":"<<itsIter<<":"<<antIndex*NDirections()*30*NChannelBlocks()+dirIndex*30*NChannelBlocks()+ch*30+itsIter<<"
          // "<<res[4].vals.size()<<","<<solutionIndex<<":"<<solutions[ch].size()<<std::arg(solutions[ch][solutionIndex])<<endl;
          _iterphases[antIndex * NDirections() * maxIter * NChannelBlocks() +
                      dirIndex * maxIter * NChannelBlocks() + ch * maxIter +
                      itsIter] = std::arg(solutions[ch][solutionIndex]);
        }
      }
      getPPValue(solutions, solutionIndex, dirIndex, avgTEC, error);
      if (error <= 0) error = 1;
      if (itsMode == "station") {
        _screenFitters[antIndex].PhaseData()[dirIndex] = avgTEC;
        _screenFitters[antIndex].WData()[dirIndex] = 1. / error;
      } else if (itsMode == "direction") {
        _screenFitters[dirIndex].PhaseData()[antIndex] = avgTEC;
        _screenFitters[dirIndex].WData()[antIndex] = 1. / error;
      } else if (itsMode == "full") {
        _screenFitters[0].PhaseData()[dirIndex * NAntennas() + antIndex] =
            avgTEC;
        _screenFitters[0].WData()[dirIndex * NAntennas() + antIndex] =
            1. / error;
      } else {  // csfull mode
        if (foundantcs >= 0) {
          _screenFitters[0].PhaseData()[foundantcs * NDirections() + dirIndex] =
              avgTEC;
          _screenFitters[0].WData()[foundantcs * NDirections() + dirIndex] =
              1. / error;
        } else if (foundantoth >= 0) {
          _screenFitters[foundantoth + 1].PhaseData()[dirIndex] = avgTEC;
          _screenFitters[foundantoth + 1].WData()[dirIndex] = 1. / error;
        }
      }
    }
  });

  loop.Run(0, _screenFitters.size(), [&](size_t isft, size_t /*thread*/) {
    _screenFitters[isft].doFit();
  });

  loop.Run(0, NAntennas(), [&](size_t antIndex, size_t /*thread*/) {
    int foundantcs = -999;
    int foundantoth = -999;
    if (itsMode == "csfull") {
      for (size_t iant = 0; iant < _coreAntennas.size(); iant++) {
        if (_coreAntennas[iant] == antIndex) {
          foundantcs = iant;
          break;
        }
      }
    }
    if (foundantcs < 0) {
      for (size_t iant = 0; iant < _otherAntennas.size(); iant++) {
        if (_otherAntennas[iant] == antIndex) {
          foundantoth = iant;
          break;
        }
      }
    }
    for (size_t dirIndex = 0; dirIndex < NDirections(); ++dirIndex) {
      size_t solutionIndex = antIndex * NDirections() + dirIndex;
      double avgTEC = 0;
      if (itsMode == "station")
        avgTEC = _screenFitters[antIndex].PhaseData()[dirIndex];
      else if (itsMode == "direction")
        avgTEC = _screenFitters[dirIndex].PhaseData()[antIndex];
      else if (itsMode == "full")
        avgTEC =
            _screenFitters[0].PhaseData()[dirIndex * NAntennas() + antIndex];
      else {  // csfull
        if (foundantcs >= 0)
          avgTEC = _screenFitters[0]
                       .PhaseData()[foundantcs * NDirections() + dirIndex];
        else if (foundantoth >= 0)
          avgTEC = _screenFitters[foundantoth + 1].PhaseData()[dirIndex];
      }

      for (size_t ch = 0; ch < NChannelBlocks(); ++ch)
        solutions[ch][solutionIndex] =
            std::polar<double>(1.0, avgTEC * TECtoph / itsFrequencies[ch]);

      res[3].vals[antIndex * NDirections() + dirIndex] = avgTEC;
      itsprevsol[antIndex * NDirections() + dirIndex] = avgTEC;
      for (size_t i = 0; i < 3; i++) {
        if (itsMode == "station")
          res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
              _screenFitters[antIndex].PPData()[i * NDirections() + dirIndex];
        else if (itsMode == "direction")
          res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
              _screenFitters[dirIndex].PPData()[i * NAntennas() + antIndex];
        else if (itsMode == "full")
          res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
              _screenFitters[0].PPData()[i * NDirections() * NAntennas() +
                                         dirIndex * NAntennas() + antIndex];

        else {  // csfull
          if (foundantcs >= 0)
            res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
                _screenFitters[0]
                    .PPData()[i * _coreAntennas.size() * NDirections() +
                              foundantcs * NDirections() + dirIndex];
          else if (foundantoth >= 0)
            res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
                _screenFitters[foundantoth]
                    .PPData()[i * NDirections() + dirIndex];
        }
      }
    }

    for (size_t dirIndex = 0; dirIndex < NDirections(); ++dirIndex) {
      if (itsMode == "station")
        res[2].vals[antIndex * NDirections() + dirIndex] =
            _screenFitters[antIndex].TECFitWhiteData()[dirIndex];
      else  // not implemented yet for other modes
        res[2].vals[antIndex * NDirections() + dirIndex] = 0;
    }
  });
  for (size_t i = 0; i < _screenFitters.size(); i++)
    for (size_t j = 0; j < numberofPar; j++)
      res[0].vals[i * numberofPar + j] = _screenFitters[i].ParData()[j];

  if (itsDebugMode > 0) res[4].vals = _iterphases;
  return res;
}

}  // namespace ddecal
}  // namespace dp3
