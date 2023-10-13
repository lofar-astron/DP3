// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ScreenConstraint.h"

#include <aocommon/parallelfor.h>

#include <boost/algorithm/string/case_conv.hpp>

#include <array>
#include <iostream>

namespace dp3 {
namespace ddecal {

ScreenConstraint::ScreenConstraint(const common::ParameterSet& parset,
                                   const string& prefix)
    : current_time_(0), iteration_(0) {
  beta_ = parset.getDouble(prefix + "beta", 5. / 3.);
  height_ = parset.getDouble(prefix + "height", 400e3);
  order_ = parset.getInt(prefix + "order", 3);
  r_diff_ = parset.getDouble(prefix + "rdiff", 1e3);
  mode_ = boost::to_lower_copy(parset.getString(prefix + "mode", "station"));
  avg_mode_ = boost::to_lower_copy(parset.getString(prefix + "average", "tec"));
  debug_mode_ = parset.getInt(prefix + "debug", 0);
}

void ScreenConstraint::Initialize(
    size_t nAntennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, solutions_per_direction, frequencies);
  frequencies_ = frequencies;
  previous_solution_.assign(NDirections() * NAntennas(), -999.);

  if (mode_ == "station")
    screen_fitters_.resize(NAntennas());
  else if (mode_ == "direction")
    screen_fitters_.resize(NDirections());
  else if (mode_ == "full")
    screen_fitters_.resize(1);
  else if (mode_ == "csfull")
    screen_fitters_.resize(NAntennas() - core_antennas_.size() + 1);
  else
    throw std::runtime_error("Unexpected tecscreen mode: " + mode_);

  for (size_t i = 0; i != screen_fitters_.size(); ++i) {
    screen_fitters_[i].setR0(r_diff_);
    screen_fitters_[i].setBeta(beta_);
    screen_fitters_[i].setOrder(order_);
  }
  if (debug_mode_ > 0) {
    iter_phases_.resize(NAntennas() * NDirections() * NChannelBlocks() *
                        kMaxIterations);
  }
}

void ScreenConstraint::SetCoreAntennas(const std::set<size_t>& core_antennas) {
  core_antennas_.clear();
  other_antennas_.clear();
  auto set_iter = core_antennas.begin();
  for (size_t ant = 0; ant < NAntennas(); ++ant) {
    if ((set_iter != core_antennas.end()) && (*set_iter == ant)) {
      core_antennas_.push_back(ant);
      ++set_iter;
    } else {
      other_antennas_.push_back(ant);
    }
  }

  if (mode_ == "csfull") {
    screen_fitters_.resize(NAntennas() - core_antennas_.size() + 1);
  }
}

void ScreenConstraint::InitPiercePoints(
    const std::vector<std::array<double, 3>>& antenna_pos,
    const std::vector<base::Direction>& source_directions) {
  assert(antenna_pos.size() == NAntennas());
  assert(source_directions.size() == NDirections());
  pierce_points_.resize(NAntennas());
  for (unsigned int ipos = 0; ipos < antenna_pos.size(); ipos++) {
    pierce_points_[ipos].clear();
    casacore::MPosition ant(
        casacore::MVPosition(antenna_pos[ipos][0], antenna_pos[ipos][1],
                             antenna_pos[ipos][2]),
        casacore::MPosition::ITRF);
    for (const auto& direction : source_directions) {
      casacore::MDirection src(
          casacore::MVDirection(direction.ra, direction.dec),
          casacore::MDirection::J2000);
      pierce_points_[ipos].emplace_back(ant, src, height_);
    }
  }
}

void ScreenConstraint::SetTime(double time) {
  if (current_time_ != time) {
    current_time_ = time;
    iteration_ = 0;

    CalculatePiercepoints();

    // This DynamicFor can be replaced by a StaticFor.
    // Since the ScreenFitter is unmaintained, it's not done yet.
    aocommon::DynamicFor<size_t> loop(NThreads());
    if (mode_ == "station") {
      loop.Run(0, NAntennas(), [&](size_t ipos, size_t /*thread*/) {
        screen_fitters_[ipos].calculateCorrMatrix(pierce_points_[ipos]);
      });
    } else if (mode_ == "direction") {
      loop.Run(0, NDirections(), [&](size_t idir, size_t /*thread*/) {
        std::vector<PiercePoint*> tmpV(NAntennas());
        for (size_t ipos = 0; ipos < NAntennas(); ipos++)
          tmpV[ipos] = &(pierce_points_[ipos][idir]);
        screen_fitters_[idir].calculateCorrMatrix(tmpV);
      });
    } else if (mode_ == "full") {
      std::vector<PiercePoint*> tmpV(NAntennas() * NDirections());
      size_t i = 0;
      for (size_t idir = 0; idir < NDirections(); idir++) {
        for (size_t ipos = 0; ipos < NAntennas(); ipos++)
          tmpV[i++] = &(pierce_points_[ipos][idir]);
      }
      screen_fitters_[0].calculateCorrMatrix(tmpV);
    } else if (mode_ == "csfull") {
      std::vector<PiercePoint*> tmpV(core_antennas_.size() * NDirections());
      for (size_t iant = 0; iant < core_antennas_.size(); iant++) {
        size_t ipos = core_antennas_[iant];
        for (size_t idir = 0; idir < NDirections(); idir++)
          tmpV[iant * NDirections() + idir] = &(pierce_points_[ipos][idir]);
      }
      screen_fitters_[0].calculateCorrMatrix(tmpV);
      loop.Run(0, other_antennas_.size(), [&](size_t iant, size_t /*thread*/) {
        size_t ipos = other_antennas_[iant];
        screen_fitters_[iant + 1].calculateCorrMatrix(pierce_points_[ipos]);
      });
    } else
      throw std::runtime_error("Unexpected tecscreen mode: " + mode_);

  } else
    iteration_ += 1;
}

void ScreenConstraint::CalculatePiercepoints() {
  casacore::MEpoch time(
      casacore::MVEpoch(current_time_ / (24. * 3600.)));  // convert to MJD
  for (unsigned int i = 0; i < pierce_points_.size(); i++)
    for (unsigned int j = 0; j < pierce_points_[i].size(); j++)
      pierce_points_[i][j].evaluate(time);
}

void ScreenConstraint::GetPpValue(const SolutionSpan& solutions,
                                  size_t antIndex, size_t dirIndex,
                                  double& avgTEC, double& error) const {
  assert(solutions.shape(3) == 1);  // phase only solutions
  if (avg_mode_ == "simple") {
    avgTEC = 0;
    error = 1.0;
    size_t nrch(0);
    for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
      const std::complex<double>& solution0 = solutions(ch, 0, dirIndex, 0);
      if (isfinite(solution0)) {
        double refphase = std::arg(solution0);
        // TODO: more advance frequency averaging...
        const std::complex<double>& antennaSolution =
            solutions(ch, antIndex, dirIndex, 0);
        if (isfinite(antennaSolution)) {
          avgTEC += std::arg(antennaSolution *
                             std::polar<double>(1.0, -1 * refphase)) *
                    frequencies_[ch] * kPhaseToTec;
          nrch++;
        }
      }
    }
    if (nrch > 0) avgTEC /= nrch;
  } else {
    PhaseFitter phfit;
    phfit.Initialize(frequencies_);
    double offset = 0.0;
    for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
      const std::complex<double>& solution =
          solutions(ch, antIndex, dirIndex, 0);
      if (isfinite(solution)) {
        phfit.PhaseData()[ch] = std::arg(solution);
      } else
        phfit.WeightData()[ch] = 0.0;
    }
    const size_t solutionIndex = antIndex * NDirections() + dirIndex;
    if (previous_solution_[solutionIndex] < -100) {
      avgTEC = 0.0;
      phfit.FitTEC2ModelParameters(avgTEC, offset);
      error = phfit.TEC2ModelCost(avgTEC, offset);
    } else {
      avgTEC = previous_solution_[solutionIndex] * kTecToPhase;
      phfit.FitTEC2ModelParameters(avgTEC, offset);
      error = phfit.TEC2ModelCost(avgTEC, offset);
    }

    avgTEC *= kPhaseToTec;
  }
  error = 1.0;
}

std::vector<Constraint::Result> ScreenConstraint::Apply(
    SolutionSpan& solutions, double time,
    [[maybe_unused]] std::ostream* statStream) {
  assert(solutions.shape(3) == 1);  // phase only solutions
  // check if we need to reinitialize piercepoints
  SetTime(time);
  size_t nrresults = 4;
  if (debug_mode_ > 0) nrresults = 5;
  std::vector<Result> res(nrresults);
  size_t numberofPar = screen_fitters_[0].getOrder();
  res[0].vals.resize(screen_fitters_.size() * numberofPar);
  res[0].axes = "screennr,par";
  res[0].dims.resize(2);
  res[0].dims[0] = screen_fitters_.size();
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
  if (debug_mode_ > 0) {
    res[4].vals.resize(NAntennas() * NDirections() * NChannelBlocks() *
                       kMaxIterations);
    res[4].axes = "ant,dir,freq,iter";
    res[4].dims.resize(4);
    res[4].dims[0] = NAntennas();
    res[4].dims[1] = NDirections();
    res[4].dims[2] = NChannelBlocks();
    res[4].dims[3] = kMaxIterations;
    res[4].name = "phases";
  }

  // TODOEstimate Weights

  aocommon::DynamicFor<size_t> loop(NThreads());
  loop.Run(0, NAntennas(), [&](size_t antIndex, size_t /*thread*/) {
    int foundantcs = -999;
    int foundantoth = -999;
    if (mode_ == "csfull") {
      for (size_t iant = 0; iant < core_antennas_.size(); iant++) {
        if (core_antennas_[iant] == antIndex) {
          foundantcs = iant;
          break;
        }
      }
      if (foundantcs < 0) {
        for (size_t iant = 0; iant < other_antennas_.size(); iant++) {
          if (other_antennas_[iant] == antIndex) {
            foundantoth = iant;
            break;
          }
        }
      }
    }
    for (size_t dirIndex = 0; dirIndex < NDirections(); ++dirIndex) {
      double avgTEC = 0.0, error = 1.0;
      if (debug_mode_ > 0 and iteration_ < kMaxIterations) {
        for (size_t ch = 0; ch < NChannelBlocks(); ch++) {
          iter_phases_[antIndex * NDirections() * kMaxIterations *
                           NChannelBlocks() +
                       dirIndex * kMaxIterations * NChannelBlocks() +
                       ch * kMaxIterations + iteration_] =
              std::arg(solutions(ch, antIndex, dirIndex, 0));
        }
      }
      GetPpValue(solutions, antIndex, dirIndex, avgTEC, error);
      if (error <= 0) error = 1;
      if (mode_ == "station") {
        screen_fitters_[antIndex].PhaseData()[dirIndex] = avgTEC;
        screen_fitters_[antIndex].WData()[dirIndex] = 1. / error;
      } else if (mode_ == "direction") {
        screen_fitters_[dirIndex].PhaseData()[antIndex] = avgTEC;
        screen_fitters_[dirIndex].WData()[antIndex] = 1. / error;
      } else if (mode_ == "full") {
        screen_fitters_[0].PhaseData()[dirIndex * NAntennas() + antIndex] =
            avgTEC;
        screen_fitters_[0].WData()[dirIndex * NAntennas() + antIndex] =
            1. / error;
      } else {  // csfull mode
        if (foundantcs >= 0) {
          screen_fitters_[0]
              .PhaseData()[foundantcs * NDirections() + dirIndex] = avgTEC;
          screen_fitters_[0].WData()[foundantcs * NDirections() + dirIndex] =
              1. / error;
        } else if (foundantoth >= 0) {
          screen_fitters_[foundantoth + 1].PhaseData()[dirIndex] = avgTEC;
          screen_fitters_[foundantoth + 1].WData()[dirIndex] = 1. / error;
        }
      }
    }
  });

  loop.Run(0, screen_fitters_.size(), [&](size_t isft, size_t /*thread*/) {
    screen_fitters_[isft].doFit();
  });

  loop.Run(0, NAntennas(), [&](size_t antIndex, size_t /*thread*/) {
    int foundantcs = -999;
    int foundantoth = -999;
    if (mode_ == "csfull") {
      for (size_t iant = 0; iant < core_antennas_.size(); iant++) {
        if (core_antennas_[iant] == antIndex) {
          foundantcs = iant;
          break;
        }
      }
    }
    if (foundantcs < 0) {
      for (size_t iant = 0; iant < other_antennas_.size(); iant++) {
        if (other_antennas_[iant] == antIndex) {
          foundantoth = iant;
          break;
        }
      }
    }
    for (size_t dirIndex = 0; dirIndex < NDirections(); ++dirIndex) {
      double avgTEC = 0;
      if (mode_ == "station")
        avgTEC = screen_fitters_[antIndex].PhaseData()[dirIndex];
      else if (mode_ == "direction")
        avgTEC = screen_fitters_[dirIndex].PhaseData()[antIndex];
      else if (mode_ == "full")
        avgTEC =
            screen_fitters_[0].PhaseData()[dirIndex * NAntennas() + antIndex];
      else {  // csfull
        if (foundantcs >= 0)
          avgTEC = screen_fitters_[0]
                       .PhaseData()[foundantcs * NDirections() + dirIndex];
        else if (foundantoth >= 0)
          avgTEC = screen_fitters_[foundantoth + 1].PhaseData()[dirIndex];
      }

      for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
        solutions(ch, antIndex, dirIndex, 0) =
            std::polar<double>(1.0, avgTEC * kTecToPhase / frequencies_[ch]);
      }

      res[3].vals[antIndex * NDirections() + dirIndex] = avgTEC;
      previous_solution_[antIndex * NDirections() + dirIndex] = avgTEC;
      for (size_t i = 0; i < 3; i++) {
        if (mode_ == "station")
          res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
              screen_fitters_[antIndex].PPData()[i * NDirections() + dirIndex];
        else if (mode_ == "direction")
          res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
              screen_fitters_[dirIndex].PPData()[i * NAntennas() + antIndex];
        else if (mode_ == "full")
          res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
              screen_fitters_[0].PPData()[i * NDirections() * NAntennas() +
                                          dirIndex * NAntennas() + antIndex];

        else {  // csfull
          if (foundantcs >= 0)
            res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
                screen_fitters_[0]
                    .PPData()[i * core_antennas_.size() * NDirections() +
                              foundantcs * NDirections() + dirIndex];
          else if (foundantoth >= 0)
            res[1].vals[antIndex * NDirections() * 3 + dirIndex * 3 + i] =
                screen_fitters_[foundantoth]
                    .PPData()[i * NDirections() + dirIndex];
        }
      }
    }

    for (size_t dirIndex = 0; dirIndex < NDirections(); ++dirIndex) {
      if (mode_ == "station")
        res[2].vals[antIndex * NDirections() + dirIndex] =
            screen_fitters_[antIndex].TECFitWhiteData()[dirIndex];
      else  // not implemented yet for other modes
        res[2].vals[antIndex * NDirections() + dirIndex] = 0;
    }
  });
  for (size_t i = 0; i < screen_fitters_.size(); i++)
    for (size_t j = 0; j < numberofPar; j++)
      res[0].vals[i * numberofPar + j] = screen_fitters_[i].ParData()[j];

  if (debug_mode_ > 0) res[4].vals = iter_phases_;
  return res;
}

}  // namespace ddecal
}  // namespace dp3
