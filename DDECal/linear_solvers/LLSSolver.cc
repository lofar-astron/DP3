// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "LLSSolver.h"

#include "LSMRSolver.h"
#include "QRSolver.h"
#include "SVDSolver.h"
#include "NormalEquationsSolver.h"

#include <boost/make_unique.hpp>

namespace DP3 {
namespace DPPP {

std::unique_ptr<LLSSolver> LLSSolver::Make(LLSSolverType lss_type, int m, int n,
                                           int nrhs) {
  switch (lss_type) {
    case LLSSolverType::LSMR:
#ifdef USE_LSMR
      return boost::make_unique<LSMRSolver>(m, n, nrhs);
#else
      throw std::runtime_error(
          "LSMR was requested for least-squares solving in DDECal but is not "
          "availabe: LSMR must be explicitly turned on in cmake");
#endif
      break;
    case LLSSolverType::SVD:
      return boost::make_unique<SVDSolver>(m, n, nrhs);
    case LLSSolverType::QR:
      return boost::make_unique<QRSolver>(m, n, nrhs);
    case LLSSolverType::NORMAL_EQUATIONS:
      return boost::make_unique<NormalEquationsSolver>(m, n, nrhs);
  }
  return nullptr;
}

}  // namespace DPPP
}  // namespace DP3
