// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "LLSSolver.h"

#include "QRSolver.h"
#include "SVDSolver.h"
#include "NormalEquationsSolver.h"

#include <boost/make_unique.hpp>
#include <boost/algorithm/string/case_conv.hpp>

namespace dp3 {
namespace ddecal {

std::unique_ptr<LLSSolver> LLSSolver::Make(LLSSolverType lss_type, int m, int n,
                                           int nrhs) {
  switch (lss_type) {
    case LLSSolverType::SVD:
      return boost::make_unique<SVDSolver>(m, n, nrhs);
    case LLSSolverType::QR:
      return boost::make_unique<QRSolver>(m, n, nrhs);
    case LLSSolverType::NORMAL_EQUATIONS:
      return boost::make_unique<NormalEquationsSolver>(m, n, nrhs);
  }
  return nullptr;
}

LLSSolverType LLSSolver::ParseType(const std::string& solver) {
  const std::string solver_lowercase = boost::algorithm::to_lower_copy(solver);

  if (solver_lowercase == "svd") {
    return LLSSolverType::SVD;
  } else if (solver_lowercase == "qr") {
    return LLSSolverType::QR;
  } else if (solver_lowercase == "normalequations") {
    return LLSSolverType::NORMAL_EQUATIONS;
  } else {
    throw std::runtime_error("Unknown least squares solver requested: " +
                             solver);
  }
}

}  // namespace ddecal
}  // namespace dp3
