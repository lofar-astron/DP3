// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MOCK_INPUT_H
#define MOCK_INPUT_H

#include <boost/test/unit_test.hpp>
#include "../../../DPInput.h"

namespace DP3 {
namespace DPPP {

class MockInput : public DP3::DPPP::DPInput {
 public:
  MockInput();
  ~MockInput() override;

  void getUVW(const casacore::RefRows&, double, DPBuffer& buffer) override;
  void getWeights(const casacore::RefRows&, DPBuffer& buffer) override;
  void finish() override;
  void show(std::ostream&) const override;
};
}  // namespace DPPP
}  // namespace DP3

#endif