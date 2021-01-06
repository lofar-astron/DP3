// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_PYDPSTEPIMPL_H
#define DPPP_PYDPSTEPIMPL_H

#include "PyDPStep.h"
#include "../DPPP/DPInput.h"

#include <pybind11/pybind11.h>
#include <ostream>

namespace py = pybind11;

namespace DP3 {
namespace DPPP {

class ostream {
 public:
  ostream(std::ostream& os) : os_(os) {}
  void write(std::string& s) { os_ << s; }

 private:
  std::ostream& os_;
};

class DPStepWrapper : public DPStep {
 public:
  using DPStep::DPStep;
  using DPStep::info;
  using DPStep::updateInfo;

  using DPStep::getNextStep;

  DPStep::ShPtr get_next_step() { return DPStep::ShPtr(getNextStep()); }

  bool process_next_step(const DPBuffer& dpbuffer) {
    return get_next_step()->process(dpbuffer);
  }

  int get_count() { return m_count; }
  void set_input(DPInput* input) { m_input = input; };
  void set_parset(const ParameterSet& parset) { m_parset = parset; };
  void set_name(const string& name) { m_name = name; };

  bool m_fetch_uvw = false;
  bool m_fetch_weights = false;

 protected:
  int m_count = 0;
  DPInput* m_input;
  ParameterSet m_parset;
  string m_name;
  const DPBuffer* m_dpbuffer_in;
  NSTimer m_timer;
};

class PyDPStepImpl : public DPStepWrapper {
 public:
  using DPStepWrapper::DPStepWrapper;

  virtual void show(std::ostream& os) const override;

  virtual void updateInfo(const DPInfo&) override;

  virtual bool process(const DPBuffer&) override;

  // Finish the processing of this step and subsequent steps.
  virtual void finish() override;

  void hold();

 private:
  py::object m_py_object;
};

}  // namespace DPPP
}  // namespace DP3

#endif
