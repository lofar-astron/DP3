// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_PYDPSTEPIMPL_H
#define DPPP_PYDPSTEPIMPL_H

#include "PyStep.h"
#include "../steps/InputStep.h"

#include <memory>
#include <ostream>

// Use a forward declaration and a pointer to pybind11::object in PyDPStepImpl,
// otherwise GCC gives a warning regarding visibility.
// Note: When #including <pybind11/pybind11.h> before this file, the warning
// will come back.
namespace pybind11 {
class object;
}

namespace dp3 {
namespace pythondp3 {

class ostream_wrapper {
 public:
  ostream_wrapper(std::ostream& os) : os_(os) {}
  void write(std::string& s) { os_ << s; }

 private:
  std::ostream& os_;
};

class StepWrapper : public steps::Step {
 public:
  using Step::info;
  using Step::Step;
  using Step::updateInfo;

  using Step::getNextStep;

  Step::ShPtr get_next_step() { return Step::ShPtr(getNextStep()); }

  bool process_next_step(const base::DPBuffer& dpbuffer) {
    return get_next_step()->process(dpbuffer);
  }

  int get_count() { return m_count; }
  void set_input(steps::InputStep* input) { m_input = input; };
  void set_parset(const common::ParameterSet& parset) { m_parset = parset; };
  void set_name(const string& name) { m_name = name; };

  bool m_fetch_uvw = false;
  bool m_fetch_weights = false;

 protected:
  int m_count = 0;
  steps::InputStep* m_input;
  common::ParameterSet m_parset;
  string m_name;
  const base::DPBuffer* m_dpbuffer_in;
  common::NSTimer m_timer;
};

class PyStepImpl final : public StepWrapper {
 public:
  using StepWrapper::StepWrapper;

  void show(std::ostream& os) const override;

  void updateInfo(const base::DPInfo&) override;

  bool process(const base::DPBuffer&) override;

  // Finish the processing of this step and subsequent steps.
  void finish() override;

  void hold();
  void release();

 private:
  // See the comment above near the forward declaration of pybind11::object.
  std::unique_ptr<pybind11::object> m_py_object;
};

}  // namespace pythondp3
}  // namespace dp3

#endif
