// PyDPStep.cc: "template" for the python Step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "PyStep.h"
#include "PyStepImpl.h"
#include "../base/DP3.h"

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>  // everything needed for embedding

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using dp3::steps::InputStep;
using dp3::steps::Step;

namespace dp3 {
namespace pythondp3 {

void PyStepImpl::show(std::ostream& os) const {
  pybind11::gil_scoped_acquire gil;  // Acquire the GIL while in this scope.

  // redirect sys.stdout to os
  pybind11::object sysm = pybind11::module::import("sys");
  pybind11::object stdout = sysm.attr("stdout");
  sysm.attr("stdout") = ostream_wrapper(os);

  // Try to look up the overloaded method on the Python side.
  pybind11::function overload = pybind11::get_overload(this, "show");
  if (overload) overload();  // Call the Python function.

  // restore sys.stdout
  sysm.attr("stdout") = stdout;
}

void PyStepImpl::finish() {
  pybind11::function overload = pybind11::get_overload(this, "finish");
  if (overload) overload();  // Call the Python function.

  getNextStep()->finish();
}

bool PyStepImpl::process(const DPBuffer& bufin) {
  m_count++;

  // Make a deep copy of the buffer to make the data
  // persistent across multiple process calls
  // This is not always necessary, but for python Steps
  // convenience is more important than performance
  auto dpbuffer = std::shared_ptr<DPBuffer>(new DPBuffer());
  dpbuffer->copy(bufin);

  // fetch optional data
  if (m_fetch_uvw) {
    m_input->fetchUVW(bufin, *dpbuffer, m_timer);
  }

  if (m_fetch_weights) {
    m_input->fetchWeights(bufin, *dpbuffer, m_timer);
  }

  PYBIND11_OVERLOAD_PURE(
      bool,        /* Return type */
      StepWrapper, /* Parent class */
      process,     /* Name of function in C++ (must match Python name) */
      dpbuffer     /* Argument(s) */
  );
}

void PyStepImpl::updateInfo(const DPInfo& dpinfo) {
  PYBIND11_OVERLOAD_NAME(void,          /* Return type */
                         StepWrapper,   /* Parent class */
                         "update_info", /* Name of function in Python  */
                         updateInfo,    /* Name of function in C++  */
                         dpinfo         /* Argument(s) */
  );
}

void PyStepImpl::hold() {
  m_py_object.reset(
      new pybind11::object(pybind11::cast(static_cast<Step*>(this))));
}

void PyStepImpl::release() { m_py_object.reset(); }

Step::ShPtr PyStep::create_instance(InputStep* input,
                                    const common::ParameterSet& parset,
                                    const string& prefix) {
  std::string module_name = parset.getString(prefix + "python.module");
  std::string class_name = parset.getString(prefix + "python.class");

  static pybind11::scoped_interpreter
      guard{};  // start the interpreter and keep it alive
  static pybind11::module mydpstep_module =
      pybind11::module::import(module_name.c_str());

  auto pydpstep_instance =
      mydpstep_module.attr(class_name.c_str())(parset, prefix);

  // Create a shared_ptr from the raw pointer
  // Note that the python side uses another smart pointer as holder type
  // The deleter of the shared_ptr here will only release the python object
  // If there are no other references at the python side this will trigger
  // a delete of the C++ object
  auto pydpstep_instance_ptr = std::shared_ptr<PyStepImpl>(
      pydpstep_instance.cast<PyStepImpl*>(), std::mem_fn(&PyStepImpl::release));

  pydpstep_instance_ptr->hold();
  pydpstep_instance_ptr->set_input(input);

  return pydpstep_instance_ptr;
}

}  // namespace pythondp3
}  // namespace dp3
