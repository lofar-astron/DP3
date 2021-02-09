// PyDPStep.cc: "template" for the python DPStep
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "PyDPStep.h"
#include "PyDPStepImpl.h"
#include "../DPPP/DPRun.h"

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>  // everything needed for embedding

namespace DP3 {
namespace DPPP {

void PyDPStepImpl::show(std::ostream& os) const {
  pybind11::gil_scoped_acquire gil;  // Acquire the GIL while in this scope.

  // redirect sys.stdout to os
  pybind11::object sysm = pybind11::module::import("sys");
  pybind11::object stdout = sysm.attr("stdout");
  sysm.attr("stdout") = ostream(os);

  // Try to look up the overloaded method on the Python side.
  pybind11::function overload = pybind11::get_overload(this, "show");
  if (overload) overload();  // Call the Python function.

  // restore sys.stdout
  sysm.attr("stdout") = stdout;
}

bool PyDPStepImpl::process(const DPBuffer& bufin) {
  m_count++;

  // Make a deep copy of the buffer to make the data
  // persistent across multiple process calls
  // This is not always necessary, but for python DPSteps
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
      bool,          /* Return type */
      DPStepWrapper, /* Parent class */
      process,       /* Name of function in C++ (must match Python name) */
      dpbuffer       /* Argument(s) */
  );
}

void PyDPStepImpl::finish() {
  PYBIND11_OVERLOAD_PURE(
      void,          /* Return type */
      DPStepWrapper, /* Parent class */
      finish         /* Name of function in C++ (must match Python name) */
  );
}

void PyDPStepImpl::updateInfo(const DPInfo& dpinfo) {
  PYBIND11_OVERLOAD_NAME(
      void,          /* Return type */
      DPStepWrapper, /* Parent class */
      "update_info",
      updateInfo, /* Name of function in C++ (must match Python name) */
      dpinfo      /* Argument(s) */
  );
}

void PyDPStepImpl::hold() {
  m_py_object.reset(
      new pybind11::object(pybind11::cast(static_cast<DPStep*>(this))));
}

DPStep::ShPtr PyDPStep::create_instance(DPInput* input,
                                        const ParameterSet& parset,
                                        const string& prefix) {
  std::string module_name = parset.getString(prefix + "python.module");
  std::string class_name = parset.getString(prefix + "python.class");

  static pybind11::scoped_interpreter
      guard{};  // start the interpreter and keep it alive
  static pybind11::module mydpstep_module =
      pybind11::module::import(module_name.c_str());

  auto pydpstep_instance =
      mydpstep_module.attr(class_name.c_str())(parset, prefix);
  auto pydpstep_instance_ptr =
      pydpstep_instance.cast<std::shared_ptr<PyDPStepImpl>>();
  pydpstep_instance_ptr->hold();
  pydpstep_instance_ptr->set_input(input);

  return pydpstep_instance_ptr;
}

}  // namespace DPPP
}  // namespace DP3
