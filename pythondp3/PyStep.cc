// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "PyStep.h"

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>  // everything needed for embedding

#include "PyDpBuffer.h"
#include "utils.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using dp3::steps::Step;

namespace dp3 {
namespace pythondp3 {

PyStep::PyStep() {}

void PyStep::show(std::ostream& os) const {
  pybind11::gil_scoped_acquire gil;  // Acquire the GIL while in this scope.

  // redirect sys.stdout to os
  pybind11::object sysm = pybind11::module::import("sys");
  pybind11::object stdout = sysm.attr("stdout");
  sysm.attr("stdout") = ostream_wrapper(os);

  // Try to look up the overridden method on the Python side.
  pybind11::function show_override = pybind11::get_override(this, "show");
  if (show_override) show_override();  // Call the Python function.

  // restore sys.stdout
  sysm.attr("stdout") = stdout;
}

void PyStep::finish() {
  // The default behaviour of the Python version of finish() deviates somewhat
  // from the behaviour of the C++ versions. This is to make Python
  // steps a bit more robust.
  // The differences between the C++ and Python versions of finish() are:
  // 1) In C++ finish() is pure virtual.
  //    In Python overriding the finish() method is optional
  // 2) A Python override should not call the finish() mehod of the next step.
  //    If there is a next step its finish() method will be called from here,
  //    the trampoline class.
  // Conditionally calling the finish() method of the next step allows to write
  // Python steps that can be used as end point of a pipeline
  // like for example the QueueOutput step.
  // Ticket AST-1105 aims to be bring the C++ and Python behaviour more in line.
  pybind11::function finish_override = pybind11::get_override(this, "finish");
  if (finish_override) finish_override();  // Call the Python function.
  if (getNextStep()) getNextStep()->finish();
}

bool PyStep::process(std::unique_ptr<DPBuffer> buffer) {
  PYBIND11_OVERRIDE_PURE(
      bool,    /* Return type */
      Step,    /* Parent class */
      process, /* Name of function in C++ (must match Python name) */
      PyDpBuffer(std::move(buffer)) /* Argument(s) */
  );
}

common::Fields PyStep::getRequiredFields() const {
  PYBIND11_OVERRIDE_PURE_NAME(
      common::Fields, Step, "get_required_fields", getRequiredFields,
      /* no arguments*/);
}

common::Fields PyStep::getProvidedFields() const {
  PYBIND11_OVERRIDE_PURE_NAME(
      common::Fields, Step, "get_provided_fields", getProvidedFields,
      /* no arguments*/);
}

void PyStep::updateInfo(const DPInfo& dpinfo) {
  // dpinfo is const
  // In Python constness can not be enforced
  // Passing a copy here to protect dpinfo from modications on the Python side
  PYBIND11_OVERRIDE_NAME(void,           /* Return type */
                         Step,           /* Parent class */
                         "_update_info", /* Name of function in Python  */
                         updateInfo,     /* Name of function in C++  */
                         DPInfo(dpinfo)  /* Argument(s) */
  );
}

std::shared_ptr<PyStep> PyStep::create_instance(
    const common::ParameterSet& parset, const std::string& prefix) {
  std::string module_name = parset.getString(prefix + "python.module");
  std::string class_name = parset.getString(prefix + "python.class");

  try {
    pybind11::initialize_interpreter();
  } catch (std::runtime_error& e) {
    if (strcmp(e.what(), "The interpreter is already running") != 0) throw;
  }

  pybind11::module mydpstep_module =
      pybind11::module::import(module_name.c_str());

  // Create the python step on the heap
  // This object needs to outlive the scope of this function
  // to prevent the python step from being garbage collected.
  // When using a std::unique_ptr<pybind11::object> member inside PyStep
  // instead, GCC will give a warning regarding visibility.
  auto pyobject_step = new pybind11::object(
      mydpstep_module.attr(class_name.c_str())(parset, prefix));

  // shared_ptr with custom deleter
  // The deleter deletes the pybind11:object
  // That will decrease the reference count on the python side.
  // Python's garbage collector will cleanup the actual PyStep object.
  auto pystep = std::shared_ptr<PyStep>(
      pyobject_step->cast<PyStep*>(),
      [pyobject_step](PyStep*) { delete pyobject_step; });

  return pystep;
}

}  // namespace pythondp3
}  // namespace dp3
