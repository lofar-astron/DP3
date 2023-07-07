// pydp3.cc: python bindings to create python step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

#include <dp3/base/DP3.h>
#include <dp3/steps/Step.h>

#include "PyDpBuffer.h"
#include "PyStep.h"
#include "utils.h"

using dp3::steps::Step;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

void WrapDpBuffer(py::module &m);  // Defined in PyDpBuffer.cc
void WrapDpInfo(py::module &m);    // Defined in pydpinfo.cc
void WrapFields(py::module &m);    // Defined in pyfields.cc

class PublicStep : public Step {
 public:
  using Step::info;
  using Step::infoOut;
  using Step::showTimings;
  using Step::updateInfo;
};

PYBIND11_MODULE(pydp3, m) {
  m.doc() = "pybind11 example plugin";  // optional module docstring

  WrapDpBuffer(m);
  WrapDpInfo(m);
  WrapFields(m);

  py::class_<ostream_wrapper>(m, "ostream")
      .def("write", &ostream_wrapper::write);

  m.def("make_step", &dp3::base::MakeSingleStep);
  m.def("make_main_steps",
        [](const common::ParameterSet &parset) -> std::shared_ptr<steps::Step> {
          return dp3::base::MakeMainSteps(parset);
        });
  m.def("get_chain_required_fields", &dp3::base::GetChainRequiredFields);

  py::enum_<dp3::steps::Step::MsType>(m, "MsType")
      .value("regular", dp3::steps::Step::MsType::kRegular)
      .value("bda", dp3::steps::Step::MsType::kBda);

  py::class_<dp3::steps::Step, std::shared_ptr<dp3::steps::Step>, PyStep>(
      m, "Step")
      .def(py::init<>())
      .def(
          "show", [](const dp3::steps::Step &step) { step.show(std::cout); },
          "Show step summary. When running embedded in DP3, sys.stdout will be "
          "redirected to DP3's output stream. Overrides of show() should use "
          "print() "
          "to output a summary of the step.")
      .def(
          "show_timings",
          [](const dp3::steps::Step &step, double elapsed_time) {
            std::stringstream ss;
            step.showTimings(ss, elapsed_time);
            return ss.str();
          },
          "Show the processing time of current step. Also provides the "
          "percentage of time the current step took compared to the total "
          "elapsed time given as input (in seconds).")
      .def("__str__",
           [](const dp3::steps::Step &step) {
             std::stringstream ss;
             step.show(ss);
             return ss.str();
           })
      // Step::updateInfo is protected
      // Prepending by an underscore to indicate that this mehod is not
      // supposed to be called directly
      .def("_update_info", &PublicStep::updateInfo, "Handle metadata")
      .def("set_info", &Step::setInfo, py::return_value_policy::reference,
           "Set info object. This will call _update_info() for this step and "
           "all next steps")
      // Step::getInfoIn(), Step::getInfoOut() and Step::getInfo()
      // return const references.
      // Unfortunately there is no way to enforce constness in Python.
      // Also there is no guarantee that the Step will outlive the returned
      // info object. Although not the most efficient, the safest way
      // is returning a copy here.
      .def_property_readonly(
          "info_in", &Step::getInfoIn, py::return_value_policy::copy,
          "Get a copy of the info object containing metadata of the input")
      .def_property_readonly(
          "info_out", &Step::getInfoOut, py::return_value_policy::copy,
          "Get a copy of the info object containing metadata of the output")
      // Since a python step needs to be able to adjust its info_out in its
      // _update_info() override, _info_out returns a modifyable reference.
      .def_property_readonly("_info_out", &PublicStep::infoOut,
                             py::return_value_policy::reference,
                             "Get a modifyable reference to info object "
                             "containing metadata of the output")
      // Legacy version of info_out
      .def_property_readonly("info", &PublicStep::getInfo,
                             py::return_value_policy::reference,
                             "Get a modifyable reference to info object "
                             "containing metadata of the output")
      .def(
          "process",
          [](dp3::steps::Step &step, PyDpBuffer &buffer) {
            return step.process(buffer.take());
          },
          "process buffer")
      .def("finish", &Step::finish,
           "Finish processing (nextstep->finish will be called automatically")
      .def("get_next_step", &Step::getNextStep,
           "Get a reference to the next step")
      .def("set_next_step", &Step::setNextStep,
           "Set the step that follows the current step")
      .def("get_required_fields", &Step::getRequiredFields,
           "Get the fields required by the current step")
      .def("get_provided_fields", &Step::getProvidedFields,
           "Get the fields provided by the current step");
}

}  // namespace pythondp3
}  // namespace dp3
