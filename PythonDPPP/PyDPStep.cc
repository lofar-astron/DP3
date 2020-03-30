#include "PyDPStep.h"
#include "PyDPStepImpl.h"

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h> // everything needed for embedding

namespace DP3 {
namespace DPPP {

void PyDPStepImpl ::show (std::ostream& os) const
{
    // TODO redirect sys.stdout to os

    pybind11::gil_scoped_acquire gil;  // Acquire the GIL while in this scope.
    // Try to look up the overloaded method on the Python side.
    pybind11::function overload = pybind11::get_overload(this, "show");
    if (overload) overload();  // Call the Python function.
}

bool PyDPStepImpl::process (const DPBuffer& bufin)
{
    m_count++;

    auto dpbuffer = std::shared_ptr<DPBuffer>(new DPBuffer(bufin));

    // fetch optional data
    if (m_fetch_uvw)
    {
        m_input->fetchUVW(bufin, *dpbuffer, m_timer);
    }

    if (m_fetch_weights)
    {
        m_input->fetchWeights(bufin, *dpbuffer, m_timer);
    }

    PYBIND11_OVERLOAD_PURE(
        bool,           /* Return type */
        DPStepWrapper,  /* Parent class */
        process,        /* Name of function in C++ (must match Python name) */
        dpbuffer        /* Argument(s) */
    );
}

void PyDPStepImpl::finish ()
{
    PYBIND11_OVERLOAD_PURE(
        void,           /* Return type */
        DPStepWrapper,  /* Parent class */
        finish          /* Name of function in C++ (must match Python name) */
    );
}

void PyDPStepImpl::updateInfo (const DPInfo& dpinfo)
{
    PYBIND11_OVERLOAD_NAME(
        void,          /* Return type */
        DPStepWrapper, /* Parent class */
        "update_info",
        updateInfo,    /* Name of function in C++ (must match Python name) */
        dpinfo         /* Argument(s) */
    );
}



void PyDPStepImpl::hold()
{
    m_py_object = py::cast(static_cast<DPStep*>(this));
}

DPStep::ShPtr PyDPStep::create_instance(
    DPInput* input,
    const ParameterSet& parset,
    const string& name)
{
    static py::scoped_interpreter guard{}; // start the interpreter and keep it alive
    static py::module mydpstep_module = py::module::import("dppp.idgcalstep");

    auto pydpstep_instance = mydpstep_module.attr("IDGCalStep")(parset, name);
    auto pydpstep_instance_ptr = pydpstep_instance.cast<std::shared_ptr<PyDPStepImpl>>();
    pydpstep_instance_ptr->hold();
    pydpstep_instance_ptr->set_input(input);
    pydpstep_instance_ptr->set_parset(parset);
    pydpstep_instance_ptr->set_name(name);


    return pydpstep_instance_ptr;
}

} // namespace DPPP
} // namespace DP3
