# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from . import parameterset, pydp3

"""
    DP3 Python interface

    This module provides the Python interface to DP3 functionality.

    The purpose of this module is two fold
    1) Make DP3 functionality available from Python
    2) Allow to implement DP3 steps in Python, to be used
    either directly from Python, or from the C++ code

    This module is build on top of the lower level Python bindings
    in dp3.pythondp3.
    In general it is prefered to use the higher level interface since
    it provides more checks on the inputs and better error reporting than
    the lower level bindings.
"""


################ RATIONALE ######################################

# C++ steps can segfault when provided erroneous input.
# There are various ways in Python to get a step object that has an
# underlying C++ implementation that needs to be protected from invalid usage:
#
# - instantiating a Python step. Since all Python steps need to be derived from
# the base class pydp3.Step, even the Python steps have an underlying instantiation of
# the C++ base class.
# - calling get_next_step() on an existing Python or C++ instance of Step.
# - calling make_single_step(...) or make_main_steps(...).

# In the latter case the Python side only gets control over the object after
# it has been instantiated on the C++ side. The Python wrapping code can thus
# not be implemented through inheritance. It needs to be a wrapper that
# contains the wrapped object as member.

# On the other hand Python steps can be instantiated directly from the Python side.
# To ensure the step is always wrapped, the wrapper code needs to be part of
# the inheritance chain.

# To satisfy both requirements there need to be two different classes
# implementing protections:
# 1) A StepWrapper class that contains an underlying C++ step of class pydp3.Step.
# The methods of StepWrapper check the input parameters, forward calls to C++,
# and optionally check/convert the result.
# 2) A Step class which is derived from pydp3.Step. It only wraps the methods that
# are implemented in the C++ base class. These methods check the inputs, forward
# the call to the C++ base class, and optionally check/convert the result.

############################################################################


# For now only the Step class has a wrapper
# Other functions and classes are imported directly from the pybind11 Python bindings.
# TODO As part "AST-1323 [DP3] Improve sanity checks in Python interface"
# these imports need to be replaced by wrappers that check inputs and provide other
# conveniences for Python users.

from .pydp3 import (
    DPBuffer,
    DPInfo,
    Fields,
    MsType,
    get_chain_required_fields,
    get_n_threads,
    set_n_threads,
)


class Step(pydp3.Step):
    """
    Base class for Python implemented Steps
    """

    def __init__(self, *args):
        super().__init__(*args)

        # StepWrapper has a _step member that contains the wrapped step.
        # For convenience, add a _step member in Step, too.
        self._step = self

    def get_next_step(self):
        """
        Get the next Step

        If needed the next step will be wrapped
        """
        next_step = super().get_next_step()
        return StepWrapper(next_step)

    def __eq__(self, other):
        """
        A Step and a StepWrapper compare equal if the underlying object is the same
        """

        if isinstance(other, StepWrapper):
            return self is other._step
        else:
            return self is other


class StepWrapper:
    """
    Wrapper containing an underlying pydp3.Step
    """

    def __new__(cls, step):
        """
        Instantiation of StepWrapper objects handles two special cases
        1) The underlying object is a None value. Instead of a StepWrapper instance, a None value is returned
        2) A Python step does not need wrapping, the step itself is returned
        """

        # A wrapped None value is still a None value,
        # not a StepWrapper around a None value
        if step is None:
            return None

        # Wrapping is not needed when the step derives from dp3.Step
        if isinstance(step, Step):
            return step
        else:
            obj = super().__new__(cls)
            obj._step = step
            return obj

    def __eq__(self, other):
        """
        A Step and a StepWrapper compare equal if the underlying object is the same
        """
        return self._step is other._step

    def show(self):
        self._step.show()

    def show_timings(self, elapsed_time):
        """
        Show the processing time of current step. Also provides the
        percentage of time the current step took compared to the total
        elapsed time given as input (in seconds).
        """

        return self._step.show_timings(elapsed_time)

    def __str__(self):
        return self._step.__str__()

    def process(self, buffer):
        info = self.info
        n_baselines = len(info.first_antenna_indices)
        n_channels = info.n_channels
        n_correlations = info.n_correlations
        if (n_baselines == 0) or (n_channels == 0) or (n_correlations == 0):
            msg = (
                "The step's info contains invalid data\n"
                + f"(n_baselines, n_channels, n_correlations) = ({n_baselines}, {n_channels}, {n_correlations})\n"
                + "set_info() should be called with valid info, before calling process."
            )
            raise RuntimeError(msg)

        info_in = self.info_in
        n_baselines_in = len(info_in.first_antenna_indices)
        n_channels_in = info_in.n_channels
        n_correlations_in = info_in.n_correlations

        required_fields = get_chain_required_fields(self._step)

        if required_fields.data:
            data = np.array(buffer.get_data(), copy=False)
            if not data.shape == (
                n_baselines_in,
                n_channels_in,
                n_correlations_in,
            ):
                msg = (
                    "The dimensions in this step's info_in "
                    + f"(n_baselines, n_channels, n_correlations) = ({n_baselines_in}, {n_channels_in}, {n_correlations_in})\n"
                    + f"do not match the dimensions of the data in the buffer {data.shape}"
                )
                raise RuntimeError(msg)

        if required_fields.weights:
            weights = np.array(buffer.get_weights(), copy=False)
            if not weights.shape == (
                n_baselines_in,
                n_channels_in,
                n_correlations_in,
            ):
                msg = (
                    "The dimensions in this step's info_in "
                    + f"(n_baselines, n_channels, n_correlations) = ({n_baselines_in}, {n_channels_in}, {n_correlations_in})\n"
                    + f"do not match the dimensions of the weights in the buffer {weights.shape}"
                )
                raise RuntimeError(msg)

        if required_fields.flags:
            flags = np.array(buffer.get_flags(), copy=False)
            if not flags.shape == (
                n_baselines_in,
                n_channels_in,
                n_correlations_in,
            ):
                msg = (
                    "The dimensions in this step's info_in "
                    + f"(n_baselines, n_channels, n_correlations) = ({n_baselines_in}, {n_channels_in}, {n_correlations_in})\n"
                    + f"do not match the dimensions of the flags in the buffer {flags.shape}"
                )
                raise RuntimeError(msg)

        if required_fields.uvw:
            uvw = np.array(buffer.get_uvw(), copy=False)
            if not uvw.shape == (n_baselines_in, 3):
                msg = (
                    "The dimensions in this step's info_in "
                    + f"(n_baselines, 3) = ({n_baselines_in}, 3)\n"
                    + f"do not match the dimensions of the uvw in the buffer {uvw.shape}"
                )
                raise RuntimeError(msg)

        return self._step.process(buffer)

    def set_info(self, info):
        return self._step.set_info(info)

    @property
    def info_in(self):
        "Get a copy of the info object containing metadata of the input"
        return self._step.info_in

    @property
    def info(self):
        "Get a copy of the info object containing metadata"
        return self._step.info

    def finish(self):
        return self._step.finish()

    def get_next_step(self):
        """Get a reference to the next step"""
        next_step = self._step.get_next_step()
        return StepWrapper(next_step)

    def set_next_step(self, next_step):
        "Set the step that follows the current step"

        if isinstance(next_step, StepWrapper):
            self._step.set_next_step(next_step._step)
        elif isinstance(next_step, pydp3.Step):
            self._step.set_next_step(next_step)
        else:
            raise RuntimeError(
                "next_step is not an instance of dp3.StepWrapper, dp3.Step nor dp3.pydp3.Step"
            )

    def get_required_fields(self):
        "Get the fields required by the current step"
        return self._step.get_required_fields()

    def get_provided_fields(self):
        "Get the fields provided by the current step"
        return self._step.get_provided_fields()


def make_step(type, parset, prefix, input_type):
    """
    Create a single step of a certain type
    Settings are extracted from the parset
    The prefix should include the dot
    """
    step = pydp3.make_step(type, parset, prefix, input_type)
    if step is not None:
        null_step = pydp3.make_step(
            "null", parameterset.ParameterSet(), "", MsType.regular
        )
        step.set_next_step(null_step)
    return StepWrapper(step)


def make_main_steps(parset):
    """
    Create a chain of steps from a parset
    The first step of the chain is returned
    The rest of the steps can be retrieved by calling
    the get_next_step() iteratively
    """
    step = pydp3.make_main_steps(parset)
    return StepWrapper(step)
