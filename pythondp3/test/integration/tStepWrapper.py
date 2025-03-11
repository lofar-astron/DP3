# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests are checking that StepWrapper from the python bindings
indeed catches a few easily made user errors

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tStepWrapper.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/integration/CMakeLists.txt
"""

import sys

import numpy as np
import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)
sys.path.insert(0, tcf.PYTHONMOCKDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass


def test_wrapped_average_step():
    """
    This test creates an averaging step by calling the make_step() factory function.
    This C++ step is protected by a StepWrapper from invalid input.
    """

    parset = dp3.parameterset.ParameterSet()

    parset.add("average.timestep", "10")
    parset.add("average.freqstep", "4")
    parset.add("average.minpoints", "2")
    parset.add("average.minperc", "1")

    step = dp3.make_step("averager", parset, "average.", dp3.MsType.regular)
    buffer = dp3.DPBuffer()

    # Try to run process, should fail on missing info
    with pytest.raises(RuntimeError) as e_info:
        step.process(buffer)
    assert "set_info() should be called" in str(e_info.value)

    # Fill dpinfo object

    n_correlations = 4
    dpinfo = dp3.DPInfo(n_correlations)
    antenna1 = [0, 0, 1]
    antenna2 = [1, 2, 2]
    n_baselines = len(antenna1)
    antenna_names = ["ant1", "ant2", "ant3"]
    antenna_positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    antenna_diameters = [10.0, 10.0, 10.0]
    dpinfo.set_antennas(
        antenna_names, antenna_diameters, antenna_positions, antenna1, antenna2
    )

    n_channels = 8
    dpinfo.set_channels(np.arange(n_channels), np.ones(n_channels))
    first_time = 0
    last_time = 10
    time_interval = 1.0
    dpinfo.set_times(first_time, last_time, time_interval)

    step.set_info(dpinfo)

    # The average step needs all (data, weights, flags, uvw) fields
    # First verify that it indeed needs thems all
    required_fields = step.get_required_fields()
    assert required_fields.data
    assert required_fields.weights
    assert required_fields.flags
    assert required_fields.uvw

    # Now start with an empty buffer and try to call the process method.
    # Fill one extra field at a time and try a again
    # Each call should fail with a different error, until the buffer is complete

    # Try process again, should fail on missing data
    with pytest.raises(RuntimeError) as e_info:
        step.process(buffer)
    assert "do not match the dimensions of the data" in str(e_info.value)

    # Fill data in buffer
    buffer.set_data(
        np.zeros((n_baselines, n_channels, n_correlations), np.complex64)
    )
    # Try process again, should fail on missing weights
    with pytest.raises(RuntimeError) as e_info:
        step.process(buffer)
    assert "do not match the dimensions of the weights" in str(e_info.value)

    # Fill weights
    buffer.set_weights(
        np.zeros((n_baselines, n_channels, n_correlations), np.float32)
    )
    # Try process again, should fail on missing flags
    with pytest.raises(RuntimeError) as e_info:
        step.process(buffer)
    assert "do not match the dimensions of the flags" in str(e_info.value)

    # Fill flags
    buffer.set_flags(np.zeros((n_baselines, n_channels, n_correlations), bool))
    # Try process again, should fail on missing uvw
    with pytest.raises(RuntimeError) as e_info:
        step.process(buffer)
    assert "do not match the dimensions of the uvw" in str(e_info.value)

    # Fill uvw
    buffer.set_uvw(np.zeros((n_baselines, 3), np.float64))
    # Try process again, should succeed now
    step.process(buffer)


def test_wrapped_uvw_flagger_step():
    """
    This test creates an uvw flagger step by calling the make_step() factory function.
    This C++ step is protected by a StepWrapper from invalid input.
    """

    parset = dp3.parameterset.ParameterSet()

    parset.add("uvwflagger.uvlambdamax", "1000")

    step = dp3.make_step(
        "uvwflagger", parset, "uvwflagger.", dp3.MsType.regular
    )
    buffer = dp3.DPBuffer()

    # Fill dpinfo object

    n_correlations = 4
    dpinfo = dp3.DPInfo(n_correlations)
    antenna1 = [0, 0, 1]
    antenna2 = [1, 2, 2]
    n_baselines = len(antenna1)
    antenna_names = ["ant1", "ant2", "ant3"]
    antenna_positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    antenna_diameters = [10.0, 10.0, 10.0]
    dpinfo.set_antennas(
        antenna_names, antenna_diameters, antenna_positions, antenna1, antenna2
    )

    n_channels = 8
    dpinfo.set_channels(np.arange(n_channels), np.ones(n_channels))
    first_time = 0
    last_time = 10
    time_interval = 1.0
    dpinfo.set_times(first_time, last_time, time_interval)

    step.set_info(dpinfo)

    # The uvw flagger step only needs the flags and uvw fields
    # First verify that it indeed needs only those
    required_fields = step.get_required_fields()
    assert not required_fields.data
    assert not required_fields.weights
    assert required_fields.flags
    assert required_fields.uvw

    # Now start with an empty buffer and try to call the process method.
    # Fill one extra required field at a time and try a again
    # Each call should fail with a different error, until the buffer
    # contains all required is complete

    # Try process, should fail on missing flags
    with pytest.raises(RuntimeError) as e_info:
        step.process(buffer)
    assert "do not match the dimensions of the flags" in str(e_info.value)

    # Fill flags
    buffer.set_flags(np.zeros((n_baselines, n_channels, n_correlations), bool))
    # Try process again, should fail on missing uvw
    with pytest.raises(RuntimeError) as e_info:
        step.process(buffer)
    assert "do not match the dimensions of the uvw" in str(e_info.value)

    # Fill uvw
    buffer.set_uvw(np.zeros((n_baselines, 3), np.float64))
    # Try process again, should succeed now
    step.process(buffer)
