# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests are checking the python bindings of the DP3::makeMainSteps function.

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tMakeMainSteps.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/integration/CMakeLists.txt
"""

import gc
import sys
import time

import numpy as np
import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")
import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

from utils import run_in_tmp_path, untar

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass

MSIN = "tNDPPP-generic.MS"


@pytest.fixture
def parset():
    """Provide a parameter set for the make_main_steps tests."""
    parset = dp3.parameterset.ParameterSet()
    parset.add("msin", MSIN)
    parset.add("msout", "averaged_tmp.MS")
    parset.add("msout.overwrite", "true")
    parset.add("steps", "[average]")
    parset.add("average.timestep", "2")
    parset.add("average.freqstep", "8")
    return parset


def test_make_main_steps(parset, run_in_tmp_path):
    """
    This test creates an input, averaging and output step by calling the
    make_main_steps() factory function.
    """

    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")

    first_step = dp3.make_main_steps(parset)

    first_step.set_info(dp3.DPInfo())

    step = first_step
    last_step = step
    while step is not None:
        step.show()
        last_step = step
        step = step.get_next_step()
    assert last_step.info_out.n_times == 3

    start_time = time.time()
    while first_step.process(dp3.DPBuffer()):
        print(".")
    first_step.finish()

    end_time = time.time()
    elapsed_time = end_time - start_time

    # Assert that the "show_timings" function returns the expected results.
    msreader = first_step
    averager = first_step.get_next_step()
    mswriter = averager.get_next_step()

    assert "s) MsReader" in msreader.show_timings(elapsed_time)
    assert "s) Averager" in averager.show_timings(elapsed_time)
    assert "s) MSWriter" in mswriter.show_timings(elapsed_time)


def test_make_main_steps_delete_first(parset, run_in_tmp_path):
    """
    Check the following scenario:
    - make_main_steps() creates multiple steps.
    - get_next_step() is called to get the second step.
    - The first step is deleted and garbage collected.
    In this case, the second step should remain valid.
    Running process() and finish() on the second step should remain working.
    """

    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")

    first_step = dp3.make_main_steps(parset)
    first_step.set_info(dp3.DPInfo())

    second_step = first_step.get_next_step()
    assert second_step.get_next_step()

    del first_step
    gc.collect()

    n_baselines = len(second_step.info_in.first_antenna_indices)
    n_channels = second_step.info_in.n_channels
    n_correlations = second_step.info_in.n_correlations

    for i in range(5):
        data = np.zeros(
            (n_baselines, n_channels, n_correlations), dtype=np.csingle
        )
        weights = np.ones(
            (n_baselines, n_channels, n_correlations), dtype=np.float32
        )
        flags = np.zeros(
            (n_baselines, n_channels, n_correlations), dtype=np.bool
        )
        uvw = np.zeros((n_baselines, 3), dtype=np.double)

        buffer = dp3.DPBuffer(42 + i, 1)
        buffer.set_data(data)
        buffer.set_weights(weights)
        buffer.set_flags(flags)
        buffer.set_uvw(uvw)

        second_step.process(buffer)

    second_step.finish()
