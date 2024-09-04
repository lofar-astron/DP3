# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests are checking the python bindings of the DP3::makeMainSteps function.

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tMakeMainSteps.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/integration/CMakeLists.txt
"""

import os
import sys
import time

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")
import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

from utils import assert_taql, run_in_tmp_path, untar

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass

MSIN = "tNDPPP-generic.MS"


def test_make_main_steps(run_in_tmp_path):
    """
    This test creates an input, averaging and output step by calling the
    make_main_steps() factory function.
    """

    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")

    parset = dp3.parameterset.ParameterSet()

    parset.add("msin", MSIN)
    parset.add("msout", "averaged_tmp.MS")
    parset.add("msout.overwrite", "true")
    parset.add("steps", "[average]")
    parset.add("average.timestep", "2")
    parset.add("average.freqstep", "8")

    first_step = dp3.make_main_steps(parset)

    last_info = first_step.set_info(dp3.DPInfo())
    assert last_info.n_times == 3

    step = first_step
    while step is not None:
        step.show()
        step = step.get_next_step()

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
