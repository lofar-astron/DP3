# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests are checking the python bindings of the DP3::makeMainSteps function. 

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tMakeMainSteps.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/integration/CMakeLists.txt
"""

import pytest
import os
import sys

# Append current directory to system path in order to import testconfig
sys.path.append(".")
import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

from utils import assert_taql, untar_ms

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass

MSIN = "tNDPPP-generic.MS"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env(tmp_path):
    os.chdir(tmp_path)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)


def test_make_main_steps(tmp_path):
    """
    This test creates an input, averaging and output step by calling the
    make_main_steps() factory function.
    """

    msout = str(tmp_path / "averaged_tmp.MS")

    parset = dp3.parameterset.ParameterSet()

    parset.add("msin", MSIN)
    parset.add("msout", msout)
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

    while first_step.process(dp3.DPBuffer()):
        print(".")
    first_step.finish()
