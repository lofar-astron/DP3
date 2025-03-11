# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests are checking that the python bindings step factory function.

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tMakeStep.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/integration/CMakeLists.txt
"""

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)
sys.path.insert(0, tcf.PYTHONMOCKDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass


def test_make_averaging_step():
    """
    This test creates an averaging step by calling the make_step() factory function.
    The show() method of the step is called through a conversion to a str.
    The output is compared to the expected output.
    """

    parset = dp3.parameterset.ParameterSet()

    parset.add("average.timestep", "10")
    parset.add("average.freqstep", "8")
    parset.add("average.minpoints", "2")
    parset.add("average.minperc", "1")

    step = dp3.make_step("averager", parset, "average.", dp3.MsType.regular)

    step_description = str(step)
    expected_step_description = (
        "Averager average.\n"
        "  freqstep:       8  timestep:       10\n"
        "  minpoints:      2\n"
        "  minperc:        1\n"
    )

    assert step_description == expected_step_description


def test_make_mock_pystep():
    """
    This test creates a python step by calling the make_step() factory function.

    Note that we are calling a C++ function from python which in turn calls
    python functions to create the python step.
    This shows that pybind11 wrapping and embedding can be mixed and
    objects can be passed back and forth across the Python/C++ boundary.

    The show() method of the step is called through a conversion to a str.
    The output is compared to the expected output.
    """

    parset = dp3.parameterset.ParameterSet()

    parset.add("mock.python.module", "mockpystep")
    parset.add("mock.python.class", "MockPyStep")
    parset.add("mock.datafactor", "2")
    parset.add("mock.weightsfactor", "0.5")

    step = dp3.make_step("python", parset, "mock.", dp3.MsType.regular)

    step_description = str(step)
    expected_step_description = (
        "\nMockPyStep\n" "  data factor:    2.0\n" "  weights factor: 0.5\n"
    )
    assert step_description == expected_step_description


def test_make_invalid_pystep():
    """
    This test attempts to create a step of non-existent type.
    The returned step should be None
    """

    parset = dp3.parameterset.ParameterSet()
    step = dp3.make_step("invalid", parset, "", dp3.MsType.regular)
    assert step is None
