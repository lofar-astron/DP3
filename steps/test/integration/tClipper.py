# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
from subprocess import check_call

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, get_taql_result, run_in_tmp_path, untar

"""
Tests for clipper (coarsely predict and clip bright sources).

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tClipper.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSOUT = "clipper.MS"
TEST_MAX_AMPLITUDE = 1.5


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


@pytest.fixture()
def create_model_data():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[]",
        ]
    )
    with open("clipper.skymodel", "w") as f:
        f.write(
            "FORMAT = Name, Type, Patch, Ra, Dec, I, Q, U, V, ReferenceFrequency='74e6', SpectralIndex='[]', MajorAxis, MinorAxis, Orientation\n"
        )
        f.write(", , CasAPoint, 00:00:00, +00.00.00\n")
        f.write(
            "CasAPoint_1, POINT, CasAPoint, 23:23:20.0, +58.48.26.0, 1e4, 0.0, 0.0, 0.0, 74.0e6, [-0.6]"
        )


def test_clipper(create_model_data):
    taql_command = f"UPDATE {MSOUT} set FLAG=false"
    get_taql_result(taql_command)
    taql_command = f"select from {MSOUT} where any(FLAG)"
    assert_taql(taql_command, 0)

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSOUT}",
            "msout=.",
            "steps=[clipper]",
            "clipper.sourcedb=clipper.skymodel",
            "clipper.usebeammodel=true",
            f"clipper.amplmax={TEST_MAX_AMPLITUDE}",
            "clipper.timestep=1",
            "clipper.freqstep=1",
        ]
    )
    assert_taql(taql_command, 75)
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSOUT}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            "predict.sourcedb=clipper.skymodel",
            "predict.usebeammodel=true",
        ]
    )

    # Check that high amplitudes are flagged.
    taql_command = f"select from {MSOUT} where any(not FLAG && (MODEL_DATA > {TEST_MAX_AMPLITUDE}))"
    assert_taql(taql_command, 0)


def test_freqstep(create_model_data):
    taql_command = f"UPDATE {MSOUT} set FLAG=false"
    get_taql_result(taql_command)
    taql_command = f"select from {MSOUT} where any(FLAG)"
    assert_taql(taql_command, 0)

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSOUT}",
            "msout=.",
            "steps=[clipper]",
            "clipper.sourcedb=clipper.skymodel",
            "clipper.usebeammodel=true",
            f"clipper.amplmax={TEST_MAX_AMPLITUDE}",
            "clipper.timestep=1",
            "clipper.freqstep=2",
        ]
    )
    assert_taql(taql_command, 74)
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSOUT}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            "predict.sourcedb=clipper.skymodel",
            "predict.usebeammodel=true",
        ]
    )

    # Check that most high amplitudes are flagged. Not all amplitudes are flagged,
    # as the visibilities were predicted only for half the frequency channels.
    taql_command = f"select from {MSOUT} where any(not FLAG && (MODEL_DATA > {TEST_MAX_AMPLITUDE}))"
    assert_taql(taql_command, 2)
