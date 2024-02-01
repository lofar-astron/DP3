# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import shutil
import os
import uuid
from subprocess import check_call, check_output

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, untar_ms, get_taql_result

"""
Tests for clipper (coarsely predict and clip bright sources).

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tClipper.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSOUT = "clipper.MS"
CWD = os.getcwd()
TEST_AMPLMAX = 1.5


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


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


def test_clipper(create_model_data):
    with open("clipper.skymodel", "w") as f:
        f.write(
            "FORMAT = Name, Type, Patch, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\n"
        )
        f.write(", , dummy, 01:37:41.3, +15.09.35.132\n")
        f.write(
            "src_0, POINT, dummy, 01:37:41.3, +15.09.35.132, 1000, , , , ,"
        )

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
            f"clipper.amplmax={TEST_AMPLMAX}",
            "clipper.timestep=1",
        ]
    )
    assert_taql(taql_command, 102)
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
    taql_command = f"select from {MSOUT} where any(not FLAG && (MODEL_DATA > {TEST_AMPLMAX}))"
    assert_taql(taql_command, 0)
