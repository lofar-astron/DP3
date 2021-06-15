# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import shutil
import uuid
from subprocess import check_call

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, untar_ms

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tMultiApplyCal.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
PARMDB_TGZ = "tApplyCal2.parmdb.tgz"  # Note: This archive contains tApplyCal.parmdb.
PARMDB = "tApplyCal.parmdb"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar_ms(f"{tcf.SRCDIR}/{PARMDB_TGZ}")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


@pytest.mark.parametrize(
    "applycal_cmd",
    [
        [
            f"applycal.gain.parmdb={PARMDB}",
            f"applycal.csp.parmdb={PARMDB}",
        ],
        [f"applycal.parmdb={PARMDB}"],
    ],
)
def test_applycal(applycal_cmd):
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "checkparset=1",
            "msout=.",
            "msout.datacolumn=DATA3",
            "steps=[applycal]",
            'applycal.steps="[gain,csp]"',
            "applycal.gain.correction=gain",
            "applycal.csp.correction=commonscalarphase",
            "showcounts=false",
        ]
        + applycal_cmd
    )
    taql_command = f"select from {MSIN} where not(all(DATA~=9*DATA3))"
    assert_taql(taql_command)
