# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import shutil
import uuid
from subprocess import check_call, check_output

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, untar_ms

"""
Tests for applying the beam model.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tApplyBeam.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSAPPLYBEAM = "tApplyBeam.tab"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar_ms(f"{tcf.SRCDIR}/{MSAPPLYBEAM}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


@pytest.mark.parametrize("usechannelfreq", [False, True])
def test_with_invert_and_chanfreq(usechannelfreq):
    msout = "outinv.ms"

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={msout}",
            "steps=[applybeam]",
            f"applybeam.usechannelfreq={'true' if usechannelfreq else 'false'}",
            "applybeam.invert=true",
        ]
    )
    ms_column = "DATA_ucf" if usechannelfreq else "DATA_noucf"
    taql_command = f"select from {msout} t1, {MSAPPLYBEAM} t2 where not all(near(t1.DATA,t2.{ms_column},8e-5) || (isnan(t1.DATA) && isnan(t2.{ms_column})))"
    assert_taql(taql_command)

    # Test with invert=false on the resulting outinv.ms
    check_call(
        [
            tcf.DP3EXE,
            f"msin={msout}",
            f"msout=out.ms",
            "steps=[applybeam]",
            f"applybeam.usechannelfreq={'true' if usechannelfreq else 'false'}",
            "applybeam.invert=false",
        ]
    )
    taql_command = f"select from out.ms t1, {MSIN} t2 where not all(near(t1.DATA,t2.DATA,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA)))"
    assert_taql(taql_command)


@pytest.mark.parametrize("beammode", ["ARRAY_FACTOR", "ELEMENT"])
def test_beammodes(beammode):
    msout = "outinv.ms"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={msout}",
            "steps=[applybeam]",
            "applybeam.usechannelfreq=true",
            f"applybeam.beammode={beammode}",
            "applybeam.invert=true",
        ]
    )
    taql_command = f"select from {msout} t1, {MSAPPLYBEAM} t2 where not all(near(t1.DATA,t2.DATA_{beammode},8e-5) || (isnan(t1.DATA) && isnan(t2.DATA_{beammode})))"
    assert_taql(taql_command)


def test_with_updateweights():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=.",
            f"msout.weightcolumn=NEW_WEIGHT_SPECTRUM",
            "steps=[applybeam]",
            f"applybeam.updateweights=true",
        ]
    )
    taql_command = (
        f"select from {MSIN} where all(near(WEIGHT_SPECTRUM, NEW_WEIGHT_SPECTRUM))"
    )
    assert_taql(taql_command)
