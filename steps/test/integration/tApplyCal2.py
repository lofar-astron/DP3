# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call

import pytest

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Tests for applying the calibration.

Script can be invoked in two ways:
- as standalone from the build/steps/tests/integration directory,
  using `pytest source/tApplyCal2.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""


MSIN = "tNDPPP-generic.MS"
PARMDB_TGZ = "tApplyCal2.parmdb.tgz"
PARMDB = "tApplyCal.parmdb"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar(f"{tcf.SRCDIR}/{PARMDB_TGZ}")


@pytest.mark.parametrize("updateweights", [False, True])
def test_with_and_without_weights_update(updateweights):
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=.",
            f"msout.datacolumn=DATA3",
            f"msout.weightcolumn=WEIGHTS_NEW",
            f"steps=[applycal]",
            f"applycal.parmdb={PARMDB}",
            f"showcounts=false",
            f"applycal.updateweights={updateweights}",
        ]
    )

    if updateweights:
        taql_commands = [
            f"select from {MSIN} where not(all(WEIGHTS_NEW~=81*WEIGHT_SPECTRUM))"
        ]
    else:
        taql_commands = [
            f"select from {MSIN} where not(all(DATA~=9*DATA3))",
            f"select from {MSIN} where not(all(WEIGHTS_NEW~=WEIGHT_SPECTRUM))",
        ]

    for taql_command in taql_commands:
        assert_taql(taql_command)


def test_common_scalar_phase():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=.",
            f"msout.datacolumn=DATA3",
            f"steps=[applycal]",
            f"applycal.parmdb={PARMDB}",
            f"applycal.correction=commonscalarphase",
            f"showcounts=false",
        ]
    )
    taql_command = f"select from {MSIN} where not(all(DATA~=DATA3))"
    assert_taql(taql_command)


def test_scalar_amplitude_values():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=.",
            f"msout.datacolumn=DATA3",
            f"steps=[applycal]",
            f"applycal.parmdb={PARMDB}",
            f"applycal.correction=scalaramplitude",
            f"showcounts=false",
        ]
    )
    taql_command = f"select from {MSIN} where not(all(DATA~=9*DATA3))"
    assert_taql(taql_command)


def test_rotation_angle():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=.",
            f"msout.datacolumn=DATA3",
            f"steps=[applycal]",
            f"applycal.parmdb=tApplyCal.parmdb",
            f"applycal.correction=rotationmeasure",
            f"showcounts=false",
        ]
    )


def test_full_jones():
    # Produce a trivial h5parm with full Jones solutions
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=",
            f"msout.datacolumn=DATA3",
            f"steps=[ddecal]",
            f"ddecal.h5parm=full-jones-solutions.h5",
            f"ddecal.mode=fulljones",
            f"ddecal.modeldatacolumns=[DATA]",
        ]
    )
    # apply it
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=",
            f"msout.datacolumn=DATA3",
            f"steps=[applycal]",
            f"applycal.parmdb=full-jones-solutions.h5",
            f"applycal.correction=fulljones",
            f"applycal.solset = sol000",
            f"applycal.soltab = [amplitude000,phase000]",
        ]
    )
