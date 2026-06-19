# Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from subprocess import check_call

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, get_taql_result, run_dp3, run_in_tmp_path, untar

"""
Integration tests for the Combine step.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tCombine.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN_REGULAR = (
    "tNDPPP-generic.MS"  # This MS has 8 freq channels and 6 time steps
)
MSIN_BDA = "tNDPPP-bda.MS"
BUFFERNAME = "test_buffer"
SKYMODEL = f"{tcf.RESOURCEDIR}/tNDPPP-generic-skymodel.txt"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN_REGULAR}.tgz")
    untar(f"{tcf.RESOURCEDIR}/{MSIN_BDA}.tgz")


def test_subtract():
    """
    Test whether predicting data and then subtracting the same data from a
    named buffer yields zero.
    """
    check_call(  # Predict into buffer and into named buffer, subtract the two
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN_REGULAR}",
            "msout=.",
            "steps=[predict,predict2,combine]",
            f"predict.sourcedb={SKYMODEL}",
            f"predict2.sourcedb={SKYMODEL}",
            f"predict.outputmodelname={BUFFERNAME}",
            "combine.operation=subtract",
            f"combine.buffername={BUFFERNAME}",
        ]
    )

    max_diff = 1e-4 if tcf.USE_FAST_PREDICT else 0

    taql_command = (
        f"select from {MSIN_REGULAR} where sum(abs(DATA)) > {max_diff}"
    )
    assert_taql(taql_command, 0)


def test_add():
    """
    Test whether predicting data and then subtracting the same data from a
    named buffer yields zero.
    """
    check_call(  # Predict into temporary column
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN_REGULAR}",
            "msout=.",
            "msout.datacolumn=PREDICTED_DATA",
            "steps=[predict]",
            f"predict.sourcedb={SKYMODEL}",
        ]
    )

    check_call(  # Predict into buffer and named buffer, add them
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN_REGULAR}",
            "msout=.",
            "steps=[predict,predict2,combine]",
            f"predict.sourcedb={SKYMODEL}",
            f"predict2.sourcedb={SKYMODEL}",
            f"predict.outputmodelname={BUFFERNAME}",
            "combine.operation=add",
            f"combine.buffername={BUFFERNAME}",
        ]
    )

    max_diff = 1e-4 if tcf.USE_FAST_PREDICT else 1e-10

    taql_command = f"select from {MSIN_REGULAR} where sum(abs(DATA - 2 * PREDICTED_DATA)) > {max_diff}"
    assert_taql(taql_command, 0)


def test_bda_add():
    """
    Test whether predicting data and then subtracting the same data from a
    named buffer yields zero.
    """
    run_dp3(  # Predict into temporary ms
        [
            "checkparset=1",
            f"msin={MSIN_BDA}",
            "msout=reference.ms",
            "msout.overwrite=True",
            "steps=[predict]",
            f"predict.sourcedb={SKYMODEL}",
        ]
    )

    run_dp3(  # Predict into buffer and named buffer, add them
        [
            "checkparset=1",
            f"msin={MSIN_BDA}",
            "msout=result.ms",
            "msout.overwrite=True",
            "steps=[predict,predict2,combine]",
            f"predict.sourcedb={SKYMODEL}",
            f"predict2.outputmodelname={BUFFERNAME}",
            f"predict2.sourcedb={SKYMODEL}",
            "combine.operation=add",
            f"combine.buffername={BUFFERNAME}",
        ]
    )

    max_diff = 1e-4 if tcf.USE_FAST_PREDICT else 1e-10

    taql_command = f"select from reference.ms as ref, result.ms as result where sum(abs(result.DATA - 2 * ref.DATA)) > {max_diff}"
    assert_taql(taql_command, 0)


def test_combine_with_transfer():
    """
    Test whether the the combination [transfer, combine] works.
    """
    run_dp3(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN_REGULAR}",
            "msout=avg.ms",
            "steps=[averager]",
            "averager.timestep=2",
            "averager.freqstep=4",
        ]
    )

    run_dp3(
        [
            f"msin={MSIN_REGULAR}",
            "msout=.",
            "steps=[transfer,combine]",
            "transfer.data=True",
            "transfer.source_ms=avg.ms",
            f"transfer.outputbuffername={BUFFERNAME}",
            "combine.operation=subtract",
            f"combine.buffername={BUFFERNAME}",
        ]
    )
