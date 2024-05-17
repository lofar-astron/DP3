# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
from subprocess import check_call

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Similar to tPredict.py, testing SagecalPredict using pytest.
"""

# Note we re-use the same MS for output
MSIN = "tNDPPP-generic.MS"
# But we have a different reference MS for comparison
MSPREDICT = "tSagecalPredict.tab"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar(f"{tcf.SRCDIR}/{MSPREDICT}.tgz")


def test_with_beam_replace():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            "predict.type=sagecalpredict",
            f"predict.sourcedb={MSIN}/sky",
            "predict.usebeammodel=true",
            "predict.operation=replace",
        ]
    )

    # Compare the MODEL_DATA column of the output MS with the PREDICT_beam of reference MS.
    taql_command = f"select t1.MODEL_DATA, t2.PREDICT_beam from {MSIN} t1, {MSPREDICT} t2 where not all(near(t1.MODEL_DATA,t2.PREDICT_beam,5e-2) || (isnan(t1.MODEL_DATA) && isnan(t2.PREDICT_beam)))"
    assert_taql(taql_command)


def test_with_rapthor_workflow():
    """
    This is a long test, but it tries to cover all modes of use, includeing
    in DDECal and in h5parmpredict
    1) scalar phace cal 2) complex gain cal 3) apply solutions and predict
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            "ddecal.type=ddecal",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.usebeammodel=true",
            "ddecal.beammode=array_factor",
            "ddecal.sagecalpredict=true",
            "ddecal.mode=scalarphase",
            "ddecal.maxiter=50",
            "ddecal.nchan=8",
            "ddecal.stepsize=1e-3",
            "ddecal.solint=10",
            "ddecal.h5parm=solutions_stage1.h5",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            "ddecal.type=ddecal",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.usebeammodel=true",
            "ddecal.beammode=array_factor",
            "ddecal.sagecalpredict=true",
            "ddecal.mode=complexgain",
            "ddecal.maxiter=50",
            "ddecal.nchan=8",
            "ddecal.stepsize=1e-3",
            "ddecal.solint=10",
            "ddecal.h5parm=solutions_stage2.h5",
            "ddecal.applycal.steps=[fastphase]",
            "ddecal.applycal.fastphase.correction=phase000",
            "ddecal.applycal.parmdb=solutions_stage1.h5",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            "predict.type=sagecalpredict",
            f"predict.sourcedb={MSIN}/sky",
            "predict.usebeammodel=true",
            "predict.beammode=array_factor",
            "predict.operation=replace",
            "predict.applycal.correction=phase000",
            "predict.applycal.parmdb=solutions_stage2.h5",
        ]
    )

    # Compare the MODEL_DATA column of the output MS with the PREDICT_nobeam of reference MS.
    taql_command = f"select t1.MODEL_DATA, t2.PREDICT_nobeam from {MSIN} t1, {MSPREDICT} t2 where not all(near(t1.MODEL_DATA,t2.PREDICT_nobeam,5e-2) || (isnan(t1.MODEL_DATA) && isnan(t2.PREDICT_nobeam)))"
    assert_taql(taql_command)
