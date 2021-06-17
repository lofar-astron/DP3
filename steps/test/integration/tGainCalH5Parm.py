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
Tests for applying the calibration.

Script can be invoked in two ways:
- as standalone from the build/steps/tests/integration directory,
  using `pytest source/tApplyCal2.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""


MSIN = "tNDPPP-generic.MS"
REF_DATA = "tGainCal.tab"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar_ms(f"{tcf.SRCDIR}/{REF_DATA}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


@pytest.fixture()
def create_model_data():
    # Creating MODEL_DATA so that residual can be computed
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "showprogress=false",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={MSIN}/sky",
            "predict.usebeammodel=false",
        ]
    )


def test_diagonal(create_model_data):
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-diagonal.h5",
            "gaincal.usebeammodel=false",
            "gaincal.caltype=diagonal",
            "gaincal.propagatesolutions=true",
            "gaincal.solint=1",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL",
            "steps=[applycal]",
            f"applycal.parmdb={MSIN}/inst-diagonal.h5",
            "applycal.steps=[amplitude,phase]",
            "applycal.phase.correction=phase000",
            "applycal.amplitude.correction=amplitude000",
            "applycal.amplitude.correction=amplitude000",
        ]
    )

    # 1 - Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima.
    # 2 - Checking that not everything was flagged
    taql_commands = [
        f"select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL-t1.MODEL_DATA)))) as bbsres from {MSIN} t1, {REF_DATA} t2) where dpppres>bbsres*1.02",
        f"select from {MSIN} where all(FLAG) groupby true having gcount()>100",
    ]

    for taql_command in taql_commands:
        assert_taql(taql_command)


def test_fulljones(create_model_data):
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_FULLJONES_GAINCAL",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-fulljones.h5",
            "gaincal.usebeammodel=false",
            "gaincal.caltype=fulljones",
            "gaincal.solint=1",
            "gaincal.applysolution=true",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL_NCHAN_GAINCAL",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-diagonal-nchan.h5",
            "gaincal.usebeammodel=false",
            "gaincal.caltype=diagonal",
            "gaincal.solint=4",
            "gaincal.nchan=2",
            "gaincal.applysolution=true",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL_NCHAN",
            "steps=[applycal]",
            f"applycal.parmdb={MSIN}/inst-diagonal-nchan.h5",
            "applycal.steps=[phase,amplitude]",
            "applycal.phase.correction=phase000",
            "applycal.amplitude.correction=amplitude000",
        ]
    )

    # 1 - Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima.
    # 2 - Comparing the solutions from gaincal + applycal with gaincal directly
    # 3 - Checking that not everything was flagged
    taql_commands = [
        f"select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL_NCHAN-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL_NCHAN-t1.MODEL_DATA)))) as bbsres from {MSIN} t1, {REF_DATA} t2) where dpppres>bbsres*1.02",
        f"select from {MSIN} where not(all(DPPP_DIAGONAL_NCHAN_GAINCAL ~= DPPP_DIAGONAL_NCHAN))",
        f"select from {MSIN} where all(FLAG) groupby true having gcount()>100",
    ]

    for taql_command in taql_commands:
        assert_taql(taql_command)


def test_diagonal_nchan():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL_NCHAN_7_GAINCAL",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-diagonal-nchan7.h5",
            "gaincal.usebeammodel=false",
            "gaincal.caltype=diagonal",
            "gaincal.solint=4",
            "gaincal.nchan=2",
            "gaincal.applysolution=true",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL_NCHAN_7",
            "steps=[applycal]",
            f"applycal.parmdb={MSIN}/inst-diagonal-nchan7.h5",
            "applycal.steps=[amplitude,phase]",
            "applycal.amplitude.correction=amplitude000",
            "applycal.phase.correction=phase000",
        ]
    )

    # Comparing the solutions from gaincal + applycal with gaincal directly
    taql_command = f"select from {MSIN} where not(all(DPPP_DIAGONAL_NCHAN_7_GAINCAL ~= DPPP_DIAGONAL_NCHAN_7))"
    assert_taql(taql_command)


@pytest.mark.parametrize("caltype", ["tec", "tecandphase"])
def test_caltype(caltype):
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_TEC",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-{caltype}.h5",
            f"gaincal.caltype={caltype}",
            "gaincal.solint=2",
        ]
    )


def test_filter():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=tNDPPP-filtered.MS",
            "steps=[filter,gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-filter.h5",
            "filter.baseline='!CS001HBA0&&*'",
            "gaincal.baseline='!CS002HBA1,RS305HBA&&*'",
            "gaincal.caltype=diagonal",
        ]
    )
