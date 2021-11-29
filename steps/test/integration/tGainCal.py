import pytest
import shutil
import os
import sys
import uuid
from subprocess import check_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, untar_ms

"""
Tests for gaincal (direction independent calibration).

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tGainCal.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSGAINCAL = "tGainCal.tab"
MSIN = "tNDPPP-generic.MS"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar_ms(f"{tcf.SRCDIR}/{MSGAINCAL}.tgz")

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
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={MSIN}/sky",
            "predict.usebeammodel=false",
        ]
    )


def test_caltype_diagonal(create_model_data):
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-diagonal",
            "gaincal.usebeammodel=false",
            "gaincal.caltype=diagonal",
            "gaincal.propagatesolutions=true",
            "gaincal.solint=1",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL",
            "steps=[applycal]",
            f"applycal.parmdb={MSIN}/inst-diagonal",
        ]
    )

    # Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima.
    taql_command = f"select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL-t1.MODEL_DATA)))) as bbsres from {MSIN} t1, {MSGAINCAL} t2) where dpppres>bbsres*1.02"
    assert_taql(taql_command)

    # Checking that not everything was flagged
    taql_command = (
        f"select from {MSIN} where all(FLAG) groupby true having gcount()>100"
    )
    assert_taql(taql_command)

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-diagonal-tpp",
            "gaincal.usebeammodel=false",
            "gaincal.caltype=diagonal",
            "gaincal.solint=4",
            "gaincal.timeslotsperparmupdate=1",
            "gaincal.propagatesolutions=false",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL_TPP",
            "steps=[applycal]",
            f"applycal.parmdb={MSIN}/inst-diagonal-tpp",
        ]
    )

    # "Comparing the difference between applying with timeslotsperparmupdate = default and timeslotsperparmupdate=1"
    taql_command = f"select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL-t1.MODEL_DATA)))) as bbsres from {MSIN} t1, {MSGAINCAL} t2) where dpppres>bbsres*1.02"
    assert_taql(taql_command)

    # Checking that not everything was flagged
    taql_command = (
        f"select from {MSIN} where all(FLAG) groupby true having gcount()>100"
    )
    assert_taql(taql_command)


def test_caltype_fulljones(create_model_data):
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_FULLJONES_GAINCAL",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-fulljones",
            "gaincal.usebeammodel=false",
            "gaincal.caltype=fulljones",
            "gaincal.solint=1",
            "gaincal.applysolution=true",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_FULLJONES",
            "steps=[applycal]",
            f"applycal.parmdb={MSIN}/inst-fulljones",
        ]
    )

    # Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima.
    taql_command = f"select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_FULLJONES-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_FULLJONES-t1.MODEL_DATA)))) as bbsres from {MSIN} t1, {MSGAINCAL} t2) where dpppres>bbsres*1.02"
    assert_taql(taql_command)

    # Comparing the solutions from gaincal + applycal with gaincal directly
    taql_command = f"select from tNDPPP-generic.MS where not(all(DPPP_FULLJONES_GAINCAL ~= DPPP_FULLJONES))"
    assert_taql(taql_command)

    # Checking that not everything was flagged
    taql_command = f"select from tNDPPP-generic.MS where all(FLAG) groupby true having gcount()>100"
    assert_taql(taql_command)


def test_caltype_diagonal_nchan_2(create_model_data):

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL_NCHAN_GAINCAL",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-diagonal-nchan",
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
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_DIAGONAL_NCHAN",
            "steps=[applycal]",
            f"applycal.parmdb={MSIN}/inst-diagonal-nchan",
        ]
    )

    # Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima.
    taql_command = f"select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL_NCHAN-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL_NCHAN-t1.MODEL_DATA)))) as bbsres from {MSIN} t1, {MSGAINCAL} t2) where dpppres>bbsres*1.02"
    assert_taql(taql_command)

    # Comparing the solutions from gaincal + applycal with gaincal directly
    taql_command = f"select from tNDPPP-generic.MS where not(all(DPPP_DIAGONAL_NCHAN_GAINCAL ~= DPPP_DIAGONAL_NCHAN))"
    assert_taql(taql_command)

    # Checking that not everything was flagged
    taql_command = f"select from tNDPPP-generic.MS where all(FLAG) groupby true having gcount()>100"
    assert_taql(taql_command)


@pytest.mark.parametrize("caltype", ["tec", "tecandphase"])
def test_output_parameters(caltype):

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DPPP_TEC",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-{caltype}",
            f"gaincal.caltype={caltype}",
            "gaincal.solint=2",
        ]
    )

    # For now, only testing that the right parameter names are in the output
    taql_command = f'select from {MSIN}/inst-{caltype} where (select NAME from ::NAMES)[NAMEID]=="TEC:CS001HBA0"'
    assert_taql(taql_command, 1)

    taql_command = f'select from {MSIN}/inst-{caltype} where (select NAME from ::NAMES)[NAMEID]=="CommonScalarPhase:CS001HBA0"'
    assert_taql(taql_command, 0 if caltype == "tec" else 1)


def test_filter():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=tNDPPP-filtered.MS",
            "steps=[filter,gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-filter",
            "filter.baseline='!CS001HBA0&&*'",
            "gaincal.baseline='!CS002HBA1,RS305HBA&&*'",
            "gaincal.caltype=diagonal",
        ]
    )

    taql_command = f'select from {MSIN}/inst-filter::NAMES where NAME LIKE "CS001HBA0%" OR NAME LIKE "%CS002HBA1%" OR NAME LIKE "%RS305HBA%"'
    assert_taql(taql_command)


def test_debug_output():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "numthreads=1",
            "steps=[gaincal]",
            f"gaincal.sourcedb={MSIN}/sky",
            f"gaincal.parmdb={MSIN}/inst-debug",
            "gaincal.caltype=diagonal",
            "gaincal.debuglevel=1",
        ]
    )

    assert os.path.exists("debug.h5")
