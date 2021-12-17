# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import shutil
import uuid
from subprocess import check_call, check_output, CalledProcessError, STDOUT
import numpy as np

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from testconfig import TAQLEXE
from utils import assert_taql, untar_ms

"""
Script can be invoked in two ways:
- as standalone from the build/ddecal/test/integration directory,
  using `pytest source/tDDECal.py` (extended with pytest options of your choice)
- using ctest, see DP3/ddecal/test/integration/CMakeLists.txt
"""

MSINTGZ = "tDDECal.in_MS.tgz"
MSIN = "tDDECal.MS"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSINTGZ}")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


@pytest.fixture()
def copy_data_to_model_data():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[]",
        ]
    )


@pytest.fixture()
def create_corrupted_visibilities():
    taqlcommand = f"update {MSIN} set WEIGHT_SPECTRUM=1, FLAG=False"
    check_output([TAQLEXE, "-noph", taqlcommand])

    # Use ddecal to create template h5parm
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            "ddecal.sourcedb=tDDECal.MS/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrumentcorrupted.h5",
            "ddecal.mode=complexgain",
        ]
    )

    # Modify h5 file for multiple solution intervals
    import h5py  # Don't import h5py when pytest is only collecting tests.

    h5file = h5py.File("instrumentcorrupted.h5", "r+")
    sol = h5file["sol000/amplitude000/val"]
    sol[:4, ..., 0, :] = np.sqrt(5)
    sol[4:, ..., 0, :] = np.sqrt(5 + 2)
    sol[:4, ..., 1, :] = np.sqrt(9)
    sol[4:, ..., 1, :] = np.sqrt(9 + 2)
    sol[:4, ..., 2, :] = np.sqrt(13)
    sol[4:, ..., 2, :] = np.sqrt(13 + 2)
    h5file.close()

    # Predict corrupted visibilities into DATA column
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DATA",
            "steps=[h5parmpredict]",
            f"h5parmpredict.sourcedb={MSIN}/sky",
            "h5parmpredict.applycal.parmdb=instrumentcorrupted.h5",
            "h5parmpredict.applycal.correction=amplitude000",
        ]
    )


@pytest.mark.parametrize(
    "caltype", ["complexgain", "scalarcomplexgain", "amplitudeonly", "scalaramplitude"]
)
@pytest.mark.parametrize("solint", [0, 1, 2, 4])
@pytest.mark.parametrize("nchan", [1, 2, 5])
def test(
    create_corrupted_visibilities, copy_data_to_model_data, caltype, solint, nchan
):
    # Subtract corrupted visibilities using multiple predict steps
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            f"ddecal.solint={solint}",
            f"ddecal.nchan={nchan}",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
        ]
    )

    # Calibrate on the original sources, caltype=$caltype
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3]",
            f"predict1.sourcedb={MSIN}/sky",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center,dec_off]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={MSIN}/sky",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[radec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={MSIN}/sky",
            "predict3.applycal.parmdb=instrument.h5",
            "predict3.sources=[ra_off]",
            "predict3.operation=subtract",
            "predict3.applycal.correction=amplitude000",
        ]
    )

    if solint == 0:
        tolerance = 0.15
    else:
        tolerance = 0.015
    taqlcommand_run = f"select norm_residual/norm_data FROM (select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*DATA))) as norm_data, sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*SUBTRACTED_DATA))) as norm_residual from {MSIN})"
    check_output([TAQLEXE, "-noph", taqlcommand_run])
    taql_command = f"select FROM (select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*DATA))) as norm_data, sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*SUBTRACTED_DATA))) as norm_residual from {MSIN}) where norm_residual/norm_data > {tolerance} or isinf(norm_residual/norm_data) or isnan(norm_residual/norm_data)"
    assert_taql(taql_command)


def test_h5parm_predict():
    # make calibration solutions
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode=diagonal",
        ]
    )

    # subtract using multiple predict steps
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3]",
            f"predict1.sourcedb={MSIN}/sky",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center,dec_off]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={MSIN}/sky",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[radec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={MSIN}/sky",
            "predict3.applycal.parmdb=instrument.h5",
            "predict3.sources=[ra_off]",
            "predict3.operation=subtract",
            "predict3.applycal.correction=amplitude000",
        ]
    )

    # subtract using h5parmpredict
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA_H5PARM",
            "steps=[h5parmpredict]",
            f"h5parmpredict.sourcedb={MSIN}/sky",
            "h5parmpredict.applycal.parmdb=instrument.h5",
            "h5parmpredict.operation=subtract",
            "h5parmpredict.applycal.correction=amplitude000",
        ]
    )

    # echo "Check that h5parmpredict creates the same output as multiple predict steps"
    taql_command = f"select from (select abs(gsumsqr(SUBTRACTED_DATA-SUBTRACTED_DATA)) as diff from {MSIN}) where diff>1.e-6"
    assert_taql(taql_command)


def test_pre_apply():
    # make calibration solutions
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode=diagonal",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3]",
            f"predict1.sourcedb={MSIN}/sky",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center,dec_off]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={MSIN}/sky",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[radec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={MSIN}/sky",
            "predict3.applycal.parmdb=instrument.h5",
            "predict3.sources=[ra_off]",
            "predict3.operation=subtract",
            "predict3.applycal.correction=amplitude000",
        ]
    )

    # Check that preapply runs (output is not tested)
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.applycal.parmdb=instrument.h5",
            "ddecal.applycal.steps=applyampl",
            "ddecal.applycal.applyampl.correction=amplitude000",
            "ddecal.h5parm=instrument2.h5",
            "ddecal.mode=scalarcomplexgain",
        ]
    )


def test_check_tec():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument-tec.h5",
            "ddecal.mode=tec",
        ]
    )


def test_check_tec_and_phase():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msin.baseline='!CS001HBA0'",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument-tecandphase.h5 ddecal.mode=tecandphase",
        ]
    )


@pytest.mark.parametrize(
    "solutions_per_direction", [None, [1], [2, 3, 1], [5, 5, 5], [2, 0]]
)
def test_dd_solution_intervals(solutions_per_direction):
    base_command = [
        tcf.DP3EXE,
        "checkparset=1",
        "numthreads=1",
        f"msin={MSIN}",
        "msout=.",
        "steps=[ddecal]",
        f"ddecal.sourcedb={MSIN}/sky",
        "ddecal.solint=6",
        "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
        "ddecal.h5parm=instrument-tec.h5",
        "ddecal.mode=scalaramplitude",
        "ddecal.solveralgorithm=directioniterative",
    ]
    try:
        check_output(
            base_command
            if solutions_per_direction is None
            else base_command
            + [f"ddecal.solutions_per_direction={solutions_per_direction}"],
            stderr=STDOUT,
        )
    except CalledProcessError as e:
        if solutions_per_direction == [5, 5, 5]:
            if not e.output.decode().startswith(
                "\nstd exception detected: Values in ddecal"
            ):
                raise e
        elif solutions_per_direction == [2, 0]:
            if not e.output.decode().startswith(
                "\nstd exception detected: All entries in ddecal"
            ):
                raise e
        else:
            raise (e)


def test_modelnextsteps(copy_data_to_model_data):
    # Multiply MODEL_DATA by 42
    taqlcommand_run = f"update {MSIN} set MODEL_DATA=DATA*42"
    check_output([TAQLEXE, "-noph", taqlcommand_run])

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[ddecal]",
            "ddecal.modeldatacolumns=[MODEL_DATA]",
            "ddecal.modelnextsteps.MODEL_DATA=[scaledata]",
            "ddecal.h5parm=instrument-modeldata",
            "ddecal.solint=2",
            "ddecal.nchan=3",
            "scaledata.stations='*'",
            "scaledata.scalesize=False",
            "scaledata.coeffs=1",
            "ddecal.subtract=True",
        ]
    )
    taql_command = f"select from (select abs(sumsqr(SUBTRACTED_DATA)/sumsqr(DATA)) as diff from {MSIN}) where diff>1.e-6"
    assert_taql(taql_command)
