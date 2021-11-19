# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

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
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
- using `pytest source/tBdaPredict.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""


MSIN_REGULAR = "tNDPPP-generic.MS"
MSIN_BDA = "tNDPPP-bda.MS"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN_REGULAR}.tgz")
    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN_BDA}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


@pytest.fixture()
def create_skymodel():
    with open("test.skymodel", "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation\r\n"
        )
        f.write(
            "point-0, POINT, 0.4362457236387493, 0.5287469737178224, 1.0, 0, 0, 0, , , \r\n"
        )

    check_call(
        [tcf.MAKESOURCEDBEXE, "in=test.skymodel", "out=test.sourcedb", "append=false"]
    )


@pytest.fixture()
def create_skymodel_in_phase_center():
    with open("test.skymodel", "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\r\n"
        )
        f.write(f"center, POINT, 01:37:41.299000, +33.09.35.132000, 10, , , , , \r\n")

    check_call(
        [tcf.MAKESOURCEDBEXE, "in=test.skymodel", "out=test.sourcedb", "append=false"]
    )


def test_bdapredict(create_skymodel):
    common_args = [
        "bdaaverager.timebase=20000",
        "bdaaverager.frequencybase=20000",
        "bdaaverager.maxinterval=61",
        "predict.sourcedb=test.sourcedb",
        "predict.usebeammodel=F",
    ]

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_REGULAR}",
            "msout=bdapredict0.MS",
            "steps=[predict,bdaaverager]",
        ]
        + common_args
    )
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_REGULAR}",
            "msout=bdapredict1.MS",
            "steps=[bdaaverager, predict]",
        ]
        + common_args
    )

    # Compare the DATA columns of the output MSs.
    # Because of flagging and weighting the difference can be as large
    # the difference between visibilities at the edge and centre of the interval.
    # For a scenario without flagging and with uniform weighting the tolerance can be lower.
    taql_command = "select ANTENNA1, ANTENNA2 from bdapredict0.MS t1, bdapredict1.MS t2 where not all(abs(t1.DATA-t2.DATA)<15e-2 || t1.FLAG || t1.WEIGHT_SPECTRUM==0)"
    assert_taql(taql_command)


def test_predicted_values_regular_input(create_skymodel_in_phase_center):

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_REGULAR}",
            "msout=bdapredict.MS",
            "msout.overwrite=true",
            "steps=[bdaaverager, predict]",
            "bdaaverager.timebase=600",
            "bdaaverager.frequencybase=1000",
            "predict.sourcedb=test.sourcedb",
            "predict.usebeammodel=F",
        ]
    )

    taql_command = (
        "select from bdapredict.MS "
        "where all(near(abs(DATA[,0]),10,1e-6))"
        "  and all(near(abs(DATA[,1]), 0,1e-6))"
        "  and all(near(abs(DATA[,2]), 0,1e-6))"
        "  and all(near(abs(DATA[,3]),10,1e-6))"
    )
    assert_taql(taql_command, 153)


def test_predicted_values_bda_input(create_skymodel_in_phase_center):
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_BDA}",
            "msout=bdapredict.MS",
            "msout.overwrite=true",
            "steps=[predict]",
            "predict.sourcedb=test.sourcedb",
            "predict.usebeammodel=F",
        ]
    )

    taql_command = (
        "select from bdapredict.MS "
        "where all(near(abs(DATA[,0]),10,1e-6))"
        "  and all(near(abs(DATA[,1]), 0,1e-6))"
        "  and all(near(abs(DATA[,2]), 0,1e-6))"
        "  and all(near(abs(DATA[,3]),10,1e-6))"
    )
    assert_taql(taql_command, 159)
