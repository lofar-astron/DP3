# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
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
from utils import assert_taql, untar_ms

"""
Replacement for tPredict.sh using pytest.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tPredict.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSPREDICT = "tPredict.tab"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar_ms(f"{tcf.SRCDIR}/{MSPREDICT}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


def test_with_beam_subtract():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={MSIN}/sky",
            "predict.usebeammodel=true",
            "predict.operation=subtract",
        ]
    )

    # Compare the MODEL_DATA column of the output MS with the original data minus the BBS reference output.
    taql_command = f"select t1.MODEL_DATA, t1.DATA-t2.PREDICT_beam, abs(t1.MODEL_DATA/(t1.DATA-t2.PREDICT_beam)-1) from {MSIN} t1, {MSPREDICT} t2 where not all(near(t1.MODEL_DATA,t1.DATA-t2.PREDICT_beam,5e-2) || (isnan(t1.DATA) && isnan(t2.PREDICT_beam)) || nearabs(t2.PREDICT_beam, 0, 1e-5))"
    assert_taql(taql_command)


def test_without_beam_add():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={MSIN}/sky",
            "predict.usebeammodel=false",
            "predict.operation=add",
        ]
    )

    # Compare the MODEL_DATA column of the output MS with the original data plus the BBS reference output.
    taql_command = f"select from {MSIN} t1, {MSPREDICT} t2 where not all(near(t1.MODEL_DATA,t1.DATA+t2.PREDICT_nobeam,5e-2) || (isnan(t1.DATA) && isnan(t2.PREDICT_nobeam)))"
    assert_taql(taql_command)


# This sub-test is disabled since it uses parmdbm, which is deprecated.
# This code was copied from tPredict.sh but not converted to python yet.
#
# echo; echo "Test without beam, with applycal, subtract (like peeling)"; echo
# parmdbm <<EOL
# open table="tPredict.parmdb"
# adddef Gain:0:0:Real values=3.
# adddef Gain:1:1:Real values=3.
# EOL
# cmd="$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=MODEL_DATA steps=[predict] predict.sourcedb=tNDPPP-generic.MS/sky predict.applycal.parmdb=tPredict.parmdb predict.operation=subtract"
# echo $cmd
# $cmd
## Compare the MODEL_DATA column of the output MS with the original data minus the BBS reference output.
# taqlcmd='select from tNDPPP-generic.MS t1, tPredict.tab t2 where not all(near(t1.MODEL_DATA,t1.DATA-9.0*t2.PREDICT_nobeam,5e-2) || (isnan(t1.DATA) && isnan(t2.PREDICT_nobeam)))'
# echo $taqlcmd
# $taqlexe $taqlcmd > taql.out
# diff taql.out taql.ref  ||  exit 1


@pytest.mark.parametrize("use_beam", [False, True])
def test_without_and_with_beam(use_beam):
    predict_column = "PREDICT_beam" if use_beam else "PREDICT_nobeam"

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={MSIN}/sky",
            f"predict.usebeammodel={'true' if use_beam else 'false'}",
        ]
    )
    taql_command = f"select from {MSIN} t1, {MSPREDICT} t2 where not all(near(t1.MODEL_DATA,t2.{predict_column},5e-2) || (isnan(t1.MODEL_DATA) && isnan(t2.{predict_column})))"
    assert_taql(taql_command)


@pytest.mark.parametrize("use_time_smearing", [False, True])
def test_without_and_with_time_smearing(use_time_smearing):
    sourcedb = "timesmearing.sourcedb"

    # Create a point-source model with a declination offset of 2 degrees.
    # The phase center of tNDPPP-generic.MS is 1h37m41.299, +033d09m35.132.
    dec_offset = 2  # Declination offset in degrees
    with open("timesmearing.skymodel", "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\n"
        )
        f.write(
            f"center, POINT, 01:37:41.299, +{33 + dec_offset}.09.35.132, 10, , , , ,\n"
        )

    shutil.rmtree(sourcedb, ignore_errors=True)
    check_call([tcf.MAKESOURCEDBEXE, "in=timesmearing.skymodel", f"out={sourcedb}"])
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[average,predict]",
            "average.timestep=3",
            f"predict.sourcedb={sourcedb}",
            f"msout=ts-{'on' if use_time_smearing else 'off'}.MS",
            "predict.correcttimesmearing=10" if use_time_smearing else "",
        ]
    )
    if not use_time_smearing:
        # Without time smearing, XX,XY,YX,YY should always be 10,0,0,10
        taql_command = (
            "select from ts-off.MS "
            "where all(near(abs(DATA[,0]),10,1e-6))"
            "  and all(near(abs(DATA[,1]), 0,1e-6))"
            "  and all(near(abs(DATA[,2]), 0,1e-6))"
            "  and all(near(abs(DATA[,3]),10,1e-6))"
        )
        assert_taql(taql_command, 56)
    else:
        # With time smearing, XX,XY,YX,YY  should be close to (9.34, 0, 0, 9.34)
        # for all channels, when using the longest baseline (ant1=3, ant2=5).
        taql_command = (
            "select from ts-on.MS where ANTENNA1==3 and ANTENNA2==5"
            " and all(near(abs(DATA[,0]),9.34,5e-4))"
            " and all(near(abs(DATA[,1]),   0,1e-6))"
            " and all(near(abs(DATA[,2]),   0,1e-6))"
            " and all(near(abs(DATA[,3]),9.34,5e-4))"
        )
        assert_taql(taql_command, 2)
