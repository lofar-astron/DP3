# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call

import numpy as np
import pytest

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Replacement for tPredict.sh using pytest.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tPredict.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSPREDICT = "tPredict.tab"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar(f"{tcf.SRCDIR}/{MSPREDICT}.tgz")


@pytest.fixture
def make_h5parm():
    """
    Create minimal H5Parm with scalar value. By omitting all of the dimensions
    time, ant, freq, dir, pol, the value is broadcast along all of those
    dimensions. So all times, frequencies, antennas get the value 3.
    """

    import h5py  # Don't import h5py when pytest is only collecting tests

    with h5py.File("test.h5", "w") as h5file:
        soltab = h5file.create_group("solset000/soltab000")
        soltab.attrs["TITLE"] = np.void(b"amplitude")
        dataset = soltab.create_dataset("val", shape=(1), dtype=float)
        dataset.attrs["AXES"] = np.void(b"none")
        h5file["solset000/soltab000/val"][:] = 3.0
        weights = soltab.create_dataset("weight", shape=(1), dtype=float)
        h5file["solset000/soltab000/weight"][:] = 1.0
        weights.attrs["AXES"] = np.void(b"none")


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


def test_with_applycal(make_h5parm):
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={MSIN}/sky",
            "predict.applycal.parmdb=test.h5",
            "predict.applycal.correction=soltab000",
            "predict.operation=subtract",
        ]
    )
    # Compare the MODEL_DATA column of the output MS with the original data minus the BBS reference output.
    taql_command = f"select from {MSIN} t1, {MSPREDICT} t2 where not all(near(t1.MODEL_DATA,t1.DATA-9.0*t2.PREDICT_nobeam,5e-2) || (isnan(t1.DATA) && isnan(t2.PREDICT_nobeam)))"
    assert_taql(taql_command)


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

    check_call(
        [tcf.MAKESOURCEDBEXE, "in=timesmearing.skymodel", f"out={sourcedb}"]
    )
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


def test_time_smearing_with_msupdate():
    """
    Test that DP3 can update an MS after running a Predict step that does
    time smearing, which internally changes the meta data.
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={MSIN}/sky",
            "predict.correcttimesmearing=2",
        ]
    )


@pytest.mark.parametrize("use_beam", [False, True])
def test_without_and_with_beam_parallelbaseline(use_beam):
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
            "predict.parallelbaselines=true",
        ]
    )
    taql_command = f"select from {MSIN} t1, {MSPREDICT} t2 where not all(near(t1.MODEL_DATA,t2.{predict_column},5e-2) || (isnan(t1.MODEL_DATA) && isnan(t2.{predict_column})))"
    assert_taql(taql_command)
