# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from subprocess import check_call

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
- using `pytest source/tBdaExpander.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

# These tests test that DP3 supports different bdaexpander scenerios.
# They only test that the various steps integrate well and do not check
# the output of each run

MSIN_REGULAR = "tNDPPP-generic.MS"
MSIN_BDA = "tNDPPP-bda.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN_REGULAR}.tgz")
    untar(f"{tcf.RESOURCEDIR}/{MSIN_BDA}.tgz")


@pytest.fixture()
def create_skymodel():
    with open("test.skymodel", "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\r\n"
        )
        f.write(
            "center, POINT, 16:38:28.205000, + 63.44.34.314000, 10, , , , , \r\n"
        )
        f.write(
            "ra_off, POINT, 16:38:28.205000, + 64.44.34.314000, 10, , , , , \r\n"
        )
        f.write(
            "radec_off, POINT, 16:38:28.205000, +65.44.34.314000, 10, , , , , \r\n"
        )

    check_call([tcf.MAKESOURCEDBEXE, "in=test.skymodel", "out=test.sourcedb"])


def test_only_expand():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_BDA}",
            "msout=out.MS",
            "msout.uvwcompression=false",  # TODO why is this necessary?
            "steps=[bdaexpander]",
        ]
    )


def test_expand_average():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_BDA}",
            "msout=out.MS",
            "msout.overwrite=true",
            "steps=[bdaexpander, bdaaverager]",
            "bdaaverager.frequencybase=1000",
            "bdaaverager.timebase=100",
        ]
    )

    taql_check_visibilities = f"select from(select gsumsqr(abs(t_in.DATA[isnan(t_in.DATA)]-t_out.DATA[isnan(t_out.DATA)])) as diff from {MSIN_BDA} t_in, out.MS t_out) where diff>1.e-6"
    assert_taql(taql_check_visibilities)


def test_bdaaverager_ddecal_bdaexpander(create_skymodel):
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_REGULAR}",
            "msout=out.MS",
            "msout.overwrite=true",
            "steps=[bdaaverager, ddecal, bdaexpander]",
            "ddecal.onlypredict=true",
            "ddecal.directions=[[center],[ra_off],[radec_off]]",
            "ddecal.sourcedb=test.sourcedb",
        ]
    )


def test_bdaexpander_ddecal(create_skymodel):
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_BDA}",
            "msout=out.MS",  # msout=. does not work -> documented
            "msout.overwrite=true",
            "msout.uvwcompression=false",  # TODO why is this necessary?
            "steps=[bdaexpander, ddecal]",
            "ddecal.directions=[[center],[ra_off],[radec_off]]",
            "ddecal.sourcedb=test.sourcedb",
        ]
    )


def test_regular_buffer_writing():
    # Checks that after writing a regular MS after the
    # bdaexpander step, DP3 can read that regular MS.
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            f"msin={MSIN_BDA}",
            "msout=regular_buffer.MS",
            "msout.uvwcompression=false",  # TODO why is this necessary?
            "steps=[bdaexpander]",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=true",
            "msin=regular_buffer.MS",
            "msout=out.MS",
            "msout.uvwcompression=false",
            "steps=[]",
        ]
    )
