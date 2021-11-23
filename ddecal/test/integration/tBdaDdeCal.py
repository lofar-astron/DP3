# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import shutil
import uuid
import numpy as np
from subprocess import check_call, check_output, run

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, untar_ms

"""
Script can be invoked in two ways:
- as standalone from the build/ddecal/test/integration directory,
  using `pytest source/tBdaDdeCal.py` (extended with pytest options of your choice)
- using ctest, see DP3/ddecal/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-bda.MS"
CWD = os.getcwd()
CORRUPTIONS = 3, 4, 7 # Corruption gain factors per antenna


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


@pytest.fixture()
def create_skymodel():
    with open("test.skymodel", "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\r\n"
        )
        f.write("center, POINT, 16:38:28.205000, +63.44.34.314000, 1, , , , , \r\n")
        f.write("ra_off, POINT, 16:58:28.205000, +63.44.34.314000, 1, , , , , \r\n")
        f.write("radec_off, POINT, 16:38:28.205000, +65.44.34.314000, 1, , , , , \r\n")


@pytest.fixture()
def create_corrupted_data():
    with open("test_corrupted.txt", "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\r\n"
        )
        f.write(f"center, POINT, 16:38:28.205000, +63.44.34.314000, {CORRUPTIONS[0] * CORRUPTIONS[0]}, , , , , \r\n")
        f.write(f"ra_off, POINT, 16:58:28.205000, +63.44.34.314000, {CORRUPTIONS[1] * CORRUPTIONS[1]}, , , , , \r\n")
        f.write(f"radec_off, POINT, 16:38:28.205000, +65.44.34.314000, {CORRUPTIONS[2] * CORRUPTIONS[2]}, , , , , \r\n")

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=corrupted.MS",
            "steps=[predict]",
            "predict.sourcedb=test_corrupted.txt",
            "numthreads=1",
        ]
    )

    check_call(
        [
            tcf.TAQLEXE,
            "update corrupted.MS set WEIGHT_SPECTRUM=1"
        ]
    )


def test_only_predict(create_skymodel):
    """Test that the right patches are summed with predict_only"""

    common_args = [
        "checkparset=1",
        f"msin={MSIN}",
        "msout.overwrite=true",
        "numthreads=1",
    ]

    predict_args = [
        "steps=[predict]",
        "predict.sourcedb=test.skymodel",
    ]

    check_call(
        [
            tcf.DP3EXE,
            "msout=BDADDECal_onlypredict.MS",
            "steps=[ddecal]",
            "ddecal.onlypredict=true",
            "ddecal.directions=[[center, dec_off],[ra_off],[radec_off]]",
            "ddecal.sourcedb=test.skymodel",
        ]
        + common_args
    )

    check_call(
        [
            tcf.DP3EXE,
            "msout=PREDICT_DIR_1.MS",
            "predict.sources=[center, dec_off]"
        ]
        + common_args
        + predict_args
    )

    check_call(
        [
            tcf.DP3EXE,
            "msout=PREDICT_DIR_2.MS",
            "predict.sources=[ra_off]"
        ]
        + common_args
        + predict_args
    )

    check_call(
        [
            tcf.DP3EXE,
            "msout=PREDICT_DIR_3.MS",
            "predict.sources=[radec_off]",
        ]
        + common_args
        + predict_args
    )

    taql_check_visibilities = (
        "select from ("
        "select gsumsqr(abs("
        "t_all.DATA[isnan(t_all.DATA)] - "
        "(t_1.DATA[isnan(t_1.DATA)] + t_2.DATA[isnan(t_2.DATA)] + t_3.DATA[isnan(t_3.DATA)])"
        ")) as diff from BDADDECal_onlypredict.MS t_all,"
        "PREDICT_DIR_1.MS t_1, PREDICT_DIR_2.MS t_2, PREDICT_DIR_3.MS t_3 )"
        "where diff > 1.e-6"
    )
    assert_taql(taql_check_visibilities)

    # ddecal should output the same weights it received
    taql_check_weights = (
        "select from"
        "(select gsumsqr(abs("
        "t_out.WEIGHT_SPECTRUM[isnan(t_out.WEIGHT_SPECTRUM)] - "
        "t_in.WEIGHT_SPECTRUM[isnan(t_in.WEIGHT_SPECTRUM)]"
        f")) as diff from BDADDECal_onlypredict.MS t_out, {MSIN} t_in)"
        "where diff > 1.e-6"
    )
    assert_taql(taql_check_weights)

# Only test a limited set of caltype + nchannels combinations, since testing
# all combinations does not have much extra value.
# The beginning # of the caltype_nchan string is the caltype.
# The number at the end is the number of channels.
@pytest.mark.parametrize(
    "caltype_nchan",
    [
        "scalar0",
        "scalarphase3",
        "scalaramplitude1",
        "diagonal2",
        "diagonalphase3",
        "diagonalamplitude2",
        "tec0",
        "tec3",
        "tecandphase1",
        # "tecscreen", # Requires Armadillo
        # "fulljones", # not implemented for BDA
        # "rotation", # part of fulljones -> not implemented
        # "rotation+diagonal", # part of fulljones -> not implemented
    ],
)

def test_caltype(create_skymodel, create_corrupted_data, caltype_nchan):
    """Test calibration for different calibration types"""
    caltype = caltype_nchan[:-1]
    nchan = int(caltype_nchan[-1])

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin=corrupted.MS",
            "msout=out.MS",
            "steps=[ddecal]",
            "ddecal.directions=[[center], [ra_off], [radec_off]]",
            "ddecal.h5parm=solutions.h5",
            "ddecal.sourcedb=test.skymodel",
            f"ddecal.mode={caltype}",
            "ddecal.solint=2",
            f"ddecal.nchan={nchan}",
            "numthreads=1",
        ]
    )

    import h5py # Don't import h5py when pytest is only collecting tests.
    h5 = h5py.File("solutions.h5", "r")

    if caltype in ["scalar", "diagonal", "scalaramplitude", "diagonalamplitude" ]:
        amplitude_solutions = h5['sol000/amplitude000/val']

        if caltype.startswith('scalar'):
            assert amplitude_solutions.attrs['AXES'] == b'time,freq,ant,dir'
        else:
            assert amplitude_solutions.attrs['AXES'] == b'time,freq,ant,dir,pol'

        if nchan == 0:
            assert amplitude_solutions.shape[1] == 1
        else:
            assert amplitude_solutions.shape[1] == 8 // nchan

        for corruption_num in range(3):
            np.testing.assert_array_almost_equal(
                amplitude_solutions[:,:,:,corruption_num],
                CORRUPTIONS[corruption_num],
                decimal=3
            )

    if caltype in ["scalar", "diagonal", "scalarphase", "diagonalphase"]:
        np.testing.assert_array_almost_equal(
            h5['sol000/phase000/val'][:,:,:,:],
            0.,
            decimal=5
        )

    if caltype in ["tec", "tecandphase"]:
        tec = h5['sol000/tec000/val']
        assert tec.attrs['AXES'] == b'time,ant,dir,freq'

    if caltype in ["tecandphase"]:
        phase = h5['sol000/phase000/val']
        assert phase.attrs['AXES'] == b'time,ant,dir,freq'

def test_subtract(create_skymodel, create_corrupted_data):
    """Test subtraction"""
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin=corrupted.MS",
            "msout=out.MS",
            "steps=[ddecal]",
            "ddecal.directions=[[center], [ra_off], [radec_off]]",
            "ddecal.h5parm=solutions.h5",
            "ddecal.sourcedb=test.skymodel",
            "ddecal.mode=diagonal",
            "ddecal.solint=2",
            "ddecal.nchan=3",
            "ddecal.subtract=true",
            "numthreads=1",
        ]
    )

    residual = float(check_output(
        [
            tcf.TAQLEXE,
            "-nopr",
            "-noph",
            "select gmax(abs(DATA)) from out.MS"
        ]
    ))

    assert residual < 0.01

    # ddecal should output the same weights it received
    taql_check_weights = (
        "select from"
        "(select gsumsqr(abs("
        "corrupted.WEIGHT_SPECTRUM[isnan(corrupted.WEIGHT_SPECTRUM)] - "
        "out.WEIGHT_SPECTRUM[isnan(out.WEIGHT_SPECTRUM)]"
        ")) as diff from corrupted.MS corrupted, out.MS out)"
        " where diff > 0"
    )
    assert_taql(taql_check_weights)
