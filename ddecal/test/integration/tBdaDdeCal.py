# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
from subprocess import CalledProcessError, check_call, check_output, run

import numpy as np
import pytest

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_dp3, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/ddecal/test/integration directory,
  using `pytest source/tBdaDdeCal.py` (extended with pytest options of your choice)
- using ctest, see DP3/ddecal/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-bda.MS"
MSIN_REGULAR = "tNDPPP-generic.MS"
CORRUPTIONS = 3, 4, 7  # Corruption gain factors per antenna


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar(f"{tcf.RESOURCEDIR}/{MSIN_REGULAR}.tgz")


@pytest.fixture()
def skymodel_filename():
    """Create a skymodel file for tests and return its filename."""
    filename = "test.skymodel"
    with open(filename, "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\r\n"
            "center, POINT, 16:38:28.205000, +63.44.34.314000, 1, , , , , \r\n"
            "ra_off, POINT, 16:58:28.205000, +63.44.34.314000, 1, , , , , \r\n"
            "radec_off, POINT, 16:38:28.205000, +65.44.34.314000, 1, , , , , \r\n"
        )
    return filename


@pytest.fixture()
def create_corrupted_data():
    with open("test_corrupted.txt", "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\r\n"
        )
        f.write(
            f"center, POINT, 16:38:28.205000, +63.44.34.314000, {CORRUPTIONS[0] * CORRUPTIONS[0]}, , , , , \r\n"
        )
        f.write(
            f"ra_off, POINT, 16:58:28.205000, +63.44.34.314000, {CORRUPTIONS[1] * CORRUPTIONS[1]}, , , , , \r\n"
        )
        f.write(
            f"radec_off, POINT, 16:38:28.205000, +65.44.34.314000, {CORRUPTIONS[2] * CORRUPTIONS[2]}, , , , , \r\n"
        )

    run_dp3(
        [
            f"msin={MSIN}",
            "msout=corrupted.MS",
            "steps=[predict]",
            "predict.sourcedb=test_corrupted.txt",
        ]
    )

    check_call([tcf.TAQLEXE, "update corrupted.MS set WEIGHT_SPECTRUM=1"])


@pytest.fixture()
def create_corrupted_data_from_regular():
    """
    Use a lower frequency averaging to ensure the condition on the frequency
    channels is satisfied: channels per chan block >= max averaged channels
    """

    with open("test_corrupted.txt", "w") as f:
        f.write(
            "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\r\n"
        )
        f.write(
            f"center, POINT, 16:38:28.205000, +63.44.34.314000, {CORRUPTIONS[0] * CORRUPTIONS[0]}, , , , , \r\n"
        )
        f.write(
            f"ra_off, POINT, 16:58:28.205000, +63.44.34.314000, {CORRUPTIONS[1] * CORRUPTIONS[1]}, , , , , \r\n"
        )
        f.write(
            f"radec_off, POINT, 16:38:28.205000, +65.44.34.314000, {CORRUPTIONS[2] * CORRUPTIONS[2]}, , , , , \r\n"
        )

    run_dp3(
        [
            f"msin={MSIN_REGULAR}",
            "msout=corrupted.MS",
            "steps=[bdaaverager, predict]",
            "bdaaverager.frequencybase=40",
            "bdaaverager.timebase=100",
            "predict.sourcedb=test_corrupted.txt",
        ]
    )

    check_call([tcf.TAQLEXE, "update corrupted.MS set WEIGHT_SPECTRUM=1"])


def test_only_predict(skymodel_filename):
    """Test that the right patches are summed with predict_only"""

    common_args = [
        f"msin={MSIN}",
        "msout.overwrite=true",
    ]

    predict_args = [
        "steps=[predict]",
        f"predict.sourcedb={skymodel_filename}",
    ]

    run_dp3(
        [
            "msout=BDADDECal_onlypredict.MS",
            "steps=[ddecal]",
            "ddecal.onlypredict=true",
            "ddecal.directions=[[center, dec_off],[ra_off],[radec_off]]",
            f"ddecal.sourcedb={skymodel_filename}",
        ]
        + common_args
    )

    run_dp3(
        [
            "msout=PREDICT_DIR_1.MS",
            "predict.sources=[center, dec_off]",
        ]
        + common_args
        + predict_args
    )

    run_dp3(
        ["msout=PREDICT_DIR_2.MS", "predict.sources=[ra_off]"]
        + common_args
        + predict_args
    )

    run_dp3(
        [
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


def test_uvwflagger(skymodel_filename, create_corrupted_data_from_regular):
    """Test that uvwflagger settings lead to the right amount of NaNs in the solution file"""

    run_dp3(
        [
            f"msin={MSIN}",
            "msout=out.MS",
            "steps=[ddecal]",
            "ddecal.directions=[[center], [ra_off], [radec_off]]",
            "ddecal.h5parm=solutions.h5",
            f"ddecal.sourcedb={skymodel_filename}",
            "ddecal.mode=scalar",
            "ddecal.solint=2",
            "ddecal.nchan=10",
            "ddecal.uvlambdamin=12000.0",
        ]
    )

    import h5py  # Don't import h5py when pytest is only collecting tests.

    h5 = h5py.File("solutions.h5", "r")
    amplitude_solutions = h5["sol000/amplitude000/val"]
    phase_solutions = h5["sol000/phase000/val"]

    # When uvw flagging is disabled, the NaNs in the solution file are only 9
    expected_flagged_solutions = 54

    assert (
        np.count_nonzero(np.isnan(amplitude_solutions))
        == expected_flagged_solutions
    )
    assert (
        np.count_nonzero(np.isnan(phase_solutions))
        == expected_flagged_solutions
    )


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
def test_caltype(
    skymodel_filename, create_corrupted_data_from_regular, caltype_nchan
):
    """Test calibration for different calibration types"""
    caltype = caltype_nchan[:-1]
    nchan = int(caltype_nchan[-1])
    print(nchan)

    run_dp3(
        [
            f"msin=corrupted.MS",
            "msout=out.MS",
            "steps=[ddecal]",
            "ddecal.h5parm=solutions.h5",
            f"ddecal.sourcedb={skymodel_filename}",
            f"ddecal.mode={caltype}",
            "ddecal.solint=2",
            f"ddecal.nchan={nchan}",
        ]
    )

    import h5py  # Don't import h5py when pytest is only collecting tests.

    h5 = h5py.File("solutions.h5", "r")

    if caltype in [
        "scalar",
        "diagonal",
        "scalaramplitude",
        "diagonalamplitude",
    ]:
        amplitude_solutions = h5["sol000/amplitude000/val"]

        if caltype.startswith("scalar"):
            assert amplitude_solutions.attrs["AXES"] == b"time,freq,ant,dir"
        else:
            assert (
                amplitude_solutions.attrs["AXES"] == b"time,freq,ant,dir,pol"
            )

        if nchan == 0:
            assert amplitude_solutions.shape[1] == 1
        else:
            assert amplitude_solutions.shape[1] == 8 // nchan

        for corruption_num in range(3):
            np.testing.assert_array_almost_equal(
                amplitude_solutions[:, :, :, corruption_num],
                CORRUPTIONS[corruption_num],
                decimal=3,
            )

    if caltype in ["scalar", "diagonal", "scalarphase", "diagonalphase"]:
        np.testing.assert_array_almost_equal(
            h5["sol000/phase000/val"][:, :, :, :], 0.0, decimal=5
        )

    if caltype in ["tec", "tecandphase"]:
        tec = h5["sol000/tec000/val"]
        assert tec.attrs["AXES"] == b"time,ant,dir,freq"

    if caltype in ["tecandphase"]:
        phase = h5["sol000/phase000/val"]
        assert phase.attrs["AXES"] == b"time,ant,dir,freq"


def test_subtract(skymodel_filename, create_corrupted_data):
    """Test subtraction"""
    run_dp3(
        [
            f"msin=corrupted.MS",
            "msout=out.MS",
            "steps=[ddecal]",
            # Use explicit directions in this test.
            "ddecal.directions=[[center], [ra_off], [radec_off]]",
            "ddecal.h5parm=solutions.h5",
            f"ddecal.sourcedb={skymodel_filename}",
            "ddecal.mode=diagonal",
            "ddecal.solint=2",
            "ddecal.nchan=8",
            "ddecal.subtract=true",
        ]
    )

    residual = float(
        check_output(
            [
                tcf.TAQLEXE,
                "-nopr",
                "-noph",
                "select gmax(abs(DATA)) from out.MS",
            ]
        )
    )

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


def test_invalid_input(skymodel_filename):
    """Assert that exception is thrown when an incompatible value of solint or nchan is given"""

    common_args = [
        f"msin={MSIN}",
        "msout=out.MS",
        "steps=[ddecal]",
        "ddecal.h5parm=solutions.h5",
        f"ddecal.sourcedb={skymodel_filename}",
        "ddecal.mode=diagonal",
        "ddecal.subtract=true",
    ]

    with pytest.raises(CalledProcessError):
        run_dp3(["ddecal.solint=1", "ddecal.nchan=4"] + common_args)

    with pytest.raises(CalledProcessError):
        run_dp3(["ddecal.solint=2", "ddecal.nchan=1"] + common_args)


def test_reuse_model_data(skymodel_filename):
    # Apply ddecal directly and generate reference output.
    run_dp3(
        [
            f"msin={MSIN}",
            "steps=[ddecal]",
            f"ddecal.sourcedb={skymodel_filename}",
            "ddecal.h5parm=cal_ref.h5",
            "ddecal.directions=[[center,dec_off,ra_off,radec_off]]",
            "ddecal.subtract=true",
            "ddecal.solint=2",
            "ddecal.nchan=10",
            "msout=ref.ms",
        ]
    )

    # Run ddecal twice where the second ddecal reuses model data.
    run_dp3(
        [
            f"msin={MSIN}",
            "steps=[ddecal1,ddecal2]",
            f"ddecal1.sourcedb={skymodel_filename}",
            "ddecal1.directions=[[center,dec_off,ra_off,radec_off]]",
            "ddecal1.onlypredict=true",
            "ddecal1.keepmodel=true",
            "ddecal2.reusemodel=[ddecal1.center]",
            "ddecal2.h5parm=cal_reuse.h5",
            "ddecal2.subtract=true",
            "ddecal2.solint=2",
            "ddecal2.nchan=10",
            "msout=reused.ms",
        ]
    )

    assert_taql(
        "select from (select abs(gsumsqr(reused.DATA - ref.DATA)) as diff "
        "from reused.ms reused, ref.ms ref) where diff>1.0e-6"
    )
