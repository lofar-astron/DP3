# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
from subprocess import STDOUT, CalledProcessError, check_call, check_output

import numpy
import pytest
from numpy.random import PCG64, Generator

""" Append current directory to system path in order to import testconfig """
sys.path.append(".")
import testconfig
from utils import run_in_tmp_path

MSIN = f"{testconfig.BINDIR}/test_data/MWA-single-timeslot.ms"
SKYMODEL = f"{testconfig.BINDIR}/test_data/sky.txt"
H5FILE = f"{testconfig.BINDIR}/test_data/solutions.h5"
selected_antenna = 5


@pytest.fixture(autouse=True)
def skymodel(run_in_tmp_path):
    with open(SKYMODEL, "w") as sky_model_file:
        print(
            "Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='150000000', MajorAxis, MinorAxis, Orientation\n"
            "source-0,POINT,08h19m59.999s,-42d45m00s,1.0,[],false,150000000,,,",
            file=sky_model_file,
        )


def calldp3(parameters):
    print("DP3 " + " ".join(parameters))
    check_call([testconfig.DP3EXE] + parameters)


def run_with_diagonal_mode(diagonal_mode, true_solutions):
    # Perform the solve
    calldp3(
        [
            "checkparset=1",
            f"msin=test.ms",
            "msout=",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.h5parm={H5FILE}",
            "ddecal.mode=faradayrotation",
            "ddecal.maxiter=250",
            "ddecal.tolerance=1e-6",
            f"ddecal.faradaydiagonalmode={diagonal_mode}",
        ]
    )

    import h5py  # Don't import h5py when pytest is only collecting tests.

    with h5py.File(H5FILE, "r+") as h5file:
        rm_solutions = h5file["sol000/rotationmeasure000/val"]

        reference_antenna = 0
        selected_value = (
            rm_solutions[0, selected_antenna]
            - rm_solutions[0, reference_antenna]
        )
        true_value = (
            true_solutions[0, selected_antenna]
            - true_solutions[0, reference_antenna]
        )
        assert selected_value == pytest.approx(true_value, abs=1.0e-3)
        weights = h5file["sol000/rotationmeasure000/weight"]

        if diagonal_mode not in ["diagonalphase", "scalarphase", "rotation"]:
            amplitude_solutions = h5file["sol000/amplitude000/val"]
            assert amplitude_solutions[:, reference_antenna] == pytest.approx(
                1.0, abs=1.0e-3
            )
            assert amplitude_solutions[:, selected_antenna] == pytest.approx(
                1.0, abs=1.0e-3
            )


def test_faraday_constraint():

    # A filter step is added because the MS has fully flagged antennas, which
    # causes problems. The filter step removes that antenna.
    calldp3(
        [
            "checkparset=1",
            f"msin={MSIN}",
            "msout=test.ms",
            "msout.overwrite=true",
            "steps=[filter]",
            "filter.remove=True",
        ]
    )
    taqlcommand_run = (
        f"update test.ms set DATA=1, FLAG=False, WEIGHT_SPECTRUM=1"
    )
    check_output([testconfig.TAQLEXE, taqlcommand_run])
    # Generate a template H5.
    calldp3(
        [
            "checkparset=1",
            "msin=test.ms",
            "msout=",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.h5parm={H5FILE}",
            "ddecal.mode=faradayrotation",
        ]
    )

    # Modify h5 file so test the rotationmeasure
    import h5py  # Don't import h5py when pytest is only collecting tests.

    with h5py.File(H5FILE, "r+") as h5file:
        solutions = h5file["sol000/rotationmeasure000/val"]
        # Select a specific seed to have some reproducability
        rng = numpy.random.default_rng(0)
        true_solutions = rng.uniform(0.0, 1.0, solutions[:].shape)
        solutions[:] = true_solutions
        weights = h5file["sol000/rotationmeasure000/weight"]
        weights[:] = 1.0

    # Predict data with the given solutions
    calldp3(
        [
            "checkparset=1",
            f"msin=test.ms",
            "msout=",
            "steps=[predict]",
            f"predict.sourcedb={SKYMODEL}",
            f"predict.applycal.parmdb={H5FILE}",
            "predict.applycal.correction=rotationmeasure000",
        ]
    )

    run_with_diagonal_mode("scalaramplitude", true_solutions)
    run_with_diagonal_mode("diagonalamplitude", true_solutions)
    run_with_diagonal_mode("scalar", true_solutions)
    run_with_diagonal_mode("diagonal", true_solutions)
    run_with_diagonal_mode("diagonalphase", true_solutions)
    run_with_diagonal_mode("scalarphase", true_solutions)
