# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from subprocess import STDOUT, CalledProcessError, check_output

import numpy as np
import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import casacore.tables
import testconfig as tcf
from utils import (
    COMMON_DP3_ARGUMENTS,
    assert_taql,
    run_dp3,
    run_in_tmp_path,
    untar,
)

"""
Script can be invoked in two ways:
- as standalone from the build/ddecal/test/integration directory,
  using `pytest source/tDDECal.py` (extended with pytest options of your choice)
- using ctest, see DP3/ddecal/test/integration/CMakeLists.txt
"""

IDG_RESOURCES = "idg-fits-sources.tbz2"
MSINTGZ = "tDDECal.in_MS.tgz"
MSIN = "tDDECal.MS"
SKYMODEL = f"{MSIN}/sky.txt"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSINTGZ}")
    # The test table has an invalid keyword, which needs to be removed
    # to make Dp3 behave properly.
    table = casacore.tables.table(MSIN, readonly=False)
    table.removekeyword("sky")
    table.close()


@pytest.fixture()
def idgpredict_env():
    """Extract data for testing DDECal with IDGPredict"""
    untar(f"{tcf.DDECAL_RESOURCEDIR}/{IDG_RESOURCES}")


@pytest.fixture()
def copy_data_to_model_data():
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[]",
        ]
    )


@pytest.fixture()
def create_corrupted_visibilities():
    taqlcommand = f"update {MSIN} set WEIGHT_SPECTRUM=1, FLAG=False"
    check_output([tcf.TAQLEXE, "-noph", taqlcommand])

    # Use ddecal to create template h5parm
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
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
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DATA",
            "steps=[h5parmpredict]",
            f"h5parmpredict.sourcedb={SKYMODEL}",
            "h5parmpredict.applycal.parmdb=instrumentcorrupted.h5",
            "h5parmpredict.applycal.correction=amplitude000",
        ]
    )


def test_modeldata_different_shape(copy_data_to_model_data):
    base_command = [
        tcf.DP3EXE,
        f"msin={MSIN}",
        "ddecal.modeldatacolumns=[MODEL_DATA]",
        "msout=.",
        "steps=[filter,ddecal]",
        "filter.startchan=3",
    ] + COMMON_DP3_ARGUMENTS

    with pytest.raises(CalledProcessError) as e:
        check_output(base_command, stderr=STDOUT)

        error_string = "The column MODEL_DATA has shape NCORR: 4, NCHAN: 8, "
        "while the input buffer for the current step has shape NCORR: 4, "
        "NCHAN: 5.\nAny operation which alters the number of channels, "
        "correlations and baselines can cause a shape mismatch (example: "
        "filter, average, ...).\n\n"
        assert error_string in e.output.decode()


@pytest.mark.parametrize(
    "caltype",
    ["complexgain", "scalarcomplexgain", "amplitudeonly", "scalaramplitude"],
)
@pytest.mark.parametrize("solint", [0, 1, 2, 4])
@pytest.mark.parametrize("nchan", [1, 2, 5])
def test(
    create_corrupted_visibilities,
    copy_data_to_model_data,
    caltype,
    solint,
    nchan,
):
    # Subtract corrupted visibilities using multiple predict steps
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.solint={solint}",
            f"ddecal.nchan={nchan}",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
        ]
    )

    # Calibrate on the original sources, caltype=$caltype
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3]",
            f"predict1.sourcedb={SKYMODEL}",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center,dec_off]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={SKYMODEL}",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[radec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={SKYMODEL}",
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
    check_output([tcf.TAQLEXE, "-noph", taqlcommand_run])
    taql_command = f"select FROM (select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*DATA))) as norm_data, sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*SUBTRACTED_DATA))) as norm_residual from {MSIN}) where norm_residual/norm_data > {tolerance} or isinf(norm_residual/norm_data) or isnan(norm_residual/norm_data)"
    assert_taql(taql_command)


@pytest.mark.parametrize(
    "solutions_per_direction", [[1, 1, 1], [1, 3, 6], [3, 3, 3], [6, 6, 6]]
)
@pytest.mark.parametrize(
    "caltype",
    ["scalaramplitude", "scalar"],
)
def test_calibration_with_dd_intervals(
    create_corrupted_visibilities,
    copy_data_to_model_data,
    solutions_per_direction,
    caltype,
):
    # This test checks that the calibration with different solution intervals per direction gives the same result as the corruption applied in the fixture "create_corrupted_visibilities".

    # Run calibration on corrupted dataset and store solutions in instrument.h5
    sol_int = 6
    n_timeslots_in_ms = 6

    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.solint={sol_int}",
            "ddecal.nchan=1",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
            "ddecal.solveralgorithm=directioniterative",
            f"ddecal.solutions_per_direction={solutions_per_direction}",
        ]
    )

    # Verify that he values in instrument.h5 are close to the values in instrumentcorrupted.h5
    import h5py  # Don't import h5py when pytest is only collecting tests.

    ddecal_solutions = h5py.File("instrument.h5", "r")
    reference_solutions = h5py.File("instrumentcorrupted.h5", "r")

    assert (
        ddecal_solutions["sol000/amplitude000/val"].attrs["AXES"]
        == b"time,freq,ant,dir"
    )

    for direction_index in range(
        0, len(solutions_per_direction)
    ):  # loop over number of directions
        for solint_index in range(
            0, int(np.ceil(n_timeslots_in_ms / sol_int))
        ):  # loop over number of solution intervals
            for interval_in_direction_index in range(
                0, solutions_per_direction[direction_index]
            ):  # loop over dd-intervals within one solution interval
                values_ddecal = ddecal_solutions["sol000/amplitude000/val"][
                    solint_index * np.max(solutions_per_direction)
                    + interval_in_direction_index,
                    :,
                    :,
                    direction_index,
                ]
                corresponding_index = solint_index * np.max(
                    solutions_per_direction
                ) + int(sol_int / solutions_per_direction[direction_index])

                if (
                    corresponding_index
                    >= reference_solutions["sol000/amplitude000/val"].shape[0]
                ):
                    corresponding_index = (
                        reference_solutions["sol000/amplitude000/val"].shape[0]
                        - 1
                    )

                values_reference = reference_solutions[
                    "sol000/amplitude000/val"
                ][corresponding_index, :, :, direction_index, 0]

                assert (abs(values_ddecal - values_reference) < 1.2).all()


@pytest.mark.xfail(reason="Will be fixed in AST-924")
@pytest.mark.parametrize(
    "solutions_per_direction",
    [[1, 1, 1], [1, 2, 4], [2, 2, 2], [4, 4, 4]],
)
@pytest.mark.parametrize(
    "caltype",
    ["scalaramplitude", "scalar", "diagonal", "diagonalamplitude"],
)
def test_bug_ast_924(
    create_corrupted_visibilities,
    copy_data_to_model_data,
    solutions_per_direction,
    caltype,
):
    # This is the same test as "test_calibration_with_dd_intervals" , but with a different solution interval:
    # In this case sol_int is not an integer divisior of n_timeslots_in_ms

    sol_int = 4
    n_timeslots_in_ms = 6
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.solint={sol_int}",
            "ddecal.nchan=1",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
            "ddecal.solveralgorithm=directioniterative",
            f"ddecal.solutions_per_direction={solutions_per_direction}",
        ]
    )

    # Verify that he values in instrument.h5 are close to the values in instrumentcorrupted.h5
    import h5py  # Don't import h5py when pytest is only collecting tests.

    ddecal_solutions = h5py.File("instrument.h5", "r")
    reference_solutions = h5py.File("instrumentcorrupted.h5", "r")

    for direction_index in range(
        0, len(solutions_per_direction)
    ):  # loop over number of directions
        for solint_index in range(
            0, int(np.ceil(n_timeslots_in_ms / sol_int))
        ):  # loop over number of solution intervals
            for interval_in_direction_index in range(
                0, solutions_per_direction[direction_index]
            ):  # loop over dd-intervals within one solution interval
                values_ddecal = ddecal_solutions["sol000/amplitude000/val"][
                    solint_index * np.max(solutions_per_direction)
                    + interval_in_direction_index,
                    :,
                    :,
                    direction_index,
                ]
                corresponding_index = solint_index * np.max(
                    solutions_per_direction
                ) + int(sol_int / solutions_per_direction[direction_index])

                if (
                    corresponding_index
                    >= reference_solutions["sol000/amplitude000/val"].shape[0]
                ):
                    corresponding_index = (
                        reference_solutions["sol000/amplitude000/val"].shape[0]
                        - 1
                    )

                values_reference = reference_solutions[
                    "sol000/amplitude000/val"
                ][corresponding_index, :, :, direction_index, 0]

                assert (abs(values_ddecal - values_reference) < 1.2).all()


@pytest.mark.xfail(reason="Will be implemented in AST-920")
@pytest.mark.parametrize(
    "solutions_per_direction", [[1, 1, 1], [1, 3, 6], [3, 3, 3], [6, 6, 6]]
)
@pytest.mark.parametrize(
    "caltype",
    [
        "amplitudeonly",
        "scalaramplitude",
        "scalar",
        "diagonal",
        "diagonalamplitude",
    ],
)
def test_subtract_with_dd_intervals(
    create_corrupted_visibilities,
    copy_data_to_model_data,
    solutions_per_direction,
    caltype,
):
    # This test checks that the subtraction operation works correctly
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DURING_DDECAL",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.solint=6",
            "ddecal.nchan=2",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
            "ddecal.solveralgorithm=directioniterative",
            f"ddecal.solutions_per_direction={solutions_per_direction}",
            "ddecal.subtract=true",
        ]
    )

    common_predict_command = [
        f"msin={MSIN}",
        "msout=.",
        "msout.datacolumn=SUBTRACTED_DURING_PREDICT",
        "steps=[predict1,predict2,predict3]",
        f"predict1.sourcedb={SKYMODEL}",
        "predict1.applycal.parmdb=instrument.h5",
        "predict1.sources=[center,dec_off]",
        "predict1.operation=subtract",
        f"predict2.sourcedb={SKYMODEL}",
        "predict2.applycal.parmdb=instrument.h5",
        "predict2.sources=[radec_off]",
        "predict2.operation=subtract",
        f"predict3.sourcedb={SKYMODEL}",
        "predict3.applycal.parmdb=instrument.h5",
        "predict3.sources=[ra_off]",
        "predict3.operation=subtract",
    ]
    corrections_amplitude_only = [
        "predict1.applycal.correction=amplitude000",
        "predict2.applycal.correction=amplitude000",
        "predict3.applycal.correction=amplitude000",
    ]

    corrections_amplitude_and_phase = [
        "predict1.applycal.steps=[amplitude,phase]",
        "predict1.applycal.amplitude.correction=amplitude000",
        "predict1.applycal.phase.correction=phase000",
        "predict2.applycal.steps=[amplitude,phase]",
        "predict2.applycal.amplitude.correction=amplitude000",
        "predict2.applycal.phase.correction=phase000",
        "predict3.applycal.steps=[amplitude,phase]",
        "predict3.applycal.amplitude.correction=amplitude000",
        "predict3.applycal.phase.correction=phase000",
    ]

    if caltype == "scalar" or caltype == "diagonal":
        run_dp3(common_predict_command + corrections_amplitude_and_phase)

    else:
        run_dp3(common_predict_command + corrections_amplitude_only)

    # Quantify the difference between the subtraction in the DDECal step and in the Predict step
    predict_residual = float(
        check_output(
            [
                tcf.TAQLEXE,
                "-nopr",
                "-noph",
                f"select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*(SUBTRACTED_DURING_DDECAL- SUBTRACTED_DURING_PREDICT)))) from {MSIN}",
            ]
        )
    )

    tolerance = 2
    assert predict_residual < tolerance


def test_h5parm_predict():
    # make calibration solutions
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            # By omitting the direction, the directions are determined from the
            # sourcedb. This gives the same results as:
            # "ddecal.directions=[[center],[dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode=diagonal",
        ]
    )

    # subtract using multiple predict steps
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3,predict4]",
            f"predict1.sourcedb={SKYMODEL}",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={SKYMODEL}",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[dec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={SKYMODEL}",
            "predict3.applycal.parmdb=instrument.h5",
            "predict3.sources=[radec_off]",
            "predict3.operation=subtract",
            "predict3.applycal.correction=amplitude000",
            f"predict4.sourcedb={SKYMODEL}",
            "predict4.applycal.parmdb=instrument.h5",
            "predict4.sources=[ra_off]",
            "predict4.operation=subtract",
            "predict4.applycal.correction=amplitude000",
        ]
    )

    # subtract using h5parmpredict
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA_H5PARM",
            "steps=[h5parmpredict]",
            f"h5parmpredict.sourcedb={SKYMODEL}",
            "h5parmpredict.applycal.parmdb=instrument.h5",
            "h5parmpredict.operation=subtract",
            "h5parmpredict.applycal.correction=amplitude000",
        ]
    )

    # echo "Check that h5parmpredict creates the same output as multiple predict steps"
    taql_command = f"select from (select abs(gsumsqr(SUBTRACTED_DATA-SUBTRACTED_DATA_H5PARM)) as diff from {MSIN}) where diff>1.e-6"
    assert_taql(taql_command)


def test_oneapplycal_from_buffer():

    # Apply ddecal and oneapplycal sequentially
    run_dp3(
        [
            f"msin={MSIN}",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.h5parm=cal.h5",
            "ddecal.directions=[[center,dec_off,ra_off,radec_off]]",
            "msout=.",
            "msout.datacolumn=DATA_REF_BUF",
        ]
    )
    run_dp3(
        [
            f"msin={MSIN}",
            "msin.datacolumn=DATA_REF_BUF",
            "steps=[applycal]",
            "applycal.parmdb=cal.h5",
            "applycal.steps=[phase,ampl]",
            "applycal.phase.correction=phase000",
            "applycal.ampl.correction=amplitude000",
            "msout=.",
            "msout.datacolumn=DATA_REF_BUF",
        ]
    )

    # DDECal and immediate application
    run_dp3(
        [
            f"msin={MSIN}",
            "steps=[ddecal,applycal]",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.directions=[[center,dec_off,ra_off,radec_off]]",
            "ddecal.storebuffer=True",
            "ddecal.h5parm=unused_output.h5",
            "msout=.",
            "msout.datacolumn=DATA_NEW_BUF",
            "applycal.parmdb=",
        ]
    )

    # echo "Check that h5parmpredict creates the same output as multiple predict steps"
    taql_command = f"select from (select abs(gsumsqr(DATA_REF_BUF-DATA_NEW_BUF)) as diff from {MSIN}) where diff>1.e-6"
    assert_taql(taql_command)


def test_pre_apply():
    # make calibration solutions
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode=diagonal",
        ]
    )

    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3]",
            f"predict1.sourcedb={SKYMODEL}",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center,dec_off]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={SKYMODEL}",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[radec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={SKYMODEL}",
            "predict3.applycal.parmdb=instrument.h5",
            "predict3.sources=[ra_off]",
            "predict3.operation=subtract",
            "predict3.applycal.correction=amplitude000",
        ]
    )

    # Check that preapply runs (output is not tested)
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.applycal.parmdb=instrument.h5",
            "ddecal.applycal.steps=applyampl",
            "ddecal.applycal.applyampl.correction=amplitude000",
            "ddecal.h5parm=instrument2.h5",
            "ddecal.mode=scalarcomplexgain",
        ]
    )


def test_check_tec():
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument-tec.h5",
            "ddecal.mode=tec",
        ]
    )


def test_check_tec_and_phase():
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msin.baseline='!CS001HBA0'",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
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
        f"msin={MSIN}",
        "msout=.",
        "steps=[ddecal]",
        f"ddecal.sourcedb={SKYMODEL}",
        "ddecal.solint=6",
        "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
        "ddecal.h5parm=instrument-tec.h5",
        "ddecal.mode=scalaramplitude",
        "ddecal.solveralgorithm=directioniterative",
    ] + COMMON_DP3_ARGUMENTS
    try:
        check_output(
            (
                base_command
                if solutions_per_direction is None
                else base_command
                + [f"ddecal.solutions_per_direction={solutions_per_direction}"]
            ),
            stderr=STDOUT,
        )
    except CalledProcessError as e:
        # AST-1184: The output may contain other warnings, e.g., libgcov version
        # mismatch errors. Therefore check if the output contains the exception.
        if solutions_per_direction == [5, 5, 5]:
            if (
                not "\nstd exception detected: Values in solutions_per_direction should be integer divisors of solint"
                in e.output.decode()
            ):
                raise e
        elif solutions_per_direction == [2, 0]:
            if (
                not "\nstd exception detected: All entries in solutions_per_direction should be > 0.\n"
                in e.output.decode()
            ):
                raise e
        else:
            raise (e)


def test_modelnextsteps(copy_data_to_model_data):
    import h5py  # Don't import h5py when pytest is only collecting tests.

    # Multiply MODEL_DATA by 42
    taqlcommand_run = f"update {MSIN} set MODEL_DATA=DATA*42"
    check_output([tcf.TAQLEXE, "-noph", taqlcommand_run])

    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[ddecal]",
            "ddecal.modeldatacolumns=[MODEL_DATA]",
            "ddecal.modelnextsteps.MODEL_DATA=[scaledata]",
            "ddecal.h5parm=instrument-modeldata.h5",
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

    with h5py.File("instrument-modeldata.h5", "r") as h5file:
        sol = h5file["sol000/amplitude000/val"]
        # First check one element to get a nice error message
        assert sol[0, 0, 0, 0, 0] == pytest.approx(
            1.0 / np.sqrt(42), abs=1.0e-3
        )
        assert np.all(
            np.isclose(sol[np.isfinite(sol)], 1 / np.sqrt(42), atol=1.0e-3)
        )


def test_bda_constraints():
    import h5py  # Don't import h5py when pytest is only collecting tests.

    common_test_arguments = [
        f"msin={MSIN}",
        "msout.overwrite=true",
        "msin.datacolumn=DATA",
        "solve.type=ddecal",
        "solve.mode=complexgain",
        "solve.usebeammodel=True",
        "solve.beammode=array_factor",
        "solve.maxiter=150",
        "solve.directions=[[center,dec_off],[ra_off],[radec_off]]",
        "numthreads=6",
        "solve.onebeamperpatch=False",
        "solve.propagatesolutions=True",
        "solve.smoothnessconstraint=4000000.0",
        "solve.solint=37",
        "solve.nchan=10",
        "solve.solveralgorithm=hybrid",
        f"solve.sourcedb={SKYMODEL}",
        "solve.stepsize=0.2",
        "solve.tolerance=0.005",
    ]

    run_dp3(
        [
            "msout=tmp_bda.ms",
            "steps=[avg,solve]",
            "avg.type=bdaaverager",
            "avg.timebase=20000",
            "avg.frequencybase=10000",
            "avg.maxinterval=37",
            "solve.h5parm=test_bda.h5parm",
        ]
        + common_test_arguments
    )

    run_dp3(
        [
            "msout=tmp_no_bda.ms",
            "steps=[solve]",
            "solve.h5parm=test_no_bda.h5parm",
        ]
        + common_test_arguments
    )

    f_bda = h5py.File("test_bda.h5parm", "r")
    f_no_bda = h5py.File("test_no_bda.h5parm", "r")

    ampl_bda = f_bda["sol000/amplitude000/val"]
    ampl_no_bda = f_no_bda["sol000/amplitude000/val"]
    phase_bda = f_bda["sol000/phase000/val"]
    phase_no_bda = f_no_bda["sol000/phase000/val"]

    np.testing.assert_allclose(
        ampl_bda, ampl_no_bda, rtol=0.05, atol=0, equal_nan=True
    )
    np.testing.assert_allclose(
        phase_bda, phase_no_bda, rtol=0.3, atol=0, equal_nan=True
    )


@pytest.mark.parametrize(
    "caltype",
    [
        "scalar",
        "diagonal",
        "fulljones",
    ],
)
def test_station_with_auto_correlation_only(caltype):
    """
    AST-1146: Test passing input to DDECal that contains a station that
    only has visibilities for its auto-correlation, and no other visibilities.
    Since DDECal skips auto-correlations, it ends up with a station without
    visibilities, which can cause issues like empty matrices / null pointers.
    """

    # Add a station with an auto-correlation to the MS.
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=StationAdded.MS",
            "msout.uvwcompression=false",
            "steps=[stationadder]",
            "stationadder.stations={AUTO:[RS106HBA,RS208HBA]}",
            "stationadder.autocorr=true",
            "stationadder.sumauto=false",
        ]
    )

    # Filter out the non-auto-correlations for the added station and run DDECal.
    # Since the Filter step uses the input MS when parsing the filter, we cannot
    # combine the DP3 call above with this call.
    run_dp3(
        [
            f"msin=StationAdded.MS",
            "msout=Test.MS",
            "msout.uvwcompression=false",
            "steps=[filter,ddecal]",
            # This filter keeps the auto-correlations for station 'AUTO'
            # and removes all other visibilities for 'AUTO'.
            "filter.baseline=!AUTO",
            f"ddecal.mode={caltype}",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.solint=2",
            "ddecal.nchan=0",
            "ddecal.h5parm=solutions.h5",
        ]
    )


def test_uvwflagging():
    """
    Test whether DDECal
    1) Does not use the visibilities in the solver that are excluded by a UVW selection
    2) Does not propagate the flags used for UVW selection to the next step

    This test is nearly identical to test_uvwflagging() in steps/test/integration/tGainCal.py
    """

    # Clear flags
    taqlcommand_run = f"UPDATE {MSIN} SET FLAG=FALSE"
    check_output([tcf.TAQLEXE, "-noph", taqlcommand_run])

    # UVW Flagging
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=",
            "steps=[uvwflagger]",
            "uvwflagger.uvmmax=10000",
        ]
    )

    # Verify that the test scenario is sane and does indeed contain flagged UVW values
    # Check that 84 rows have been flagged
    assert_taql(f"SELECT FLAG FROM {MSIN} WHERE ANY(FLAG)", 84)

    # Fill DATA column diagonal with valid visibilities of value 4.0, for unflagged data,
    # or corrupted visibilities with value 8.0, for flagged data
    taqlcommand_run = f"UPDATE {MSIN} SET DATA=iif(FLAG, 8.0, 4.0)*RESIZE([1,0,0,1],SHAPE(DATA),1)"
    check_output([tcf.TAQLEXE, "-noph", taqlcommand_run])

    # Clear flags
    taqlcommand_run = f"UPDATE {MSIN} SET FLAG=FALSE"
    check_output([tcf.TAQLEXE, "-noph", taqlcommand_run])

    # Create a sky model
    import casacore.tables  # Don't import casacore.tables when pytest is only collecting tests.

    t_field = casacore.tables.table(f"{MSIN}::FIELD")
    ra, dec = t_field[0]["PHASE_DIR"][0]
    with open("skymodel.txt", "w") as sky_model_file:
        print(
            "Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='150000000', MajorAxis, MinorAxis, Orientation",
            file=sky_model_file,
        )
        print(
            f"source-0,POINT,{ra},{dec},1.0,[],false,150000000,,,",
            file=sky_model_file,
        )

    # Calibrate
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal1,ddecal2]",
            "ddecal1.sourcedb=skymodel.txt",
            "ddecal1.solint=0",
            "ddecal1.nchan=0",
            "ddecal1.h5parm=instrument1.h5",
            "ddecal1.mode=complexgain",
            "ddecal1.uvmmax=10000",  # When this line is removed the corrupted visibilities will be used
            # and the test will fail
            "ddecal2.sourcedb=skymodel.txt",
            "ddecal2.solint=0",
            "ddecal2.nchan=0",
            "ddecal2.h5parm=instrument2.h5",
            "ddecal2.mode=complexgain",
        ]
    )

    # Check solutions
    import h5py  # Don't import h5py when pytest is only collecting tests.

    with h5py.File("instrument1.h5", "r") as h5file:
        sol = h5file["sol000/amplitude000/val"]
        # Because of the flags some antennas will have no solution, indicated by a NaN.
        # Check whether the finite solutions are as expected.
        # The equation that is solved for is:
        #   g_i * model_visibility * g_j = visibility.
        # The visibilities (DATA column) are 4.0, the model visibilities are 1.0
        # so the gain solutions should be 2.0.
        # First check that at least 70% of the solutions are valid
        assert np.sum(np.isfinite(sol)) / sol.size > 0.7
        assert np.all(np.isclose(sol[np.isfinite(sol)], 2.0, atol=1.0e-3))

    with h5py.File("instrument2.h5", "r") as h5file:
        sol = h5file["sol000/amplitude000/val"]
        # The second step has no uvwflaggin, so the result should be incorrect
        # First check that at least 70% of the solutions are valid
        assert np.sum(np.isfinite(sol)) / sol.size > 0.7
        assert not np.all(np.isclose(sol[np.isfinite(sol)], 2.0, atol=1.0e-3))

    # Check flags
    # The flags used internally by DDECal to flag unwanted uvw values should not propagate
    # to the next step/msout step
    assert_taql(f"SELECT FLAG FROM {MSIN} WHERE ANY(FLAG)")


def test_minvisratio(copy_data_to_model_data):
    """
    Test whether DDECal
    1) Applies the "minvisratio" setting correctly.
    2) Does not modify flags and weights in that case.
    """

    # Set data values, flags and weights. For the second channel,
    # set 25% of the flags by flagging the second correlation.
    taqlcommand_run = f"UPDATE {MSIN} SET DATA=42, MODEL_DATA=42, WEIGHT_SPECTRUM=4, FLAG=FALSE, FLAG[1,1]=TRUE"
    check_output([tcf.TAQLEXE, "-noph", taqlcommand_run])

    # Calibrate
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            # Force writing flags and weights.
            # Since DDECal does not have them in its provided fields, the output
            # step will not write them by default. Force-writing them allows
            # checking if DDECal indeed does not change the flags and weights.
            "msout.flagcolumn=DDEFLAG",
            "msout.weightcolumn=DDEWEIGHTS",
            "steps=[ddecal]",
            "ddecal.modeldatacolumns=MODEL_DATA",
            "ddecal.subtract=true",
            # With a minvisratio of 80%, DDECal should flag the second channel,
            # since only 75% is unflagged.
            "ddecal.minvisratio=0.8",
            "ddecal.h5parm=solutions.h5",
        ]
    )

    # Check output data, flags, and weights.
    # The data should be NaN for the second channel and zero for the other channels.
    # The model data, flags and weights should be equal to the input values.
    assert_taql(f"SELECT FROM {MSIN} WHERE ANY(!ISNAN(DATA[1,]))")
    assert_taql(f"SELECT FROM {MSIN} WHERE ANY(!NEARABS(DATA[0,], 0, 1e-5))")
    assert_taql(f"SELECT FROM {MSIN} WHERE ANY(!NEARABS(DATA[2:,], 0, 1e-5))")
    assert_taql(f"SELECT FROM {MSIN} WHERE ANY(MODEL_DATA!=42)")
    # Since MsReader flags all correlations if a single correlation is flagged,
    # the output flags do not match the input flags in this test.
    # TODO(AST-1280): Convert this test into a C++ unit test, once DDECal can
    # read model visibilities from a DPBuffer.
    assert_taql(
        f"SELECT FROM {MSIN} WHERE NTRUE(DDEFLAG)!=4 OR ANY(DDEFLAG[1,]=FALSE)"
    )
    assert_taql(f"SELECT FROM {MSIN} WHERE ANY(DDEWEIGHTS!=4)")


def test_extra_data_columns():
    """
    Test whether DDECal can reuse model data from extra data columns.
    """
    run_dp3(
        [
            f"msin={MSIN}",
            "msin.extradatacolumns=[DATA]",
            "steps=[ddecal]",
            "ddecal.reusemodel=[DATA]",
            "ddecal.h5parm=solutions.h5",
            "msout=.",
        ]
    )


def test_reuse_model_data():
    # Apply ddecal directly and generate reference output.
    run_dp3(
        [
            f"msin={MSIN}",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.h5parm=cal_ref.h5",
            "ddecal.directions=[[center,dec_off,ra_off,radec_off]]",
            "ddecal.subtract=true",
            "msout=.",
            "msout.datacolumn=DATA_REF",
        ]
    )

    # Run ddecal twice where the second ddecal reuses model data.
    run_dp3(
        [
            f"msin={MSIN}",
            "steps=[ddecal1,ddecal2]",
            f"ddecal1.sourcedb={SKYMODEL}",
            "ddecal1.directions=[[center,dec_off,ra_off,radec_off]]",
            "ddecal1.onlypredict=true",
            "ddecal1.keepmodel=true",
            "ddecal2.reusemodel=[ddecal1.center]",
            "ddecal2.h5parm=cal_reuse.h5",
            "ddecal2.subtract=true",
            "msout=.",
            "msout.datacolumn=DATA_REUSE",
        ]
    )

    assert_taql(
        "select from (select abs(gsumsqr(DATA_REUSE - DATA_REF)) as diff "
        f"from {MSIN}) where diff>1.0e-6"
    )


def test_all_model_sources(idgpredict_env, copy_data_to_model_data):
    """Test DDECal with all model data source types enabled"""

    # Multiply MODEL_DATA by 42
    taqlcommand_run = f"update {MSIN} set MODEL_DATA=DATA*42"
    check_output([tcf.TAQLEXE, "-noph", taqlcommand_run])

    # ddecal2 uses the following model data sources:
    # - From the skymodel {SKYMODEL}: Directions 'ra_off' and 'dec_off'.
    # - From foursources-model.fits, using IDGPredict: The directions as
    #   specified in foursources.reg.
    # - From the input MS: The visibilities in the 'MODEL_DATA' column.
    # - From DPBuffer: The 'center' visibilities which ddecal1 added.
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal1,ddecal2]",
            f"ddecal1.sourcedb={SKYMODEL}",
            "ddecal1.directions=center",
            "ddecal1.onlypredict=true",
            "ddecal1.keepmodel=true",
            f"ddecal2.sourcedb={SKYMODEL}",
            "ddecal2.directions=[[ra_off],[dec_off]]",
            f"ddecal2.idg.regions={tcf.DDECAL_RESOURCEDIR}/foursources.reg",
            "ddecal2.idg.images=[foursources-model.fits]",
            "ddecal2.modeldatacolumns=[MODEL_DATA]",
            "ddecal2.reusemodel=[ddecal1.center]",
            "ddecal2.solutions_per_direction=[1,2,4,4,4,1,4,2]",
            "ddecal2.solint=4",
            "ddecal2.h5parm=ddecal2.h5parm",
            "ddecal2.solveralgorithm=directioniterative",
        ]
    )


@pytest.mark.parametrize(
    "mode,soltab,soltab2",
    [
        ("fulljones", "amplitude000", "phase000"),
        ("diagonal", "phase000", ""),
        ("scalarphase", "phase000", ""),
        ("tec", "tec000", ""),
        ("tecandphase", "tec000", "phase000"),
    ],
)
def test_h5parm_initial_solutions(mode, soltab, soltab2):
    """Test using initial solutions from H5Parm."""

    use_approximate_tec = "tec" in mode

    soltabs = f"[{soltab}]" if not soltab2 else f"[{soltab},{soltab2}]"

    # Generate H5 with initial solutions
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.h5parm=initial-solutions.h5",
            f"ddecal.mode={mode}",
            "ddecal.approximatetec=True" if use_approximate_tec else "",
        ]
    )

    # Reuse those solutions as initial conditions
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            "ddecal.h5parm=solutions.h5",
            f"ddecal.mode={mode}",
            "ddecal.initialsolutions.h5parm=initial-solutions.h5",
            f"ddecal.initialsolutions.soltab={soltabs}",
            "ddecal.approximatetec=True" if use_approximate_tec else "",
        ]
    )


def test_h5parm_initial_solutions_with_dd_intervals():
    """
    Test using initial solutions from H5Parm, with direction-dependent solution intervals.
    """

    # Generate H5 with initial solutions
    run_dp3(
        [
            f"msin={MSIN}",
            f"msout=",
            f"steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.h5parm=initial_solutions.h5",
            "ddecal.mode=diagonal",
        ]
    )

    # Reuse those solutions as initial conditions
    run_dp3(
        [
            f"msin={MSIN}",
            f"msout=",
            f"steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.h5parm=solutions.h5",
            "ddecal.mode=scalarphase",
            f"ddecal.initialsolutions.h5parm=initial_solutions.h5",
            f"ddecal.initialsolutions.soltab=[phase000]",
            "ddecal.solutions_per_direction=[1,4,1,1]",
            "ddecal.solint=4",
            "ddecal.solveralgorithm=directioniterative",
        ]
    )
