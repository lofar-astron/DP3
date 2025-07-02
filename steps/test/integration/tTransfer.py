# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from subprocess import check_call

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, get_taql_result, run_in_tmp_path, untar

"""
Integration tests for the Transfer step.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tFlagTransfer.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"  # This MS has 8 freq channels and 6 time steps
MSAVG = "avg_out.ms"
MSOUT = "out.ms"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


@pytest.fixture
def ms_with_fewer_baselines():
    check_call(  # Filter, average & flag 1st and 6th channel of tNDPPP-generic.MS
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msin.useflag=false",
            f"msout={MSAVG}",
            "steps=[filter,averager,preflagger]",
            "filter.baseline=CS*&",
            "filter.remove=true",
            "averager.timestep=3",
            "preflagger.chan=[0, 5]",
        ]
    )


def test_flag_transfer():
    """
    Check whether flags of a single time slot in a low-resolution dataset are correctly
    transfered to three time slots in a higher-resolution dataset.
    """
    check_call(  # Average & flag 2nd time slot of tNDPPP-generic.MS
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msin.useflag=false",
            f"msout={MSAVG}",
            "steps=[averager,preflagger]",
            "averager.timestep=3",
            "averager.freqstep=4",
            "preflagger.timeslot=[1]",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[transfer]",
            f"transfer.source_ms={MSAVG}",
            "transfer.flags=True",
        ]
    )

    # Assert that flags have been transfered to 3 out of 6 time slots
    taql_command = f"select distinct TIME from {MSOUT} where all(FLAG=True)"
    assert_taql(taql_command, 3)
    taql_command = f"select distinct TIME from {MSOUT}"
    assert_taql(taql_command, 6)


def test_data_transfer():
    """
    Check whether data of a single time slot in a low-resolution dataset are correctly
    transfered to three time slots in a higher-resolution dataset.
    """
    check_call(  # Average & flag 2nd time slot of tNDPPP-generic.MS
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msin.useflag=false",
            f"msout={MSAVG}",
            "steps=[averager]",
            "averager.timestep=2",
            "averager.freqstep=4",
        ]
    )

    # Ensure that 1 out of 3 averaged time slots has non-zero data
    check_call([tcf.TAQLEXE, f"update {MSAVG} set DATA=0"])
    check_call(
        [
            tcf.TAQLEXE,
            f"update {MSAVG} set DATA=1 where TIME > [select distinct TIME from {MSAVG}][1]",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[transfer]",
            f"transfer.source_ms={MSAVG}",
            "transfer.data=True",
        ]
    )

    # Assert that data have been transfered to 3 out of 6 time slots
    taql_command = f"select distinct TIME from {MSOUT} where all(DATA=1)"
    assert_taql(taql_command, 2)
    taql_command = f"select distinct TIME from {MSOUT}"
    assert_taql(taql_command, 6)


def test_simultaneous_data_and_flag_transfer():
    """
    Check whether data of a single time slot in a low-resolution dataset are correctly
    transfered to three time slots in a higher-resolution dataset.
    """
    check_call(  # Average & flag 2nd time slot of tNDPPP-generic.MS
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msin.useflag=false",
            f"msout={MSAVG}",
            "steps=[averager,preflagger]",
            "averager.timestep=2",
            "averager.freqstep=4",
            "preflagger.timeslot=[0,2]",
        ]
    )

    # Ensure that 1 out of 3 averaged time slots has non-zero data
    check_call([tcf.TAQLEXE, f"update {MSAVG} set DATA=0"])
    check_call(
        [
            tcf.TAQLEXE,
            f"update {MSAVG} set DATA=1 where TIME > [select distinct TIME from {MSAVG}][1]",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[transfer]",
            f"transfer.source_ms={MSAVG}",
            "transfer.data=True",
            "transfer.flags=True",
        ]
    )

    # Assert that data have been transfered to 2 out of 6 time slots and flags to 4 out of 6 time slots
    taql_command = f"select distinct TIME from {MSOUT} where all(DATA=1)"
    assert_taql(taql_command, 2)
    taql_command = f"select distinct TIME from {MSOUT} where all(FLAG=True)"
    assert_taql(taql_command, 4)
    taql_command = f"select distinct TIME from {MSOUT}"
    assert_taql(taql_command, 6)


def test_chan_flag_transfer():
    """
    Check whether flags of a single channel are correctly transfered from a
    low-resolution to higher-resolution dataset.
    """
    check_call(  # Average & flag 1st channel of tNDPPP-generic.MS
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msin.useflag=false",
            f"msout={MSAVG}",
            "steps=[averager,preflagger]",
            "averager.timestep=3",
            "averager.freqstep=4",
            "preflagger.chan=[0]",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[transfer]",
            f"transfer.source_ms={MSAVG}",
            "transfer.flags=True",
        ]
    )

    # Assert that flags have been transfered to the first 4 out of 8 channels
    taql_command = f"select gntrues(FLAG[,0]) from out.ms"
    flag_counts_str = get_taql_result(taql_command).split("\n")[1]

    expected_flag_counts_str = "[168, 168, 168, 168, 0, 0, 0, 0]"
    assert flag_counts_str == expected_flag_counts_str


def test_flag_transfer_without_freq_avg():
    """
    Check whether flags are properly transfered when frequency averaging is not used.
    """
    check_call(  # Average & flag 1st and 6th channel of tNDPPP-generic.MS
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msin.useflag=false",
            f"msout={MSAVG}",
            "steps=[averager,preflagger]",
            "averager.timestep=3",
            "preflagger.chan=[0, 5]",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[transfer]",
            f"transfer.source_ms={MSAVG}",
            "transfer.flags=True",
        ]
    )

    # Assert that flags have been transfered to the 1st and 6th channels
    taql_command = f"select gntrues(FLAG[,0]) from out.ms"
    flag_counts_str = get_taql_result(taql_command).split("\n")[1]

    expected_flag_counts_str = "[168, 0, 0, 0, 0, 168, 0, 0]"
    assert flag_counts_str == expected_flag_counts_str


def test_flag_transfer_with_fewer_baselines(ms_with_fewer_baselines):
    """
    Check whether flags are properly transfered when frequency averaging is not used.
    """
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[transfer]",
            f"transfer.source_ms={MSAVG}",
            "transfer.flags=True",
            "transfer.baseline=CS*&",
        ]
    )


def test_exception_with_fewer_baselines(ms_with_fewer_baselines):
    """
    Check whether an exception is thrown when there's fewer baselines in source_ms.
    """
    with pytest.raises(Exception):
        check_call(
            [
                tcf.DP3EXE,
                "checkparset=1",
                f"msin={MSIN}",
                f"msout={MSOUT}",
                "steps=[transfer]",
                f"transfer.source_ms={MSAVG}",
                "transfer.flags=True",
            ]
        )
