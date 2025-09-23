# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import subprocess
from subprocess import Popen, check_call, check_output, run

import pytest
from testconfig import DP3EXE, TAQLEXE

COMMON_DP3_ARGUMENTS = ["checkparset=1", "numthreads=1"]


def assert_taql(command, expected_rows=0):
    result = (
        subprocess.check_output([TAQLEXE, "-noph", "-nopa", command])
        .decode()
        .strip()
    )
    assert result == f"select result of {expected_rows} rows"


def untar(source):
    if not os.path.isfile(source):
        raise IOError(
            f"Not able to find {source} containing test input files."
        )
    subprocess.check_call(["tar", "xf", source])


def get_taql_result(command):
    """Get the output of a taql command"""
    result = (
        subprocess.check_output([TAQLEXE, "-noph", "-nopr", command])
        .decode()
        .strip()
    )
    return result


@pytest.fixture()
def run_in_tmp_path(tmp_path):
    """
    Creates a temporary directory, runs the test in it, and removes the
    directory.
    """
    # 'tmp_path' is a base fixture from Pytest that already
    # does everything else, including cleaning up.
    os.chdir(tmp_path)


def run_dp3(arguments):
    """
    Run DP3 with the given arguments.
    Prepend common DP3 test arguments (checkparset=1, numthreads=1).
    For testing DP3 with multiple threads, add numthreads=<N> to the arguments.
    """
    all_arguments = COMMON_DP3_ARGUMENTS + arguments
    print("DP3 " + " ".join(all_arguments))
    result = run([DP3EXE] + all_arguments, stdout=subprocess.PIPE)
    if result.returncode != 0:
        # The output is shortened and not well displayed in an exception, so
        # the output is explicitly printed:
        stdout_string = result.stdout.decode("utf-8")
        print(f"DP3 failed. Output was:\n{stdout_string}")
        raise subprocess.CalledProcessError(
            f"DP3 failed. Output was:\n{result.stdout}",
            cmd="DP3 " + " ".join(all_arguments),
        )
    return result.stdout


def spawn_dp3(arguments):
    """
    Run DP3 with the given arguments in the background.
    Prepend common DP3 test arguments (checkparset=1, numthreads=1).
    For testing DP3 with multiple threads, add numthreads=<N> to the arguments.
    """
    all_arguments = COMMON_DP3_ARGUMENTS + arguments
    print("DP3 " + " ".join(all_arguments))
    child = Popen(
        [DP3EXE] + all_arguments,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return child
