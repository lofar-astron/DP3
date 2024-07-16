# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
from subprocess import check_call, check_output

import pytest
from testconfig import TAQLEXE


def assert_taql(command, expected_rows=0):
    result = (
        check_output([TAQLEXE, "-noph", "-nopa", command]).decode().strip()
    )
    assert result == f"select result of {expected_rows} rows"


def untar(source):
    if not os.path.isfile(source):
        raise IOError(
            f"Not able to find {source} containing test input files."
        )
    check_call(["tar", "xf", source])


def get_taql_result(command):
    """Get the output of a taql command"""
    result = (
        check_output([TAQLEXE, "-noph", "-nopr", command]).decode().strip()
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
