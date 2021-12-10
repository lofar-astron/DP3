# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import shutil
import uuid
from subprocess import check_call

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import untar_ms, assert_taql

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tMsIn.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
CWD = os.getcwd()


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


def test_update_flags():
    """Assert that the column FLAG is written in / read from the right flagcolumn, if specified"""

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[uvwflag]",
            "uvwflag.ulambdamax=0",
            "uvwflag.ulambdamin=100000000000000",
            "msout.flagcolumn=MODIFIED_FLAGS",
        ]
    )

    # Assert that FLAG column is not changed
    taql_command = f"select from {MSIN} where all(FLAG=True)"
    assert_taql(taql_command, 52)

    # Assert that the output flags are directed to the MODIFIED_FLAGS column
    taql_command = f"select from {MSIN} where all(MODIFIED_FLAGS=True)"
    assert_taql(taql_command, 168)

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=out.MS",
            "steps=[]",
            "msin.flagcolumn=MODIFIED_FLAGS",
        ]
    )

    # Assert that the right flag column is read from the input
    taql_command = f"select from out.MS where all(FLAG=True)"
    assert_taql(taql_command, 168)


def test_create_flag_column():
    """Assert that the output FLAG is written in the right flagcolumn also when no steps are performed"""

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[]",
            "msout.flagcolumn=DUPLICATED_FLAGS",
        ]
    )

    # Assert that the right flag column is created
    taql_command = f"select from {MSIN} where all(DUPLICATED_FLAGS=True)"
    assert_taql(taql_command, 52)


def test_unflag_all():
    """Assert that the output FLAG is written in the FLAG column, if flagcolumn is not specified"""

    check_call([tcf.TAQLEXE, f"update {MSIN} set FLAG=False"])

    # Assert that all FLAG are set to false
    taql_command = f"select from {MSIN} where all(FLAG=False)"
    assert_taql(taql_command, 168)

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[preflag1,preflag2]",
            "preflag1.corrtype=cross",
            "preflag2.corrtype=auto",
        ]
    )

    # Assert that all FLAG are set to True
    taql_command = f"select from {MSIN} where all(FLAG=True)"
    assert_taql(taql_command, 168)
