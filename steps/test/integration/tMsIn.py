# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import shutil

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call

import pytest

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tMsIn.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


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


def test_autoweight():
    """Check that DP3 runs successfully when autoweight is set and there is a step which requires weights"""

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=out.ms",
            "steps=[average]",
            "msin.autoweight=true",
            "msin.forceautoweight=true",
            "average.timestep=5",
        ]
    )


@pytest.fixture
def make_multiple_bands():
    # Split input MS into four frequency chunks (to test the MultiMsReader)
    for startchan in range(0, 8, 2):
        check_call(
            [
                tcf.DP3EXE,
                f"msin={MSIN}",
                f"msout=chunk{startchan//2}.MS",
                f"msin.startchan={startchan}",
                f"msin.nchan=2",
                "steps=[]",
            ]
        )


def test_multimsin_missingdata(make_multiple_bands):
    shutil.rmtree("chunk2.MS")

    check_call(
        [
            tcf.DP3EXE,
            "msin=[chunk0.MS,chunk1.MS,DOESNOTEXIST,chunk3.MS]",
            "steps=[]",
            "msout=joined.MS",
            "msin.missingdata=true",
            "msin.orderms=false",
        ]
    )

    taql_command = f"select from joined.MS::SPECTRAL_WINDOW where shape(CHAN_FREQ)[0] == 8"
    assert_taql(taql_command, 1)


@pytest.mark.parametrize("order_ms", [True, False])
def test_multimsin_orderms(make_multiple_bands, order_ms):
    check_call(
        [
            tcf.DP3EXE,
            "msin=[chunk0.MS,chunk2.MS,chunk1.MS,chunk3.MS]",
            "steps=[]",
            "msout=joined.MS",
            f"msin.orderms={order_ms}",
        ]
    )

    taql_command = f"select from joined.MS::SPECTRAL_WINDOW WHERE CHAN_FREQ[4] > CHAN_FREQ[3]"
    if order_ms:
        assert_taql(taql_command, 1)
    else:
        assert_taql(taql_command, 0)
