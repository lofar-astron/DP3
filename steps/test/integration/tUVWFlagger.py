# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import re
from subprocess import check_call, check_output

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as config
from utils import assert_taql, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tUVWFlagger..py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{config.RESOURCEDIR}/{MSIN}.tgz")


def test_update_flags_inplace():
    """Assert that updating the FLAG column works.

    The summary of the set flags is tested (AST-862).
    """

    count_flags_set = f"select from {MSIN} where all(FLAG=True)"

    check_call([config.TAQLEXE, "update", MSIN, "set", "FLAG=False"])
    assert_taql(count_flags_set)

    # The first run all visibilites are flagged.
    result = check_output(
        [
            config.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[uvwflag]",
            "uvwflag.ulambdamax=0",
            "uvwflag.ulambdamin=100000000000000",
        ]
    )
    assert_taql(count_flags_set, 168)
    assert re.search(
        b"\nTotal flagged:   100.000%   \\(1344 out of 1344 visibilities\\)\n\n\n",
        result,
    )

    # The second run nothing is flagged, since everything is already flagged
    result = check_output(
        [
            config.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[uvwflag]",
            "uvwflag.ulambdamax=0",
            "uvwflag.ulambdamin=100000000000000",
        ]
    )
    assert_taql(count_flags_set, 168)
    assert re.search(
        b"\nTotal flagged:     0.000%   \\(0 out of 1344 visibilities\\)\n\n\n",
        result,
    )


def test_update_flags_new_table():
    """Assert that updating the flags in a different output column works."""
    count_flags_set = f"select from {MSIN} where all(FLAG=True)"

    check_call([config.TAQLEXE, "update", MSIN, "set", "FLAG=False"])
    assert_taql(count_flags_set, 0)

    result = check_output(
        [
            config.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[uvwflag]",
            "uvwflag.ulambdamax=0",
            "uvwflag.ulambdamin=100000000000000",
            "msout.flagcolumn=MODIFIED_FLAGS",
        ]
    )

    # Assert that FLAG column is not changed
    assert_taql(count_flags_set, 0)

    assert re.search(
        b"\nTotal flagged:   100.000%   \\(1344 out of 1344 visibilities\\)\n\n\n",
        result,
    )
    count_flags_set = f"select from {MSIN} where all(MODIFIED_FLAGS=True)"
    assert_taql(count_flags_set, 168)

    # The second run compares the flagged changes based on the FLAG column
    # in the input. So it still will again flag 100% of the original data.

    result = check_output(
        [
            config.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[uvwflag]",
            "uvwflag.ulambdamax=0",
            "uvwflag.ulambdamin=100000000000000",
            "msout.flagcolumn=MODIFIED_FLAGS",
        ]
    )
    assert re.search(
        b"\nTotal flagged:   100.000%   \\(1344 out of 1344 visibilities\\)\n\n\n",
        result,
    )
    count_flags_set = f"select from {MSIN} where all(MODIFIED_FLAGS=True)"
    assert_taql(count_flags_set, 168)
