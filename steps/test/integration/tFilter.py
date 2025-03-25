# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call

import pytest

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Integration tests for the Filter step.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/Filter.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


def test_filter_added_station():
    """
    Adds two stations and filters out the first added station.
    """
    MSOUT = "out.ms"
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[stationadder,filter]",
            "stationadder.stations={CSNEW:[CS001HBA0,CS002HBA0],RSNEW:[RS106HBA,RS208HBA]}",
            "stationadder.average=False",  # Necessary to keep uvw consistent for uvw compression
            "filter.baseline=!CSNEW",
            "filter.remove=true",
        ]
    )

    taql_command = f"select from {MSOUT}/ANTENNA where NAME='CSNEW'"
    assert_taql(taql_command, 0)

    taql_command = f"select from {MSOUT}/ANTENNA where NAME='RSNEW'"
    assert_taql(taql_command, 1)
