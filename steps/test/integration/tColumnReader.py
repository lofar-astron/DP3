# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import shutil
import sys
from subprocess import check_call

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tColumnReader.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


def test_filter_and_columnreader():
    """Read MS, keep only one baseline, then use columnreader"""

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=.",
            "steps=[filter,columnreader]",
            "columnreader.column=DATA",
            "filter.baseline='0&1'",
            "msout=",
        ]
    )
