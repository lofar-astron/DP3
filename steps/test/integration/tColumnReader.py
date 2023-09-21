# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import shutil
from subprocess import check_call

import sys

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import untar_ms

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tColumnReader.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(tmpdir_factory):
    tmpdir = str(tmpdir_factory.mktemp("data"))
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    shutil.rmtree(tmpdir)


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
