# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call

import pytest

sys.path.append(".")

import testconfig as tcf
from utils import run_in_tmp_path, untar

# FIXME: consider merging into tGainCal?

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tReadOnly.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSGAINCAL = "tGainCal.tab"


@pytest.fixture(autouse=True)
def source_env():
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar(f"{tcf.SRCDIR}/{MSGAINCAL}.tgz")


def test_read_only():
    check_call(["chmod", "a-w", "-R", MSIN])
    # Sky table needs to be writeable (to acquire a lock?)
    check_call(["chmod", "u+w", "-R", f"{MSIN}/sky"])
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[gaincal]",
            f"gaincal.parmdb={MSIN}/inst-diagonal",
            f"gaincal.sourcedb={MSIN}/sky",
            "gaincal.caltype=diagonal",
            "gaincal.parmdb=gaincal.h5",
            "msout=",
        ]
    )
    # Reset permissions to clean up tmp directory
    check_call(["chmod", "a+w", "-R", MSIN])
