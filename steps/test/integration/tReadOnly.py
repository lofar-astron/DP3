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
from utils import untar_ms

# FIXME: consider merging into tGainCal?

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tReadOnly.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSGAINCAL = "tGainCal.tab"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar_ms(f"{tcf.SRCDIR}/{MSGAINCAL}.tgz")

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


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
