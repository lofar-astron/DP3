# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import re
import shutil
import uuid
from subprocess import check_call

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import untar_ms, get_taql_result, check_output

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


def test_chunking():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=chunktest.ms",
            "msout.chunkduration=18",
            "steps=[]",
        ]
    )
    # The integration time of this set is 10 seconds. With chunking after 18 sec, this
    # means 2 timeslots per ms. There are 6 timeslots in the ms, hence there should be
    # 3 output measurement sets:
    assert os.path.exists("chunktest-000.ms")
    assert os.path.exists("chunktest-001.ms")
    assert os.path.exists("chunktest-002.ms")
    
    # Each should have two timesteps:
    taql_command = f"select unique TIME from chunktest-000.ms"
    result = get_taql_result(taql_command)
    assert result == 'Unit: s\n29-Mar-2013/13:59:53.007\n29-Mar-2013/14:00:03.021'
    
    taql_command = f"select unique TIME from chunktest-001.ms"
    result = get_taql_result(taql_command)
    assert result == 'Unit: s\n29-Mar-2013/14:00:13.035\n29-Mar-2013/14:00:23.049'

    taql_command = f"select unique TIME from chunktest-002.ms"
    result = get_taql_result(taql_command)
    assert result == 'Unit: s\n29-Mar-2013/14:00:33.063\n29-Mar-2013/14:00:43.076'
