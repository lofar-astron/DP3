# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import re

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call

import pytest

sys.path.append(".")

import testconfig as tcf
from utils import (
    assert_taql,
    check_output,
    get_taql_result,
    run_in_tmp_path,
    untar,
)

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


def test_chunking():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=chunktest.ms",
            "msout.chunkduration=18",
            "steps=[filter]",
            "filter.baseline=[CR]S*&",
            "filter.remove=True",
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
    assert (
        result == "Unit: s\n29-Mar-2013/13:59:53.007\n29-Mar-2013/14:00:03.021"
    )

    taql_command = f"select unique TIME from chunktest-001.ms"
    result = get_taql_result(taql_command)
    assert (
        result == "Unit: s\n29-Mar-2013/14:00:13.035\n29-Mar-2013/14:00:23.049"
    )

    taql_command = f"select unique TIME from chunktest-002.ms"
    result = get_taql_result(taql_command)
    assert (
        result == "Unit: s\n29-Mar-2013/14:00:33.063\n29-Mar-2013/14:00:43.076"
    )


def test_write_thread_enabled():
    """Assert that the threaded writing is enabled.

    This requires the output step to be the last step."""

    result = check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=out.MS",
            "steps=[]",
        ]
    )

    assert re.search(b"use thread: *true", result)
    assert re.search(
        b"(1[0-9]| [ 0-9])[0-9]\\.[0-9]% \\([ 0-9]{5} [m ]s\\) Creating task\n",
        result,
    )
    assert re.search(
        b"(1[0-9]| [ 0-9])[0-9]\\.[0-9]% \\([ 0-9]{5} [m ]s\\) Writing \\(threaded\\)\n",
        result,
    )


def test_write_thread_disabled():
    """Assert that the threaded writing is disabled.

    This requires the output step not to be the last step."""

    # Use a split so the default msout step won't be created.
    result = check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[split]",
            "split.steps=[out,average]",
            "split.replaceparms=[out.name]",
            "out.name=[out.MS]",
        ]
    )

    assert re.search(b"use thread: *false", result)
    assert (
        re.search(
            b"(1[0-9]| [ 0-9])[0-9]\\.[0-9]% \\([ 0-9]{5} [m ]s\\) Creating task\n",
            result,
        )
        == None
    )
    assert re.search(
        b"(1[0-9]| [ 0-9])[0-9]\\.[0-9]% \\([ 0-9]{5} [m ]s\\) Writing\n",
        result,
    )


def get_directory_size(p):
    """
    Computes the sum of the file sizes in directory 'p'.
    Does not traverse 'p' recursively.
    """
    return sum(
        os.path.getsize(p + "/" + f)
        for f in os.listdir(p)
        if os.path.isfile(p + "/" + f)
    )


def test_dysco():
    # This just checks if it runs without errors
    check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[]",
            "msout=dysco.MS",
            "msout.storagemanager=Dysco",
        ]
    )


def test_uvw_compression():
    check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[]",
            "msout=out1.MS",
            "msout.uvwcompression=True",
        ]
    )

    check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[]",
            "msout=out2.MS",
            "msout.uvwcompression=False",
        ]
    )

    taql_command = "select from out1.MS AS out1, out2.MS AS out2 where not all(out1.UVW~=out2.UVW)"
    assert_taql(taql_command)
    # I measured these sizes:
    # - Without uvw compression: 1911382 bytes
    # - With uvw compression: 1813862 bytes
    # Hence comression saves 97520 bytes. The test is a bit more flexible to allow some changes to occur.
    assert get_directory_size("out1.MS") + 87500 <= get_directory_size(
        "out2.MS"
    )


def test_scalar_flags():
    check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[]",
            "msout=out1.MS",
            "msout.scalarflags=True",
        ]
    )

    check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[]",
            "msout=out2.MS",
            "msout.scalarflags=False",
        ]
    )

    taql_command = "select from out1.MS AS out1, out2.MS AS out2 where not all(out1.FLAG==out2.FLAG)"
    assert_taql(taql_command)
    # Currently it saves 131169 bytes. The test is a bit more flexible to allow some changes to occur.
    assert get_directory_size("out1.MS") + 125000 < get_directory_size(
        "out2.MS"
    )


def test_antenna_compression():
    check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[]",
            "msout=out1.MS",
            "msout.antennacompression=True",
        ]
    )

    check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[]",
            "msout=out2.MS",
            "msout.antennacompression=False",
        ]
    )

    taql_command = "select from out1.MS AS out1, out2.MS AS out2 where not all(out1.ANTENNA1==out2.ANTENNA1) or not all(out1.ANTENNA2==out2.ANTENNA2)"
    assert_taql(taql_command)
    # I measured these sizes:
    # - Without antennacompression compression: 1911382 bytes
    # - With antennacompression compression: 1845518 bytes
    # Hence compression saves 65864 bytes. The test is a bit more flexible to allow some changes to occur.
    assert get_directory_size("out1.MS") + 64500 <= get_directory_size(
        "out2.MS"
    )
