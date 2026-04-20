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
from utils import run_dp3, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tColumnReader.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSIN_BDA = "tNDPPP-bda.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar(f"{tcf.RESOURCEDIR}/{MSIN_BDA}.tgz")


def test_filter_and_columnreader():
    """Read MS, keep only one baseline, then use columnreader"""

    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[filter,columnreader]",
            "columnreader.column=DATA",
            "filter.baseline='0&1'",
            "msout=",
        ]
    )


def test_bda_ddecal_with_column_readers():
    """
    Perform a basic test of reading in a BDA model data column inside the DDECal step. This test only checks
    if this doesn't crash, the results are not validated.
    """

    # Make a dirty image as template
    try:
        check_call(
            [
                "wsclean",
                "-reorder",
                "-size",
                "1024",
                "1024",
                "-scale",
                "1amin",
                "-no-dirty",
                MSIN_BDA,
            ]
        )
    except:
        pytest.skip(reason="WSClean is not available")

    os.rename("wsclean-image.fits", "wsclean-model.fits")

    check_call(
        [
            "wsclean",
            "-predict",
            "-reorder",
            "-model-column",
            # This doesn't work yet in WSClean, so disable for now:
            # "-model-storage-manager sisco",
            "DIRECTION1",
            MSIN_BDA,
        ]
    )

    check_call(
        [
            "wsclean",
            "-predict",
            "-reorder",
            "-model-column",
            # "-model-storage-manager sisco",
            "DIRECTION2",
            MSIN_BDA,
        ]
    )

    run_dp3(
        [
            f"msin={MSIN_BDA}",
            "msout=",
            "steps=[ddecal]",
            "ddecal.nchan=0",
            "ddecal.modeldatacolumns=[DIRECTION1,DIRECTION2]",
            "ddecal.solint=1000",
            "ddecal.h5parm=solutions.h5",
        ]
    )
