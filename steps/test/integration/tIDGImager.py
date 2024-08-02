import os
import sys
from subprocess import check_call

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import get_taql_result, run_in_tmp_path, untar

"""
Tests for testing IDGImager

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tIDGImager.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    get_taql_result(f"UPDATE {MSIN} SET DATA=1")


def test_idg_image_generation_cpu():
    """
    The test checks if images are generated and
    some standard input parameters are parsed correctly
    such as the file name.
    The measurement set contains 5 time samples.
    """
    check_call(
        [tcf.DP3EXE, f"msin={MSIN}", "msout=out.MS", "steps=[idgimager]"]
    )
    for k in range(5):
        assert os.path.exists(f"image_t{k}.fits")

    test_name = "testname_f"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=out.MS",
            "msout.overwrite=True",
            f"idgimager.image_name={test_name}%t.fits",
            "steps=[idgimager]",
        ]
    )
    for k in range(5):
        assert os.path.exists(f"{test_name}{k}.fits")
