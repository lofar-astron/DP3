import pytest
import shutil
import os
import sys
import uuid
from subprocess import check_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import get_taql_result, untar_ms

"""
Tests for testing IDGImager

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tIDGImager.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(tmpdir_factory):
    tmpdir = str(tmpdir_factory.mktemp("data"))
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    get_taql_result(f"UPDATE {MSIN} SET DATA=1")

    # Tests are executed here
    yield

    # Post-test: clean up
    shutil.rmtree(tmpdir)


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
