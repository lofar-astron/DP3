import sys
from subprocess import check_call, check_output

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import run_in_tmp_path, untar

"""
Tests for the showsourcedb tool

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tShowSourceDb.py` (extended with pytest options of your choice)
- using ctest, see DP3/parmdb/test/integration/CMakeLists.txt
"""

MSIN_GENERIC = "tNDPPP-generic.MS"
TEXT_MODEL = open(f"{tcf.RESOURCEDIR}/tNDPPP-generic-skymodel.txt").read()


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN_GENERIC}.tgz")


def test_showsourcedb():
    # showsourcedb prints it's version on the first line, so skip this
    skymodel_output = (
        check_output(
            [tcf.SHOWSOURCEDBEXE, "in=tNDPPP-generic.MS/sky", "mode=skymodel"]
        )
        .decode("utf-8")
        .split("\n", 1)[-1]
    )

    assert skymodel_output == TEXT_MODEL
