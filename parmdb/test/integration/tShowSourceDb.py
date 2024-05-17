import pytest
import sys
from subprocess import check_call, check_output

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

MSIN = "tDemix.in_MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


def test_skymodel_sourcedb_roundtrip():
    """Check that skymodel in default format is reproduced after makesourcedb and showsourcedb"""

    # sky.txt is not in the default format, create a skymodel in default format by going through sourcedb
    check_call([tcf.MAKESOURCEDBEXE, "in=tDemix_tmp/sky.txt", "out=sourcedb"])
    # The first line of showsourcedb is the application's announcement
    skymodel_defaultformat_input = (
        check_output([tcf.SHOWSOURCEDBEXE, "in=sourcedb", "mode=skymodel"])
        .decode("utf-8")
        .split("\n", 1)[-1]
    )
    with open("tDemix_tmp/sky_defaultformat.txt", "w") as f:
        f.write(skymodel_defaultformat_input)

    # Now do the roundtrip test: make sourcedb and print the result in the default format
    check_call(
        [
            tcf.MAKESOURCEDBEXE,
            "in=tDemix_tmp/sky_defaultformat.txt",
            "out=sourcedb_defaultformat",
        ]
    )
    skymodel_defaultformat_output = (
        check_output(
            [tcf.SHOWSOURCEDBEXE, "in=sourcedb_defaultformat", "mode=skymodel"]
        )
        .decode("utf-8")
        .split("\n", 1)[-1]
    )

    assert skymodel_defaultformat_input == skymodel_defaultformat_output
