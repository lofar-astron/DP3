import pytest
import shutil
import os
import sys
import uuid
from subprocess import check_call
from subprocess import check_output

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import untar_ms

"""
Tests for the showsourcedb tool

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tShowSourceDb.py` (extended with pytest options of your choice)
- using ctest, see DP3/parmdb/test/integration/CMakeLists.txt
"""

MSIN = "tDemix.in_MS"
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

def test_skymodel_sourcedb_roundtrip():
    check_call([tcf.MAKESOURCEDBEXE, "in=tDemix_tmp/sky.txt", "out=sourcedb"])
    # The first line of showsourcedb is the application's announcement
    out = check_output([tcf.SHOWSOURCEDBEXE, "in=sourcedb", "mode=skymodel"]).decode('utf-8').split('\n',1)[-1]
    data = open("tDemix_tmp/sky.txt", "r").read()
    assert data == out
