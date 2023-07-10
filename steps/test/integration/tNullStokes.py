import pytest
import shutil
import os
import sys
import uuid
from subprocess import check_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, untar_ms

"""
Tests for nullstokes (zeroing out the Stokes parameters Q and/or U).

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tNullStokes.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSOUT = "test_data.MS"
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


@pytest.fixture()
def create_model_data():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[]",
        ]
    )


def test_stokes_Qzero(create_model_data):
    check_call(
        [
            tcf.DP3EXE,
            "msin=test_data.MS",
            "msout=.",
            "msout.datacolumn=STOKES_QDATA",
            "steps=[nullstokes]",
            "nullstokes.modify_q=true",
        ]
    )

    # Checking that all Stokes Q values are equal to zero.
    taql_command = f"select from (select sum(abs(mscal.stokes(STOKES_QDATA, 'Q'))) as Qtot from {MSOUT}) where Qtot > 1e-5"
    assert_taql(taql_command)

    # Checking that Stokes I, U, V have not been changed
    taql_command = f"select from (select mscal.stokes(STOKES_QDATA, 'I') as a, mscal.stokes(DATA, 'I') as b from {MSOUT}) where not all(near(a, b))"
    assert_taql(taql_command)

    taql_command = f"select from (select mscal.stokes(STOKES_QDATA, 'U') as a, mscal.stokes(DATA, 'U') as b from {MSOUT}) where not all(near(a, b))"
    assert_taql(taql_command)

    taql_command = f"select from (select mscal.stokes(STOKES_QDATA, 'V') as a, mscal.stokes(DATA, 'V') as b from {MSOUT}) where not all(near(a, b))"
    assert_taql(taql_command)


def test_stokes_Uzero(create_model_data):
    check_call(
        [
            tcf.DP3EXE,
            "msin=test_data.MS",
            "msout=.",
            "msout.datacolumn=STOKES_UDATA",
            "steps=[nullstokes]",
            "nullstokes.modify_u=true",
        ]
    )

    # Checking that all Stokes U values are equal to zero.
    taql_command = f"select from (select sum(abs(mscal.stokes(STOKES_UDATA, 'U'))) as Utot from {MSOUT}) where Utot > 1e-5"
    assert_taql(taql_command)

    # Checking that Stokes I, Q, V have not been changed
    taql_command = f"select from (select mscal.stokes(STOKES_UDATA, 'I') as a, mscal.stokes(DATA, 'I') as b from {MSOUT}) where not all(near(a, b))"
    assert_taql(taql_command)

    taql_command = f"select from (select mscal.stokes(STOKES_UDATA, 'Q') as a, mscal.stokes(DATA, 'Q') as b from {MSOUT}) where not all(near(a, b))"
    assert_taql(taql_command)

    taql_command = f"select from (select mscal.stokes(STOKES_UDATA, 'V') as a, mscal.stokes(DATA, 'V') as b from {MSOUT}) where not all(near(a, b))"
    assert_taql(taql_command)
