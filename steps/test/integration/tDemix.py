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
Tests for applying the beam model.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tDemix.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tDemix.in_MS"
CWD = os.getcwd()

common_args = [
    "msin=tDemix_tmp/tDemix.MS",
    "msout=tDemix_out.MS",
    "msout.overwrite=True",
    "msout.tilesize=1",
    "msin.datacolumn=DATA",
    "msout.datacolumn=DATA",
    "steps=[demix]",
    "demix.type=demixer",
    "demix.corrtype=cross",
    "demix.baseline='CS00[0-9]HBA0&'",
    "demix.freqstep=64",
    "demix.timestep=10",
    "demix.demixfreqstep=64",
    "demix.demixtimestep=10",
    "demix.instrumentmodel='tDemix_tmp/instrument'",
    "demix.subtractsources=[CasA]",
]

skymodel_arg="demix.skymodel='tDemix_tmp/{}'"

@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    check_call([tcf.MAKESOURCEDBEXE, "in=tDemix_tmp/sky.txt", "out=tDemix_tmp/sourcedb"])

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


@pytest.mark.parametrize("skymodel", ['sky.txt', 'sourcedb'])
def test_without_target(skymodel):

    check_call([tcf.DP3EXE, "demix.ignoretarget=true", skymodel_arg.format(skymodel)] + common_args)

    # Compare some columns of the output MS with the reference output.
    taql_command = f"select from tDemix_out.MS t1, tDemix_tmp/tDemix_ref1.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  not all(t1.LOFAR_FULL_RES_FLAG = t2.LOFAR_FULL_RES_FLAG)  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME"
    assert_taql(taql_command)

@pytest.mark.parametrize("skymodel", ['sky.txt', 'sourcedb'])
def test_with_target_projected_away(skymodel):
    check_call([tcf.DP3EXE, "demix.ignoretarget=false", skymodel_arg.format(skymodel)] + common_args)

    # Compare some columns of the output MS with the reference output.
    taql_command = f"select from tDemix_out.MS t1, tDemix_tmp/tDemix_ref2.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  not all(t1.LOFAR_FULL_RES_FLAG = t2.LOFAR_FULL_RES_FLAG)  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME"
    assert_taql(taql_command)

@pytest.mark.parametrize("skymodel", ['sky.txt', 'sourcedb'])
def test_with_target(skymodel):
    check_call(
        [
            tcf.DP3EXE,
            "demix.target=CIZA.SP1A.FITS.pbcor_patch_s537",
            "demix.freqstep=32",
            "demix.timestep=5",
            skymodel_arg.format(skymodel),
        ]
        + common_args
    )

    # Compare some columns of the output MS with the reference output.
    taql_command = f"select from tDemix_out.MS t1, tDemix_tmp/tDemix_ref2.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  not all(t1.LOFAR_FULL_RES_FLAG = t2.LOFAR_FULL_RES_FLAG)  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME"
    assert_taql(taql_command)
