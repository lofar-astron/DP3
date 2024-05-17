# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import sys
from subprocess import check_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Tests for applying the beam model.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tDemix.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN_DEMIX = "tDemix.in_MS"
MSIN_GENERIC = "tNDPPP-generic.MS"

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
    "demix.demixfreqstep=64",
    "demix.demixtimestep=10",
    "demix.instrumentmodel='tDemix_tmp/instrument'",
    "demix.subtractsources=[CasA]",
]

skymodel_arg = "demix.skymodel='tDemix_tmp/{}'"


@pytest.fixture()
def demix_ms(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN_DEMIX}.tgz")
    check_call(
        [
            tcf.MAKESOURCEDBEXE,
            "in=tDemix_tmp/sky.txt",
            "out=tDemix_tmp/sourcedb",
        ]
    )


@pytest.fixture()
def generic_ms(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN_GENERIC}.tgz")


@pytest.mark.parametrize("skymodel", ["sky.txt", "sourcedb"])
def test_without_target(demix_ms, skymodel):
    check_call(
        [
            tcf.DP3EXE,
            "demix.ignoretarget=true",
            "demix.freqstep=64",
            "demix.timestep=10",
            skymodel_arg.format(skymodel),
        ]
        + common_args
    )

    # Compare some columns of the output MS with the reference output.
    taql_command = f"select from tDemix_out.MS t1, tDemix_tmp/tDemix_ref1.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME"
    assert_taql(taql_command)


@pytest.mark.parametrize("skymodel", ["sky.txt", "sourcedb"])
def test_with_target_projected_away(demix_ms, skymodel):
    check_call(
        [
            tcf.DP3EXE,
            "demix.ignoretarget=false",
            "demix.freqstep=64",
            "demix.timestep=10",
            skymodel_arg.format(skymodel),
        ]
        + common_args
    )

    # Compare some columns of the output MS with the reference output.
    taql_command = f"select from tDemix_out.MS t1, tDemix_tmp/tDemix_ref2.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME"
    assert_taql(taql_command)


@pytest.mark.parametrize("skymodel", ["sky.txt", "sourcedb"])
def test_with_target(demix_ms, skymodel):
    check_call(
        [
            tcf.DP3EXE,
            "demix.targetsource=CIZA.SP1A.FITS.pbcor_patch_s537",
            "demix.freqstep=32",
            "demix.timestep=5",
            "demix.maxiter=100",
            skymodel_arg.format(skymodel),
        ]
        + common_args
    )

    # Compare some columns of the output MS with the reference output.
    taql_command = f"select from tDemix_out.MS t1, tDemix_tmp/tDemix_ref3.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME"
    assert_taql(taql_command)


def test_time_freq_resolution(demix_ms):
    check_call(
        [
            tcf.DP3EXE,
            "demix.ignoretarget=true",
            "demix.freqstep=64",
            "demix.timestep=10",
            "demix.demixfreqresolution=200kHz",
            "demix.demixtimeresolution=10.0139",
            "demix.skymodel=tDemix_tmp/sky.txt",
        ]
        + common_args
    )

    # Compare some columns of the output MS with the reference output.
    taql_command = f"select from tDemix_out.MS t1, tDemix_tmp/tDemix_ref1.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME"
    assert_taql(taql_command)


@pytest.mark.skip(reason="Issue #268 is not fixed yet")
def test_multiple_baseline_selection(generic_ms):
    # Test for https://git.astron.nl/RD/DP3/-/issues/268, where the baseline
    # selection for demix was already applied using msin and filter.
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN_GENERIC}",
            "msin.baseline=CS*&",
            "steps=[filter,demix,null]",
            "filter.baseline=CS*&",
            "filter.remove=true",
            "demix.baseline=CS*&",
            f"demix.skymodel={MSIN_GENERIC}/sky",
        ]
    )
