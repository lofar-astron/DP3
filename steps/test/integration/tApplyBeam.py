# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
from subprocess import check_call

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Tests for applying the beam model.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tApplyBeam.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
DISH_MSIN = "tDish.MS"
MSAPPLYBEAM = "tApplyBeam.tab"

"""
The tDish.MS is a reduced version of a MEERKAT dataset, generated with the following command:

DP3 \
msin=MKT_CDFS_2_5_856_963MHz.ms \
msout=tDish.MS \
msout.overwrite=true \
steps=[filter] \
msin.starttime='29-Jun-2019/05:08:00.581' \
msin.ntimes=5 \
msin.nchan=5 \
filter.blrange="[0,100,0,100]"
"""


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar(f"{tcf.RESOURCEDIR}/{DISH_MSIN}.tgz")
    untar(f"{tcf.SRCDIR}/{MSAPPLYBEAM}.tgz")


@pytest.mark.parametrize("usechannelfreq", [False, True])
def test_with_invert_and_chanfreq(usechannelfreq):
    msout = "outinv.ms"

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={msout}",
            "steps=[applybeam]",
            f"applybeam.usechannelfreq={'true' if usechannelfreq else 'false'}",
            "applybeam.invert=true",
        ]
    )
    ms_column = "DATA_ucf" if usechannelfreq else "DATA_noucf"
    taql_command = f"select from {msout} t1, {MSAPPLYBEAM} t2 where not all(near(t1.DATA,t2.{ms_column},8e-5) || (isnan(t1.DATA) && isnan(t2.{ms_column})))"
    assert_taql(taql_command)

    # Test with invert=false on the resulting outinv.ms
    check_call(
        [
            tcf.DP3EXE,
            f"msin={msout}",
            f"msout=out.ms",
            "steps=[applybeam]",
            f"applybeam.usechannelfreq={'true' if usechannelfreq else 'false'}",
            "applybeam.invert=false",
        ]
    )
    taql_command = f"select from out.ms t1, {MSIN} t2 where not all(near(t1.DATA,t2.DATA,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA)))"
    assert_taql(taql_command)


def test_skip_stations():
    msout = "outinv.ms"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={msout}",
            "steps=[applybeam]",
            "applybeam.usechannelfreq=true",
            "applybeam.beammode=element",
            "applybeam.invert=true",
            "applybeam.skipstations=[RS208HBA,RS305HBA]",
        ]
    )

    # Check computed values (for included stations) against reference values for not skipped stations
    taql_command = f"""
        select from {msout} t1, {MSAPPLYBEAM} t2 where not (
            all(near(t1.DATA, t2.DATA_ELEMENT, 8e-5)) ||
            mscal.ant1name() in ["RS208HBA", "RS305HBA"] ||
            mscal.ant2name() in ["RS208HBA", "RS305HBA"]
        )
       """
    assert_taql(taql_command.replace("\n", " "))

    # Check that nothing changed for skipped stations
    taql_command = f"""
        select from {msout} t1, {MSIN} t2 where
            not(all(near(t1.DATA, t2.DATA, 8e-5))) &&
            mscal.ant1name() in ["RS208HBA", "RS305HBA"] &&
            mscal.ant2name() in ["RS208HBA", "RS305HBA"]
        """
    assert_taql(taql_command.replace("\n", " "))


@pytest.mark.parametrize("beammode", ["ARRAY_FACTOR", "ELEMENT"])
def test_beammodes(beammode):
    msout = "outinv.ms"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={msout}",
            "steps=[applybeam]",
            "applybeam.usechannelfreq=true",
            f"applybeam.beammode={beammode}",
            "applybeam.invert=true",
        ]
    )

    # Check computed values (for included stations) against reference values
    taql_command = f"select from {msout} t1, {MSAPPLYBEAM} t2 where not all(near(t1.DATA,t2.DATA_{beammode},8e-5) || (isnan(t1.DATA) && isnan(t2.DATA_{beammode})))"
    assert_taql(taql_command.replace("\n", " "))


def test_with_updateweights():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=.",
            f"msout.weightcolumn=NEW_WEIGHT_SPECTRUM",
            "steps=[applybeam]",
            f"applybeam.updateweights=true",
        ]
    )
    taql_command = f"select from {MSIN} where all(near(WEIGHT_SPECTRUM, NEW_WEIGHT_SPECTRUM))"
    assert_taql(taql_command)


from testconfig import TAQLEXE


def test_dish_beam():
    # Apply the beam with an offset of [0.01,-0.02] to the phase center
    check_call(
        [
            tcf.DP3EXE,
            f"msin={DISH_MSIN}",
            f"msout=beam_applied.ms",
            "msout.overwrite=true",
            "steps=[applybeam]",
            "applybeam.usechannelfreq=true",
            "applybeam.direction=[0.91848969rad,-0.50149271rad]",
        ]
    )

    # Assert that there are 25 elements which are zeroes before and after applying the beam
    zero_elements = f"select t1.DATA[0,0] from {DISH_MSIN} t1, beam_applied.ms t2 where abs(t1.DATA[0,0])==0.0 AND abs(t2.DATA[0,0])==0.0"
    assert_taql(zero_elements, 25)

    # Assert that all other elements have changed after applying the beam
    different_elements = f"select t1.DATA[0,0] from {DISH_MSIN} t1, beam_applied.ms t2 where (abs(t1.DATA[0,0]-t2.DATA[0,0])>1e-1) AND abs(t2.DATA[0,0])!=0.0"
    assert_taql(different_elements, 425)

    # Assert that the beam value is correct per each channel
    assert_taql(
        f"select t1.DATA[0,0]/t2.DATA[0,0] from {DISH_MSIN} t1, beam_applied.ms t2 where abs(t1.DATA[0,0] / t2.DATA[0,0] - 0.146382) > 1e-6 AND abs(t2.DATA[0,0])!=0.0"
    )
    assert_taql(
        f"select t1.DATA[1,0]/t2.DATA[1,0] from {DISH_MSIN} t1, beam_applied.ms t2 where abs(t1.DATA[1,0] / t2.DATA[1,0] - 0.146003) > 1e-6 AND abs(t2.DATA[1,0])!=0.0"
    )
    assert_taql(
        f"select t1.DATA[2,0]/t2.DATA[2,0] from {DISH_MSIN} t1, beam_applied.ms t2 where abs(t1.DATA[2,0] / t2.DATA[2,0] - 0.145628) > 1e-6 AND abs(t2.DATA[2,0])!=0.0"
    )
    assert_taql(
        f"select t1.DATA[3,0]/t2.DATA[3,0] from {DISH_MSIN} t1, beam_applied.ms t2 where abs(t1.DATA[3,0] / t2.DATA[3,0] - 0.145258) > 1e-6 AND abs(t2.DATA[3,0])!=0.0"
    )
    assert_taql(
        f"select t1.DATA[4,0]/t2.DATA[4,0] from {DISH_MSIN} t1, beam_applied.ms t2 where abs(t1.DATA[4,0] / t2.DATA[4,0] - 0.144933) > 1e-6 AND abs(t2.DATA[4,0])!=0.0"
    )
