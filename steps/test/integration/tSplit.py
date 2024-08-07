# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call

import pytest

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tSplit.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSAPPLYBEAM = "tApplyBeam.tab"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")
    untar(f"{tcf.SRCDIR}/{MSAPPLYBEAM}.tgz")


def test_split():
    msout_list = ["splitout1.ms", "splitout2.ms"]
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[split]",
            "split.steps=[applybeam,out]",
            "split.replaceparms=[out.name,applybeam.usechannelfreq]",
            f"out.name=[{msout_list[0]},{msout_list[1]}]",
            f"applybeam.usechannelfreq=[false,true]",
            "applybeam.invert=true",
        ]
    )

    for msout in msout_list:
        data_column = "DATA_noucf" if msout == "splitout1.ms" else "DATA_ucf"
        taql_command = f"select from {msout} t1, {MSAPPLYBEAM} t2 where not all(near(t1.DATA,t2.{data_column},8e-5) || (isnan(t1.DATA) && isnan(t2.{data_column})))"
        assert_taql(taql_command)
