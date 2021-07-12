# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import shutil
import uuid
from subprocess import check_call, check_output, run

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from testconfig import TAQLEXE
from utils import assert_taql, untar_ms

"""
Script can be invoked in two ways:
- as standalone from the build/ddecal/test/integration directory,
  using `pytest source/tBdaDdeCal.py` (extended with pytest options of your choice)
- using ctest, see DP3/ddecal/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-bda.MS"
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
def create_skymodel():
    with open("test.skymodel","w") as f:
        f.write("FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\r\n")
        f.write("center, POINT, 16:38:28.205000, + 63.44.34.314000, 10, , , , , \r\n")
        f.write("ra_off, POINT, 16:38:28.205000, + 64.44.34.314000, 10, , , , , \r\n")
        f.write("radec_off, POINT, 16:38:28.205000, +65.44.34.314000, 10, , , , , \r\n")

    check_call(
            [
                tcf.MAKESOURCEDBEXE, 
                "in=test.skymodel", 
                "out=test.sourcedb"
            ]
        )


def test_only_predict(create_skymodel):
    
    check_call(
            [
                tcf.DP3EXE,
                f"msin={MSIN}",
                "msout=BDADDECal_onlypredict.MS",
                "steps=[ddecal]",
                "ddecal.onlypredict=true",
                "ddecal.directions=[[center, dec_off],[ra_off],[radec_off]]",
                "ddecal.sourcedb=test.sourcedb",
                "msout.overwrite=true", 
                "bda=true"
            ]
        )

    check_call(
            [
                tcf.DP3EXE,
                f"msin={MSIN}",
                "msout=PREDICT_DIR_1.MS",
                "steps=[predict]",
                "predict.sources=[center, dec_off]",
                "predict.sourcedb=test.sourcedb",
                "msout.overwrite=true", 
                "bda=true"
            ]
        )

    check_call(
            [
                tcf.DP3EXE,
                f"msin={MSIN}",
                "msout=PREDICT_DIR_2.MS",
                "steps=[predict]",
                "predict.sources=[ra_off]",
                "predict.sourcedb=test.sourcedb",
                "msout.overwrite=true", 
                "bda=true"
            ]
        )


    check_call(
            [
                tcf.DP3EXE,
                f"msin={MSIN}",
                "msout=PREDICT_DIR_3.MS",
                "steps=[predict]",
                "predict.sources=[radec_off]",
                "predict.sourcedb=test.sourcedb",
                "msout.overwrite=true", 
                "bda=true"
            ]
        )

    taql_command = "select from(select gsumsqr(abs(t_all.DATA[isnan(t_all.DATA)]-(t_1.DATA[isnan(t_1.DATA)]+t_2.DATA[isnan(t_2.DATA)]+t_3.DATA[isnan(t_3.DATA)]))) as diff from BDADDECal_onlypredict.MS t_all, PREDICT_DIR_1.MS t_1,  PREDICT_DIR_2.MS t_2,  PREDICT_DIR_3.MS t_3) where diff>1.e-6"
    assert_taql(taql_command)
