# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
from subprocess import check_call, check_output
from envbash import load_envbash

"""
Replacement for tPredict.sh using pytest.

Script can be invoked in two ways:
- as standalone from the build/steps/integration directory,
  using `pytest source/tPredict.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
MSPREDICT = "tPredict.tab"


def assert_taql(taqlexe, command):
    result = check_output([taqlexe, "-noph", command]).decode().strip()
    assert result == "select result of 0 rows"


def untar_ms(source, outdir):
    if not os.path.isfile(source):
        raise IOError(f"Not able to find {source} containing the reference solutions.")

    # Untar if needed
    check_call(
        ["tar", "xf", source, "-C", outdir, "--skip-old-files",]
    )


@pytest.fixture(autouse=True)
def source_env():
    if not os.path.isfile("testInit.sh"):
        raise IOError(
            "Not able to find testInit.sh file. This file should be located in build/steps/test/integration."
        )
    load_envbash("testInit.sh")

    # Untar tNDPPP-generic.MS.tgz in binary dir if needed
    untar_ms(f"{os.environ['resourcedir']}/{MSIN}.tgz", os.environ["bindir"])

    # Untar tPredict.tab.tgz in binary dir if needed
    untar_ms(f"{os.environ['srcdir']}/{MSPREDICT}.tgz", os.environ["bindir"])


def test_with_beam_subtract():
    msin = os.path.join(os.environ["bindir"], MSIN)
    mspredict = os.path.join(os.environ["bindir"], MSPREDICT)

    check_call(
        [
            os.environ["dp3exe"],
            f"msin={msin}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={msin}/sky",
            "predict.usebeammodel=true",
            "predict.operation=subtract",
        ]
    )

    # Compare the MODEL_DATA column of the output MS with the original data minus the BBS reference output.
    taql_command = f"select t1.MODEL_DATA, t1.DATA-t2.PREDICT_beam, abs(t1.MODEL_DATA/(t1.DATA-t2.PREDICT_beam)-1) from {msin} t1, {mspredict} t2 where not all(near(t1.MODEL_DATA,t1.DATA-t2.PREDICT_beam,5e-2) || (isnan(t1.DATA) && isnan(t2.PREDICT_beam)) || nearabs(t2.PREDICT_beam, 0, 1e-5))"
    assert_taql(os.environ["taqlexe"], taql_command)


def test_without_beam_add():
    msin = os.path.join(os.environ["bindir"], MSIN)
    mspredict = os.path.join(os.environ["bindir"], MSPREDICT)

    check_call(
        [
            os.environ["dp3exe"],
            f"msin={msin}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={msin}/sky",
            "predict.usebeammodel=false",
            "predict.operation=add",
        ]
    )

    # Compare the MODEL_DATA column of the output MS with the original data plus the BBS reference output.
    taql_command = f"select from {msin} t1, {mspredict} t2 where not all(near(t1.MODEL_DATA,t1.DATA+t2.PREDICT_nobeam,5e-2) || (isnan(t1.DATA) && isnan(t2.PREDICT_nobeam)))"
    assert_taql(os.environ["taqlexe"], taql_command)


@pytest.mark.parametrize("use_beam", [False, True])
def test_without_and_without_beam(use_beam):
    msin = os.path.join(os.environ["bindir"], MSIN)
    mspredict = os.path.join(os.environ["bindir"], MSPREDICT)
    predict_column = "PREDICT_beam" if use_beam else "PREDICT_nobeam"

    check_call(
        [
            os.environ["dp3exe"],
            f"msin={msin}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={msin}/sky",
            f"predict.usebeammodel={'true' if use_beam else 'false'}",
        ]
    )
    taql_command = f"select from {msin} t1, {mspredict} t2 where not all(near(t1.MODEL_DATA,t2.{predict_column},5e-2) || (isnan(t1.MODEL_DATA) && isnan(t2.{predict_column})))"
    assert_taql(os.environ["taqlexe"], taql_command)
