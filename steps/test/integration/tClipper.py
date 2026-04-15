# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call

import pytest

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, get_taql_result, run_dp3, run_in_tmp_path, untar

"""
Tests for clipper (coarsely predict and clip bright sources).

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tClipper.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"
CLIPPER_MS = "clipper.MS"
TEST_MAX_AMPLITUDE = 1.5
CLIPPER_SKYMODEL = "clipper.skymodel"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


@pytest.fixture()
def create_model_data():
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={CLIPPER_MS}",
            # "msin.useflag=false",
            "msout.scalarflags=false",
            "steps=[]",
        ]
    )
    with open(f"{CLIPPER_SKYMODEL}", "w") as f:
        f.write(
            "FORMAT = Name, Type, Patch, Ra, Dec, I, Q, U, V, ReferenceFrequency='74e6', SpectralIndex='[]', MajorAxis, MinorAxis, Orientation\n"
        )
        f.write(", , CasAPoint, 00:00:00, +00.00.00\n")
        f.write(
            "CasAPoint_1, POINT, CasAPoint, 23:23:20.0, +58.48.26.0, 1e4, 0.0, 0.0, 0.0, 74.0e6, [-0.6]"
        )


# Clipper test data: time_step, frequency_step, n_expected_flags, n_expected_unflags
test_clipper_data = [
    (1, 1, 75, 0),
    (1, 2, 74, 2),
    (2, 1, 74, 2),
    (2, 2, 74, 3),
]


@pytest.mark.parametrize(
    "time_step, frequency_step, n_expected_flags, n_expected_unflags",
    test_clipper_data,
)
@pytest.mark.parametrize(
    "use_fast_predict",
    [
        False,
        pytest.param(
            True,
            marks=pytest.mark.skipif(
                tcf.USE_FAST_PREDICT != "ON",
                reason="FastPredict is not available",
            ),
        ),
    ],
)
def test_clipper(
    time_step,
    frequency_step,
    n_expected_flags,
    n_expected_unflags,
    use_fast_predict,
    create_model_data,
):
    taql_command = f"UPDATE {CLIPPER_MS} set FLAG=false"
    get_taql_result(taql_command)
    taql_command = f"select from {CLIPPER_MS} where any(FLAG)"
    assert_taql(taql_command, 0)

    check_call(
        [
            tcf.DP3EXE,
            f"msin={CLIPPER_MS}",
            "msout=.",
            "steps=[clipper]",
            f"clipper.sourcedb={CLIPPER_SKYMODEL}",
            "clipper.usebeammodel=true",
            "clipper.flagallcorrelations=false",
            f"clipper.amplmax={TEST_MAX_AMPLITUDE}",
            f"clipper.timestep={time_step}",
            f"clipper.freqstep={frequency_step}",
            f"clipper.usefastpredict={use_fast_predict}",
        ]
    )
    taql_command = f"select from {CLIPPER_MS} where any(FLAG)"
    assert_taql(taql_command, n_expected_flags)

    check_call(
        [
            tcf.DP3EXE,
            f"msin={CLIPPER_MS}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[predict]",
            f"predict.sourcedb={CLIPPER_SKYMODEL}",
            "predict.usebeammodel=true",
            f"predict.usefastpredict={use_fast_predict}",
        ]
    )

    # Check that the expected number of high amplitudes is flagged.
    # Not all high amplitudes are flagged if clipper used a timestep
    # or freqstep unqual to 1.
    taql_command = f"select from {CLIPPER_MS} where any(not FLAG && (MODEL_DATA > {TEST_MAX_AMPLITUDE}))"
    assert_taql(taql_command, n_expected_unflags)


@pytest.mark.parametrize(
    "use_fast_predict",
    [
        False,
        pytest.param(
            True,
            marks=pytest.mark.skipif(
                tcf.USE_FAST_PREDICT != "ON",
                reason="FastPredict is not available",
            ),
        ),
    ],
)
def test_clipper_with_filter(use_fast_predict, create_model_data):
    """
    Check whether the filter sub-step indeed only clips the selected baselines,
    using a minified version of clipper.MS with only three baselines.
    """
    taql_command = f"UPDATE {CLIPPER_MS} set FLAG=false"
    get_taql_result(taql_command)
    taql_command = f"select from {CLIPPER_MS} where any(FLAG)"
    assert_taql(taql_command, 0)

    # Generate a MS with only 3 baselines:
    #   CS001HBA0 - RS106HBA
    #   CS001HBA0 - RS208HBA
    #   RS106HBA  - RS208HBA
    mini_ms = "mini-clipper.ms"
    run_dp3(
        [
            tcf.DP3EXE,
            f"msin={CLIPPER_MS}",
            "msin.ntimes=1",
            "msin.nchan=2",
            f"msout={mini_ms}",
            "steps=[filter]",
            "filter.baseline=CS001HBA0,RS106HBA,RS208HBA&",
            "filter.remove=true",
        ]
    )
    taql_command = f"select from {mini_ms} where any(FLAG)"
    assert_taql(taql_command, 0)

    # Ensure that without selection all baselines are clipped
    clip_all_ms = "clip-all-baselines.ms"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={mini_ms}",
            f"msout={clip_all_ms}",
            "steps=[clipper]",
            f"clipper.sourcedb={CLIPPER_SKYMODEL}",
            f"clipper.usefastpredict={use_fast_predict}",
            "clipper.timestep=1",
            "clipper.freqstep=1",
        ]
    )
    taql_command = f"select from {clip_all_ms} where any(FLAG)"
    assert_taql(taql_command, 3)

    # Check that both baselines that include RS106 are clipped
    core_ms = "core.ms"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={mini_ms}",
            f"msout={core_ms}",
            "steps=[clipper]",
            f"clipper.sourcedb={CLIPPER_SKYMODEL}",
            f"clipper.usefastpredict={use_fast_predict}",
            "clipper.baseline=RS106*",
            "clipper.timestep=1",
            "clipper.freqstep=1",
        ]
    )
    taql_command = f"select from {core_ms} where any(FLAG)"
    assert_taql(taql_command, 2)

    # Check that only the 3rd baseline is clipped
    remote_ms = "remote.ms"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={mini_ms}",
            f"msout={remote_ms}",
            "steps=[clipper]",
            f"clipper.sourcedb={CLIPPER_SKYMODEL}",
            f"clipper.usefastpredict={use_fast_predict}",
            "clipper.baseline=RS*&",
            "clipper.timestep=1",
            "clipper.freqstep=1",
        ]
    )
    taql_command = f"select from {remote_ms} where any(FLAG)"
    assert_taql(taql_command, 1)
