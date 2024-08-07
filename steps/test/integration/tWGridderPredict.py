# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

from subprocess import check_call, check_output

import pytest

""" Append current directory to system path in order to import testconfig """
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, get_taql_result, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tWGridderPredict.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MODEL_IMAGES = "idg-fits-sources.tbz2"
MSINTGZ = "tDDECal.in_MS.tgz"
MSIN = "tDDECal.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSINTGZ}")
    untar(f"{tcf.DDECAL_RESOURCEDIR}/{MODEL_IMAGES}")


def test_polynomial_frequency_terms():
    """
    Three images are given to wgridderpredict, each represents a polynomial
    term on frequency. So the output should be
    term0 + term1 * freq + term2 * freq**2
    Since the test involves the center pixel only, taql can check the values:
    the visibilities should all have the value of the center pixel.
    """
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[wgridderpredict]",
            f"wgridderpredict.regions={tcf.DDECAL_RESOURCEDIR}/center-center.reg",
            "wgridderpredict.images=[term0-model.fits,term1-model.fits,term2-model.fits]",
            "wgridderpredict.sumfacets=true",
            "msout.datacolumn=TERMS_DATA",
        ]
    )

    # Get frequencies from taql (parsing output, by lack of python-casacore)
    # STR is necessary to avoid exponential notation with loss of precision
    taql_command = f"SELECT STR(CHAN_FREQ, 20.12) FROM {MSIN}::SPECTRAL_WINDOW"
    freqs_str = get_taql_result(taql_command).split("\n")[0]
    freqs = [float(freq) for freq in freqs_str.strip("[]").split(", ")]

    # Check each channel one by one, ignore autocorrelations
    for ch in range(len(freqs) - 1):
        # The factors 10, 20000 and 30000 match those in tIDGPredict_ref.py
        expected_value = (
            10
            + 20000 * (freqs[ch] / freqs[0] - 1)
            + 30000 * (freqs[ch] / freqs[0] - 1) ** 2
        )
        taql_command = f"""
            SELECT FROM {MSIN}
            WHERE
                ANTENNA1 != ANTENNA2
                AND NOT near(TERMS_DATA[{ch},0], {expected_value}, 2e-5)
        """.replace(
            "\n", " "
        )
        assert_taql(taql_command)
