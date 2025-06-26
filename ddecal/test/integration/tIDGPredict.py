# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

from subprocess import check_call, check_output

import pytest

""" Append current directory to system path in order to import testconfig """
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/ddecal/test/integration directory,
  using `pytest source/tIDGPredict.py` (extended with pytest options of your choice)
- using ctest, see DP3/ddecal/test/integration/CMakeLists.txt
"""

REF_SOLUTIONS = "idg-fits-sources.tbz2"
MSINTGZ = "tDDECal.in_MS.tgz"
MSIN = "tDDECal.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSINTGZ}")
    untar(f"{tcf.DDECAL_RESOURCEDIR}/{REF_SOLUTIONS}")


def compare_results(source_name):
    """
    Ignore baselines with antennna 6 for now, since not all predictors
    generate visibilities for those baselines.
    TODO (AST-223): Investigate what the expected behavior is.
    """

    taql_command = f"select from {MSIN} where ANTENNA2 != 6 and not all(near({source_name}_DATA,MODEL_DATA,1e-3))"
    assert_taql(taql_command)


def test_input_with_four_sources():
    """
    Test an input with four sources.
    Since DP3 does not depend on wsclean, tDDECal.MS has a foursources_DATA
    column with the reference output/visibilities for foursources-model.fits.

    These commands originally generated the column:
    wsclean -use-idg -predict -name resources/foursources tDDECal.MS
    taql "alter table tDDECal.MS rename column MODEL_DATA to foursources_DATA"

    20250626: Nowadays, wsclean should write directly to the foursources_DATA
    column. Otherwise, the DP3 call below fails since MODEL_DATA_dm already
    exists in tDDECal.MS. The command would be:
    wsclean -use-idg -predict -name resources/foursources tDDECal.MS
    -model-column foursources_DATA
    """

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.idg.regions={tcf.DDECAL_RESOURCEDIR}/foursources.reg",
            "ddecal.idg.images=[foursources-model.fits]",
            "ddecal.onlypredict=True",
            "msout.datacolumn=MODEL_DATA",
            "ddecal.modelnextsteps.CygA=[uvwflag]",
        ]
    )
    compare_results("foursources")


@pytest.mark.skip(
    "This test is kept for debugging IDGPredict. It also uses WSClean, which "
    "is normally not required for running tests."
)
@pytest.mark.parametrize("source", ["center", "ra", "dec", "radec"])
@pytest.mark.parametrize("offset", ["center", "dl", "dm", "dldm"])
def test_input_with_single_sources(source, offset):
    """Test inputs that contain a single source."""

    # Many {source}-{offset}.reg files only include a "TODO" line.
    # Only use the files that contain a polygon.
    if (
        "polygon"
        not in open(f"{tcf.DDECAL_RESOURCEDIR}/{source}-{offset}.reg").read()
    ):
        pytest.skip(f"Polygon not found in {source}-{offset}.reg")

    # Generate reference predictions in the ${source}_DATA column.
    check_call(
        [
            "wsclean",
            "-use-idg",
            "-predict",
            "-name",
            source,
            "-model-column",
            f"{source}_DATA",
            MSIN,
        ]
    )

    # Predict source: $source offset: $offset using IDG
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.idg.regions={tcf.DDECAL_RESOURCEDIR}/{source}-{offset}.reg",
            f"ddecal.idg.images=[{source}-model.fits]",
            "ddecal.onlypredict=True",
            "msout.datacolumn=MODEL_DATA",
        ]
    )

    compare_results(source)


def test_result():
    """Test if IDGPredict step will have the same results as DDECal"""
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[idgpredict]",
            f"idgpredict.regions={tcf.DDECAL_RESOURCEDIR}/foursources.reg",
            "idgpredict.images=[foursources-model.fits]",
            "msout.datacolumn=MODEL_DATA",
        ]
    )
    compare_results("foursources")


def test_multiple_data_sources():
    """Test multiple data sources for DDECal"""

    common_args = [
        "checkparset=1",
        f"msin={MSIN}",
        "msout=.",
        "numthreads=1",
        "steps=[ddecal]",
        "ddecal.idg.images=[foursources-model.fits]",
        "ddecal.onlypredict=True",
        "msout.datacolumn=MODEL_DATA",
    ]

    # Create model data column with 3 sources
    check_call(
        [
            tcf.DP3EXE,
            f"ddecal.idg.regions={tcf.DDECAL_RESOURCEDIR}/threesources.reg",
        ]
        + common_args
    )

    # Run DDECal with 3 directions in the MODEL_DATA column and 1 direction using IDG
    check_call(
        [
            tcf.DP3EXE,
            f"ddecal.idg.regions={tcf.DDECAL_RESOURCEDIR}/onesource.reg",
            "ddecal.modeldatacolumns=[MODEL_DATA]",
        ]
        + common_args
    )

    # Results of a the DDECal above (3 directions modeldata and 1 direction IDG)
    # should be equal to a run with foursources (4 direction IDG).
    compare_results("foursources")


def test_polynomial_frequency_term_corrections():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.idg.regions={tcf.DDECAL_RESOURCEDIR}/center-center.reg",
            "ddecal.idg.images=[term0-model.fits,term1-model.fits,term2-model.fits]",
            "ddecal.onlypredict=True",
            "msout.datacolumn=TERMS_DATA",
        ]
    )

    # Since the test involves the center pixel only, taql can check the values.
    ch_count = int(
        check_output(
            [
                tcf.TAQLEXE,
                "-nopr",
                "-noph",
                f"select count(CHAN_FREQ) from {MSIN}::SPECTRAL_WINDOW",
            ]
        )
    )

    for ch in range(ch_count - 1):
        # The factors 10, 20000 and 30000 match those in tIDGPredict_ref.py
        taql_command = f"select from {MSIN} where not(TERMS_DATA[{ch},0]=0 or near(TERMS_DATA[{ch},0], (select 10+20000*(CHAN_FREQ[{ch}]/CHAN_FREQ[0]-1) +30000*(CHAN_FREQ[{ch}]/CHAN_FREQ[0]-1)**2 from ::SPECTRAL_WINDOW)[0], 1e-3))"
        assert_taql(taql_command)
