# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
from subprocess import check_call, check_output

import h5py
import numpy as np
import pytest

""" Append current directory to system path in order to import testconfig """
import sys

sys.path.append(".")

import testconfig as tcf
from utils import assert_taql, run_in_tmp_path, untar

"""
Script can be invoked in two ways:
- as standalone from the build/ddecal/test/integration directory,
  using `pytest source/tWGridderPredictDdeCal.py` (extended with pytest options of your choice)
- using ctest, see DP3/ddecal/test/integration/CMakeLists.txt
"""

MODEL_IMAGES = "idg-fits-sources.tbz2"
MSINTGZ = "tDDECal.in_MS.tgz"
MSIN = "tDDECal.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSINTGZ}")
    untar(f"{tcf.DDECAL_RESOURCEDIR}/{MODEL_IMAGES}")


@pytest.fixture()
def make_h5parm():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "verbosity=quiet",
            f"msin={MSIN}",
            "steps=[ddecal,null]",
            f"ddecal.sourcedb={tcf.DDECAL_RESOURCEDIR}/foursources.skymodel",
            "ddecal.h5parm=instrument.h5",
            "ddecal.maxiter=1",
            "ddecal.mode=scalaramplitude",
            "ddecal.solint=0",
            "ddecal.nchan=0",
        ]
    )


def test_wgridderpredict_with_ddecal():
    """
    Use DDECal to check that wgridderpredict writes the correct visibilities
    into model data buffer:
    - Predict four point sources into the data buffer
    - Use wgridderpredict with a fits image of those same four point sources,
      each in a facet, into four model data buffers
    - Use DDECal to calibrate each of the directions. Since the model is equal
      to the data buffer, the solutions for every facet should be 1.0.
    """
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[predict,wgridderpredict, ddecal]",
            f"predict.sourcedb={tcf.DDECAL_RESOURCEDIR}/foursources.skymodel",
            f"wgridderpredict.regions={tcf.DDECAL_RESOURCEDIR}/foursources.reg",
            "wgridderpredict.images=[foursources-model.fits]",
            "ddecal.reusemodel=[wgridderpredict.CygA, wgridderpredict.source1, wgridderpredict.source2, wgridderpredict.source3]",
            "ddecal.mode=scalaramplitude",
            "ddecal.h5parm=instrument.h5",
            "ddecal.solint=0",
            "ddecal.nchan=0",
        ]
    )

    with h5py.File("instrument.h5", "r+") as f:
        values = f["sol000"]["amplitude000"]["val"]

        # Check solutions are 1.0
        assert np.allclose(
            values[np.isfinite(values)], 1.0, rtol=0.0, atol=1.0e-4
        )


def test_wgridderpredict_and_applybeam_with_ddecal():
    """
    Use DDECal to check that ApplyBeam can apply the beam to the
    model data buffers:
    - Predict four point sources into the data buffer, using the beam model
    - Use wgridderpredict with a fits image of those same four point sources,
      each in a facet, into four model data buffers
    - Apply the beam to model data buffers
    - Use DDECal to calibrate each of the directions. Since the model is equal
      to the data buffer, the solutions for every facet should be 1.0.
    """
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[predict,wgridderpredict, applybeam, ddecal]",
            f"predict.sourcedb={tcf.DDECAL_RESOURCEDIR}/foursources.skymodel",
            "predict.usebeammodel=True",
            f"wgridderpredict.regions={tcf.DDECAL_RESOURCEDIR}/foursources.reg",
            "wgridderpredict.images=[foursources-model.fits]",
            "applybeam.usemodeldata=true",
            "applybeam.invert=false",
            "ddecal.reusemodel=[wgridderpredict.CygA, wgridderpredict.source1, wgridderpredict.source2, wgridderpredict.source3]",
            "ddecal.mode=scalaramplitude",
            "ddecal.h5parm=instrument.h5",
            "ddecal.solint=0",
            "ddecal.nchan=0",
        ]
    )

    with h5py.File("instrument.h5", "r+") as f:
        values = f["sol000"]["amplitude000"]["val"]

        # Check all solutions are 1.0
        assert np.allclose(
            values[np.isfinite(values)], 1.0, rtol=0.0, atol=1.0e-4
        )


def test_wgridderpredict_with_ddecal_h5parmvalues(make_h5parm):
    """
    Same as test_wgridderpredict_with_ddecal, but with different calibration
    solutions for each facet. The corrupted data is made with h5parmpredict:
    a different corruption is applied for each point source (direction).
    """
    with h5py.File("instrument.h5", "r+") as f:
        values = f["sol000"]["amplitude000"]["val"]

        n_antennas = len(f["sol000"]["amplitude000"]["ant"])
        n_directions = len(f["sol000"]["amplitude000"]["dir"])

        # Set solutions to new value
        values[0, 0, :, :] = (
            1.0
            + np.arange(n_directions).reshape(1, -1)
            + 1.0e-1 * np.arange(n_antennas).reshape(-1, 1)
        )
        new_values = values[:]

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[h5parmpredict,wgridderpredict,ddecal]",
            f"h5parmpredict.sourcedb={tcf.DDECAL_RESOURCEDIR}/foursources.skymodel",
            "h5parmpredict.applycal.parmdb=instrument.h5",
            "h5parmpredict.applycal.correction=amplitude000",
            f"wgridderpredict.regions={tcf.DDECAL_RESOURCEDIR}/foursources.reg",
            "wgridderpredict.images=[foursources-model.fits]",
            "ddecal.reusemodel=[wgridderpredict.CygA, wgridderpredict.source1, wgridderpredict.source2, wgridderpredict.source3]",
            "ddecal.mode=scalaramplitude",
            "ddecal.h5parm=instrument_output.h5",
            "ddecal.solint=0",
            "ddecal.nchan=0",
            "ddecal.maxiter=200",
            "ddecal.tolerance=1.0e-8",
        ]
    )

    with h5py.File("instrument_output.h5", "r+") as f:
        values = f["sol000"]["amplitude000"]["val"]

        assert np.allclose(
            values[np.isfinite(values)],
            new_values[np.isfinite(values)],
            rtol=0.0,
            atol=5.0e-3,
        )
