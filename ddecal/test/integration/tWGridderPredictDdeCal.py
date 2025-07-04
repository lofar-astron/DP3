# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import sys

import h5py
import numpy as np
import pytest

""" Append current directory to system path in order to import testconfig """
sys.path.append(".")

import testconfig as tcf
from utils import run_dp3, run_in_tmp_path, untar

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
    run_dp3(
        [
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
    run_dp3(
        [
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


def test_wgridderpredict_with_ddecal_wildcards():
    """
    Use DDECal to check that wgridderpredict writes the correct visibilities
    into model data buffer when wildcards are used:
    - Predict three out of four model point sources into the data buffer
    - Use wgridderpredict with a fits image of those same three point sources,
      each in a facet, into three model data buffers.
    - Select a subset (3) of those initial four directions using wildcards.
    - Use DDECal to calibrate each of the directions. Since the model is equal
      to the data buffer, the solutions for every facet should be 1.0.
    """
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[predict,wgridderpredict, ddecal]",
            f"predict.sourcedb={tcf.DDECAL_RESOURCEDIR}/foursources.skymodel",
            "predict.sources=[source1, source2, source3]",
            f"wgridderpredict.regions={tcf.DDECAL_RESOURCEDIR}/foursources.reg",
            "wgridderpredict.images=[foursources-model.fits]",
            "ddecal.reusemodel=[wgridderpredict.s*]",
            "ddecal.mode=scalaramplitude",
            "ddecal.h5parm=instrument.h5",
            "ddecal.solint=0",
            "ddecal.nchan=0",
        ]
    )

    with h5py.File("instrument.h5", "r+") as f:
        values = f["sol000"]["amplitude000"]["val"]
        axes = values.attrs["AXES"].decode().split(",")
        index_of_direction_axis = axes.index("dir")

        # Verify that the direction axis contains 3 elements, equal to the selected model directions.
        assert values.shape[index_of_direction_axis] == len(
            ["source1", "source2", "source3"]
        )

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
    run_dp3(
        [
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

    run_dp3(
        [
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


def test_wgridderpredict_and_applycal_with_ddecal(make_h5parm):
    """
    Use DDECal to check that ApplyCal can apply the correction to the
    model data buffers.

    - Predict four point sources into the data buffer, applying solutions from
      an h5parm file.
    - Use wgridderpredict with a fits image of those same four point sources,
      each in a facet, into four model data buffers.
    - Apply the solutions from the h5parm file to model data buffers.
    - Use DDECal to calibrate each of the directions. Since the model data is
      predicted using the same parameters as the main data buffer
      the solutions for every facet should be 1.0.
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

    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[h5parmpredict, wgridderpredict, applycal, ddecal]",
            f"h5parmpredict.sourcedb={tcf.DDECAL_RESOURCEDIR}/foursources.skymodel",
            "h5parmpredict.applycal.parmdb=instrument.h5",
            "h5parmpredict.applycal.correction=amplitude000",
            f"wgridderpredict.regions={tcf.DDECAL_RESOURCEDIR}/foursources.reg",
            "wgridderpredict.images=[foursources-model.fits]",
            "applycal.parmdb=instrument.h5",
            "applycal.invert=false",
            "applycal.correction=amplitude000",
            "applycal.usemodeldata=true",
            "ddecal.reusemodel=[wgridderpredict.CygA, wgridderpredict.source1, wgridderpredict.source2, wgridderpredict.source3]",
            "ddecal.mode=scalaramplitude",
            "ddecal.h5parm=instrument_output.h5",
            "ddecal.solint=0",
            "ddecal.nchan=0",
        ]
    )

    input_file = h5py.File("instrument.h5", "r+")
    input_values = input_file["sol000"]["amplitude000"]["val"]

    output_file = h5py.File("instrument_output.h5", "r+")
    output_values = output_file["sol000"]["amplitude000"]["val"]

    assert input_values.shape == output_values.shape

    # Check output values
    # If all data for an antenna is flagged the solution will be NaN.
    # Because data and model data are predicted applying the same h5parm file
    # all valid solutions should be 1.0.
    count = 0
    for y in np.nditer(output_values):
        if not np.isnan(y):
            assert np.isclose(y, 1.0, rtol=0.0, atol=1e-3)
            count += 1

    # One out of eight antennas has 100% of its data flagged.
    # There are four directions, so there should be 7*4=28 valid solutions.
    assert count == 28
