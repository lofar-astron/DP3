# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later……

import glob
import math
import os
import sys
from subprocess import check_call

import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import get_taql_result, run_dp3, run_in_tmp_path, untar

"""
Integration tests for the ImageStep step.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tImageStep.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""


MSIN = "tNDPPP-generic.MS"
MWA_MS = f"{tcf.BINDIR}/test_data/MWA-single-timeslot.ms"
MOCK_SKYMODEL = "mock.skymodel"
FITS_PREFIX = "imagestep"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    if tcf.HAVE_WSCLEAN != "TRUE":
        pytest.skip(reason="WSClean is not available")

    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


@pytest.fixture()
def mock_sky_model():
    """
    Write a sky model (mock-model.txt) with a single off-centre point source.
    """
    mock_model = """FORMAT = Name, Type, Patch, Ra, Dec, I, SpectralIndex='[]', LogarithmicSI, ReferenceFrequency, MajorAxis, MinorAxis, Orientation

, , point_source, 08:10:00,-42.00.00
point, POINT, point_source, 08:10:00, -42.00.00, 100.0, [], false, 150000000, 0.0, 0.0, 0.0
"""

    with open(MOCK_SKYMODEL, "w") as file:
        file.write(mock_model)


def format_wsclean_fits_file_name(prefix, image_type="image"):
    return f"{prefix}-{image_type}.fits"


def test_image_step():
    """
    Test whether the step creates the expected dirty image + image when
    the cadence encompases the entire MS.
    """
    wsclean_options = f"-name {FITS_PREFIX} -size 1024 1024 -scale 1arcmin"
    run_dp3(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=.",
            "steps=[image]",
            "image.cadence=60",
            f"image.options='{wsclean_options}'",
        ]
    )

    n_fits_images = len(
        glob.glob(format_wsclean_fits_file_name("*", image_type="*"))
    )
    assert n_fits_images == 2

    assert os.path.exists(
        format_wsclean_fits_file_name(FITS_PREFIX, image_type="dirty")
    )
    assert os.path.exists(
        format_wsclean_fits_file_name(FITS_PREFIX, image_type="image")
    )


def test_image_step_with_multiple_images():
    """
    Test whether the imager outputs new images (dirty, model, residual, image) for each
    time slot when the cadence matches the time resolution (10.0139).
    """
    wsclean_options = f"-name {FITS_PREFIX} -size 1024 1024 -scale 1arcmin -niter 1 -mgain 0.8"
    run_dp3(
        [
            f"msin={MSIN}",
            f"msout=.",
            "steps=[image]",
            "image.cadence=11",
            f"image.fitsprefix={FITS_PREFIX}",
            f"image.options='{wsclean_options}'",
        ]
    )

    n_expected_images = 6

    n_dirty_images = len(
        glob.glob(format_wsclean_fits_file_name("*", image_type="dirty"))
    )
    assert n_dirty_images == n_expected_images

    n_model_images = len(
        glob.glob(format_wsclean_fits_file_name("*", image_type="model"))
    )
    assert n_model_images == n_expected_images

    n_residual_images = len(
        glob.glob(format_wsclean_fits_file_name("*", image_type="residual"))
    )
    assert n_residual_images == n_expected_images

    n_images = len(
        glob.glob(format_wsclean_fits_file_name("*", image_type="image"))
    )
    assert n_images == n_expected_images


def test_point_source(mock_sky_model):
    """
    Verify that the imager correctly creates an image with the point source on the expected pixel.
    """
    wsclean_options = f"-name {FITS_PREFIX} -size 512 512 -scale 1arcmin -niter 666 -mgain 0.8 -auto-threshold 3"
    run_dp3(
        [
            f"msin={MWA_MS}",
            "msin.useflag=False",
            f"msout=.",
            "steps=[predict,image]",
            f"predict.sourcedb={MOCK_SKYMODEL}",
            "predict.sources=[point_source]",
            "predict.operation=replace",
            "image.cadence=5",
            f"image.options='{wsclean_options}'",
        ]
    )

    fits_image = format_wsclean_fits_file_name(FITS_PREFIX, image_type="image")
    assert os.path.exists(fits_image)

    with fits.open(fits_image) as hdu_list:
        assert len(hdu_list) == 1

        wcs = WCS(hdu_list[0].header)
        source_position = SkyCoord("08h10m00s -42d00m00s")
        source_x, source_y = wcs[0, 0, :, :].world_to_pixel(source_position)

        source_x = int(source_x)
        source_y = int(source_y)

        image_data = hdu_list[0].data[0, 0, :, :]

        # Verify the position of the point source and ensure that the intensity
        # is not far off (< 10%).
        assert math.isclose(image_data[source_y, source_x], 100, rel_tol=0.1)
        assert math.isclose(image_data[source_y, source_x], image_data.max())
