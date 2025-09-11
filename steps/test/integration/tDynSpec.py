# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later……

import os
import sys
from subprocess import check_call

import casacore.tables
import numpy as np
import pytest
from astropy.io import fits

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import get_taql_result, run_dp3, run_in_tmp_path, untar

"""
Integration tests for the DynSpec step.

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tDynSpec.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""


MSIN = "tNDPPP-generic.MS"
MSOUT = "out.ms"
SKYMODEL = f"{tcf.RESOURCEDIR}/tNDPPP-generic-skymodel.txt"
MOCK_SKYMODEL = "mock.skymodel"
MOCK_SOURCE_LIST = "mock-sources.txt"
MOCK_SOURCES = ["point_source", "constant_source"]
MOCK_H5PARM = "mock_solutions.h5"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")

    # Unset all flags
    get_taql_result(f"update {MSIN} set FLAG=False, FLAG_ROW=False")


@pytest.fixture()
def mock_sky_model():
    """
    Write a sky model (mock-model.txt) with 1 point source, 1 small Gaussian source and 1 larger Gaussian source.
    """
    mock_model = """FORMAT = Name, Type, Patch, Ra, Dec, I, SpectralIndex='[]', LogarithmicSI, ReferenceFrequency, MajorAxis, MinorAxis, Orientation

, , point_source, 02:00:00, -30.00.00
point, POINT, point_source, 02:00:00, -30.00.00, 10.0, [20.0], true, 134475000, 0.0, 0.0, 0.0

, , constant_source, 00:10:00, +30.50.00
point, POINT, constant_source, 00:10:00, +30.50.00, 100.0, [0.0], true, 134475000, 0.0, 0.0, 0.0
"""

    with open(MOCK_SKYMODEL, "w") as file:
        file.write(mock_model)


@pytest.fixture()
def mock_source_list():
    """
    Write a sky model (mock-model.txt) with 1 point source, 1 small Gaussian source and 1 larger Gaussian source.
    """
    mock_list = f"""FORMAT = Name, Ra, Dec
{MOCK_SOURCES[0]}, 02:00:00, -30.00.00
{MOCK_SOURCES[1]}, 00:10:00, +30.50.00
"""

    with open(MOCK_SOURCE_LIST, "w") as file:
        file.write(mock_list)


def format_dynspec_fits_name(source_name):
    return f"{source_name}-dynspec.fits"


def test_dynspec(mock_source_list):
    """
    Check whether the dynamic spectra are written to disk.
    """
    run_dp3(
        [
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[dynspec]",
            f"dynspec.sourcelist={MOCK_SOURCE_LIST}",
        ]
    )

    for source_name in MOCK_SOURCES:
        assert os.path.exists(format_dynspec_fits_name(source_name))


@pytest.mark.parametrize("with_beam", [True, False])
def test_dynspec_with_predicted_source(
    mock_source_list, mock_sky_model, with_beam
):
    """
    Predict a single off-centre point source, phase shift towards it
    and record it's dynamic spectrum. (Turn off source subtraction.)
    """
    run_dp3(
        [
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[predict, dynspec]",
            f"predict.sourcedb={MOCK_SKYMODEL}",
            "predict.sources=[point_source]",
            f"predict.usebeammodel={with_beam}",
            f"dynspec.sourcelist={MOCK_SOURCE_LIST}",
            f"dynspec.beamcorrection={with_beam}",
            f"dynspec.applybeam.updateweights=True" if with_beam else "",
        ]
    )

    dynspec_fits = format_dynspec_fits_name(MOCK_SOURCES[0])
    with fits.open(dynspec_fits) as hdu_list:
        stokes_i_spectrum = hdu_list[0].data[0]

        # The predicted source is time-invariable, hence the
        # spectrum should be constant over time. We verify this
        # for each channel.
        for channel_data in stokes_i_spectrum:
            assert np.allclose(channel_data, channel_data[0])

        # Predicted source has a steep spectral index, such that the
        # flux scales (exponentially) with frequency. Hence, check
        # whether the gradient is positive for each time slot.
        for time_chunk in stokes_i_spectrum.T:
            assert np.all(np.diff(time_chunk) > 0)


@pytest.mark.parametrize("with_beam", [True, False])
@pytest.mark.parametrize("with_cal", [True, False])
def test_dynspec_with_source_subtraction(
    mock_source_list, mock_sky_model, with_beam, with_cal
):
    """
    Generate solutions with DDECal and subsequently generate model visibilties with
    the H5ParmPredict substep, predict an off-centre point source, phase shift towards
    it and record a dynamic spectrum.
    """
    h5parm_name = "instrument.h5"
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={SKYMODEL}",
            f"ddecal.h5parm={h5parm_name}",
            "ddecal.mode=fulljones",
            f"ddecal.usebeammodel={with_beam}",
        ]
    )

    target_source = MOCK_SOURCES[1]
    run_dp3(
        [
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[predict, dynspec]",
            f"predict.sourcedb={MOCK_SKYMODEL}",
            f"predict.sources=[{target_source}]",
            "predict.operation=add",
            f"predict.usebeammodel={with_beam}",
            f"dynspec.h5parmpredict.sourcedb={SKYMODEL}",
            f"dynspec.h5parmpredict.applycal.parmdb={h5parm_name}",
            "dynspec.h5parmpredict.applycal.correction=fulljones",
            f"dynspec.h5parmpredict.usebeammodel={with_beam}",
            f"dynspec.sourcelist={MOCK_SOURCE_LIST}",
            f"dynspec.applycal.parmdb={h5parm_name}" if with_cal else "",
            "dynspec.applycal.correction=fulljones" if with_cal else "",
            f"dynspec.beamcorrection={with_beam}",
            f"dynspec.applybeam.updateweights=True" if with_beam else "",
        ]
    )

    dynspec_fits = format_dynspec_fits_name(target_source)
    with fits.open(dynspec_fits) as hdu_list:
        stokes_i_spectrum = hdu_list[0].data[0]

        # The predicted source has a flat spectrum (zero spectral index),
        # hence the spectrum should be roughly zero. Has a relaxed constraint,
        # since we only subtract model visibilities.
        assert stokes_i_spectrum.std() < 0.1


@pytest.mark.parametrize("with_beam", [True, False])
def test_dynspec_with_subtraction_from_model_column(
    mock_source_list,
    mock_sky_model,
    with_beam,
):
    """
    Check whether model column based source subtraction works by subtracting
    data stored in the model column. Then predict an off-centre point source,
    phase shift towards it and record it's dynamic spectrum.
    """
    model_column = "MODEL_DATA"
    run_dp3(
        [
            f"msin={MSIN}",
            "msout=.",
            f"msout.datacolumn={model_column}",
            "steps=[]",
        ]
    )

    target_source = MOCK_SOURCES[1]
    run_dp3(
        [
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[predict, dynspec]",
            f"predict.sourcedb={MOCK_SKYMODEL}",
            f"predict.sources=[{target_source}]",
            "predict.operation=add",
            f"predict.usebeammodel={with_beam}",
            f"dynspec.sourcelist={MOCK_SOURCE_LIST}",
            f"dynspec.subtractmodelcolumn={model_column}",
            f"dynspec.beamcorrection={with_beam}",
            f"dynspec.applybeam.updateweights=True" if with_beam else "",
        ]
    )

    dynspec_fits = format_dynspec_fits_name(target_source)
    with fits.open(dynspec_fits) as hdu_list:
        stokes_i_spectrum = hdu_list[0].data[0]

        # The predicted source has a flat spectrum (zero spectral index),
        # hence the spread in the spectrum should be roughly zero. The constraint
        # can be more stringent than 'test_dynspec_with_source_subtraction',
        # since we subtract a copy of the data from the data, instead of predicted
        # visibilities.
        assert stokes_i_spectrum.std() < 1e-5


def test_dynspec_with_flags(mock_source_list, mock_sky_model):
    """
    Check whether a flagged timeslot is indeed empty.
    """
    target_source = MOCK_SOURCES[1]
    run_dp3(
        [
            f"msin={MSIN}",
            f"msout={MSOUT}",
            "steps=[predict, preflagger, dynspec]",
            f"predict.sourcedb={MOCK_SKYMODEL}",
            f"predict.sources=[{target_source}]",
            "preflagger.timeslot=[0]",
            f"dynspec.sourcelist={MOCK_SOURCE_LIST}",
        ]
    )

    dynspec_fits = format_dynspec_fits_name(target_source)
    with fits.open(dynspec_fits) as hdu_list:
        stokes_i_spectrum = hdu_list[0].data[0]

        test_tolerance = 1e-5
        assert np.allclose(stokes_i_spectrum[:, 0], 0, atol=test_tolerance)
        assert not np.allclose(
            stokes_i_spectrum[:, 1:], 0, atol=test_tolerance
        )
