# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Append current directory to system path in order to import testconfig
import sys
import numpy as np
import pytest

sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass


def test_metadata():
    pybuffer = dp3.DPBuffer()

    # Test that fields are empty when using the default consructor

    assert pybuffer.get_time() == 0
    assert pybuffer.get_exposure() == 0

    assert np.array(pybuffer.get_data(), copy=False).size == 0
    assert np.array(pybuffer.get_flags(), copy=False).size == 0
    assert np.array(pybuffer.get_weights(), copy=False).size == 0
    assert np.array(pybuffer.get_uvw(), copy=False).size == 0

    pybuffer.set_exposure(55)
    assert pybuffer.get_exposure() == 55

    pybuffer.set_time(10)
    assert pybuffer.get_time() == 10


def test_constructor():
    pybuffer = dp3.DPBuffer(10, 20)

    assert pybuffer.get_time() == 10
    assert pybuffer.get_exposure() == 20


def test_weights():
    pybuffer = dp3.DPBuffer()

    weights_invalid_shape = weights = np.array(
        [[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]], np.float16
    )

    with pytest.raises(RuntimeError):
        pybuffer.set_weights(weights_invalid_shape)

    weights = np.array(
        [
            [[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]],
            [[0.07, 0.08, 0.09], [0.10, 0.11, 0.12]],
        ],
        np.float16,
    )

    pybuffer.set_weights(weights)
    assert (np.array(pybuffer.get_weights(), copy=False) == weights).all()


def test_flags():
    pybuffer = dp3.DPBuffer()

    flags_invalid_shape = np.array(
        [[True, True, True], [False, True, False]], np.bool8
    )

    with pytest.raises(RuntimeError):
        pybuffer.set_flags(flags_invalid_shape)

    flags = np.array(
        [
            [[True, True, True], [False, True, False]],
            [[False, False, False], [True, False, True]],
        ],
        np.bool8,
    )

    pybuffer.set_flags(flags)
    assert (np.array(pybuffer.get_flags(), copy=False) == flags).all()


def test_data():
    pybuffer = dp3.DPBuffer()

    with pytest.raises(RuntimeError):
        pybuffer.get_data("does not exist")

    data_invalid_shape = np.array(
        [[1.5 + 4j + 2j, 3.0 + 9j, 1.1 + 33j], [27.0 + 2j, 28.0, 29.5]],
        np.csingle,
    )

    with pytest.raises(RuntimeError):
        pybuffer.set_data(data_invalid_shape)

    data = np.array(
        [
            [[1.5 + 4j + 2j, 3.0 + 9j, 1.1 + 33j], [27.0 + 2j, 28.0, 29.5]],
            [[10.3 + 4j, 4.5 + 6j, 7.8 + 3j], [33.1, 0.0 + 1j, 6.6]],
        ],
        np.csingle,
    )

    pybuffer.set_data(data)
    assert (np.array(pybuffer.get_data(), copy=False) == data).all()


def test_uvw():
    pybuffer = dp3.DPBuffer()

    uvw = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]], np.double)

    pybuffer.set_uvw(uvw)

    assert (np.array(pybuffer.get_uvw(), copy=False) == uvw).all()

    uvw_invalid_shape = np.array([[1, 3], [4, 6], [7, 9], [10, 12]], np.double)

    with pytest.raises(RuntimeError):
        pybuffer.set_uvw(uvw_invalid_shape)


def test_process_buffer():
    pybuffer = dp3.DPBuffer()
    step = dp3.make_step(
        "null", dp3.parameterset.ParameterSet(), "", dp3.MsType.regular
    )

    # Using 'step._step' avoids StepWrapper checks.
    step._step.process(pybuffer)

    # After passing a DPBuffer to step.process, it becomes invalid.
    with pytest.raises(RuntimeError):
        pybuffer.get_time()

    with pytest.raises(RuntimeError):
        step._step.process(pybuffer)


def test_extra_data():
    pybuffer = dp3.DPBuffer()

    n_baselines = 1
    n_channels = 3
    n_correlations = 4
    data = np.ones((n_baselines, n_channels, n_correlations), dtype=np.csingle)

    pybuffer.set_data(data)

    modeldata = 5 * np.ones(
        (n_baselines, n_channels, n_correlations), dtype=np.csingle
    )

    # Test that the "model_data" column cannot be filled if it does not exist yet
    with pytest.raises(RuntimeError) as e:
        pybuffer.set_extra_data("model_data", modeldata)
    assert (
        "No data named 'model_data' is found in the current DPBuffer"
        in str(e.value)
    )

    # Add and fill in the "model_data" column
    pybuffer.add_data("model_data")
    pybuffer.set_extra_data("model_data", modeldata)

    assert (pybuffer.get_data() == data).all()
    assert (pybuffer.get_data("model_data") == modeldata).all()

    # Remove the "model_data" column and check that it is correctly removed
    pybuffer.remove_data("model_data")
    with pytest.raises(RuntimeError) as e:
        pybuffer.get_data("model_data")
    assert "Buffer has no data named 'model_data'" in str(e.value)
