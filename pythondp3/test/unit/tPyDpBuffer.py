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
