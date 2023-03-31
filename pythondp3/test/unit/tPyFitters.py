# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests are checking that the python bindings for the Constraint classes
behave correctly.

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/unit directory,
  using `pytest source/tPyFitters.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/unit/CMakeLists.txt
"""

import numpy as np
import sys

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3.fitters
except ImportError:
    pass

TEC_CONSTANT = -8.44797245e9
N_FREQUENCIES = 42
N_ANTENNAS = 5


def test_fit_single_antenna():
    frequencies = np.linspace(100e6, 110e6, N_FREQUENCIES)
    tec_value = 0.42
    clean_gains = np.exp((0 + 1j) * tec_value * TEC_CONSTANT / frequencies)

    np.random.seed(42)
    noise_amplitude = 0.1
    noise = np.exp(
        (0 + 1j)
        * np.random.normal(scale=noise_amplitude, size=(N_FREQUENCIES))
    )
    noisy_gains = clean_gains * noise
    assert np.isclose(np.abs(noisy_gains), 1).all()

    tec_result, phase_result, error_result = dp3.fitters.fit_tec(
        noisy_gains, frequencies
    )

    assert tec_result.name == "tec"
    assert phase_result.name == "phase"
    assert error_result.name == "error"

    assert tec_result.axes == "ant"
    assert phase_result.axes == "ant"
    assert error_result.axes == "ant"

    assert tec_result.shape == [1]
    assert phase_result.shape == [1]
    assert error_result.shape == [1]

    assert str(tec_result) == "[Result name:tec axes:ant shape:1]"
    assert str(phase_result) == "[Result name:phase axes:ant shape:1]"
    assert str(error_result) == "[Result name:error axes:ant shape:1]"

    assert (
        repr(tec_result)
        == "<dp3.fitters.Result name=tec axes=ant values=[0.41872,] weights=[42,]>"
    )
    assert (
        repr(phase_result)
        == "<dp3.fitters.Result name=phase axes=ant values=[-0.121689,] weights=[42,]>"
    )
    assert (
        repr(error_result)
        == "<dp3.fitters.Result name=error axes=ant values=[0.074679,] weights=[42,]>"
    )

    assert np.isclose(tec_result.values[0], tec_value, rtol=1e-2)
    assert np.isclose(phase_result.values[0], 0.0, atol=1.0)  # Tolerance: 1 Hz
    assert np.isclose(error_result.values[0], 0.0, atol=0.1)

    assert tec_result.weights == [N_FREQUENCIES]
    assert phase_result.weights == [N_FREQUENCIES]
    assert error_result.weights == [N_FREQUENCIES]

    assert np.isclose(np.abs(noisy_gains), 1).all()
    diff_angles = np.angle(clean_gains) - np.angle(noisy_gains)
    atol_angle = 0.03

    for diff_angle in np.ndarray.flatten(diff_angles):
        # The phase difference should be close to 0, 2*PI or -2*PI.
        assert (
            np.isclose(diff_angle, 0, atol=atol_angle)
            or np.isclose(diff_angle, 2 * np.pi, atol=atol_angle)
            or np.isclose(diff_angle, -2 * np.pi, atol=atol_angle)
        )


def test_fit_tec_only():
    frequencies = np.linspace(120e6, 150e6, N_FREQUENCIES)
    tec_values = np.linspace(0, N_ANTENNAS - 1, N_ANTENNAS)
    clean_phases = tec_values * TEC_CONSTANT / frequencies[:, np.newaxis]
    clean_gains = np.exp((0 + 1j) * clean_phases)

    np.random.seed(42)
    noise_amplitude = 0.1
    noise = np.exp(
        (0 + 1j)
        * np.random.normal(
            scale=noise_amplitude, size=(N_FREQUENCIES, N_ANTENNAS)
        )
    )
    noisy_gains = clean_gains * noise
    assert np.isclose(np.abs(noisy_gains), 1).all()

    tec_result, error_result = dp3.fitters.fit_tec(
        noisy_gains, frequencies, mode="tec_only"
    )

    assert tec_result.name == "tec"
    assert np.isclose(tec_result.values, tec_values, rtol=1e-3).all()
    assert tec_result.weights == [N_FREQUENCIES] * N_ANTENNAS
    assert tec_result.axes == "ant"
    assert tec_result.shape == [N_ANTENNAS]
    assert str(tec_result) == f"[Result name:tec axes:ant shape:{N_ANTENNAS}]"

    assert error_result.name == "error"
    assert np.isclose(error_result.values, 0.0, atol=0.15).all()
    assert error_result.weights == [N_FREQUENCIES] * N_ANTENNAS
    assert error_result.axes == "ant"
    assert error_result.shape == [N_ANTENNAS]
    assert (
        str(error_result) == f"[Result name:error axes:ant shape:{N_ANTENNAS}]"
    )

    assert np.isclose(np.abs(noisy_gains), 1).all()
    diff_angles = np.angle(clean_gains) - np.angle(noisy_gains)
    atol_angle = 0.05

    for diff_angle in np.ndarray.flatten(diff_angles):
        # The phase difference should be close to 0, 2*PI or -2*PI.
        assert (
            np.isclose(diff_angle, 0, atol=atol_angle)
            or np.isclose(diff_angle, 2 * np.pi, atol=atol_angle)
            or np.isclose(diff_angle, -2 * np.pi, atol=atol_angle)
        )
