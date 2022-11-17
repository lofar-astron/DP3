# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests check that the python bindings for the DPInfo class behave correctly.

This script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tPyDpInfo.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/integration/CMakeLists.txt
"""

import pytest
import sys

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass


def test_antenna_properties():
    info = dp3.DPInfo()

    # Check default values.
    assert info.n_antenna == 0
    assert info.antenna_names == []
    assert info.antenna_positions == []
    assert info.first_antenna_indices == []
    assert info.second_antenna_indices == []

    # Check that properties are read-only.
    with pytest.raises(AttributeError):
        info.n_antenna = 3
    with pytest.raises(AttributeError):
        info.antenna_names = ["very", "nice", "names"]
    with pytest.raises(AttributeError):
        info.antenna_positions = [[0, 0, 0], [0, 42, 0], [42, 0, 0]]
    with pytest.raises(AttributeError):
        info.first_antenna_indices = [0, 0, 0]
    with pytest.raises(AttributeError):
        info.second_antenna_indices = [0, 1, 2]

    # Check that set_antennas() yields new property values.
    names = ["name_1", "name_2"]
    diameters = [42, 43]
    positions = [[1, 2, 3], [4, 5, 6]]
    first_indices = [0]
    second_indices = [1]
    info.set_antennas(
        names, diameters, positions, first_indices, second_indices
    )
    assert info.n_antenna == 2
    assert info.antenna_names == names
    assert info.antenna_positions == positions
    assert info.first_antenna_indices == first_indices
    assert info.second_antenna_indices == second_indices


def test_channel_properties():
    info = dp3.DPInfo()

    # Check default values.
    assert info.n_channels == 0
    assert info.channel_frequencies == []
    assert info.bda_channel_frequencies == [[]]
    assert info.channel_widths == []
    assert info.bda_channel_widths == [[]]

    # Check that properties are read-only.
    with pytest.raises(AttributeError):
        info.n_channels = 3
    with pytest.raises(AttributeError):
        info.channel_frequencies = [10.0e6, 11.0e6, 12.0e6]
    with pytest.raises(AttributeError):
        info.bda_channel_frequencies = [[13.0e6, 14.0e6, 15.0e6]]
    with pytest.raises(AttributeError):
        info.channel_frequencies = [2.0e6, 2.0e6, 2.0e6]
    with pytest.raises(AttributeError):
        info.bda_channel_frequencies = [[3.0e6, 3.0e6, 3.0e6]]

    # Check that set_channels() yields new property values.
    frequencies = [42.0e6, 43.0e6]
    widths = [1.0e6, 1.0e6]
    info.set_channels(frequencies, widths)
    assert info.n_channels == 2
    assert info.channel_frequencies == frequencies
    assert info.bda_channel_frequencies == [frequencies]
    assert info.channel_widths == widths
    assert info.bda_channel_widths == [widths]
