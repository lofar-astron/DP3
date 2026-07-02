# Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Tests for the EveryBeam CreateTelescope binding exposed through DP3.

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/unit directory,
  using `pytest source/tPyCreateTelescope.py` (extended with pytest options of
  your choice)
- using ctest, see pythondp3/test/unit/CMakeLists.txt
"""

import os
import subprocess
import sys

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass


def test_everybeam_namespace_is_exported():
    import dp3.everybeam as everybeam

    assert everybeam is dp3.everybeam
    assert dp3.everybeam is dp3.pydp3.everybeam
    assert callable(everybeam.create_telescope)
    assert everybeam.create_telescope is dp3.pydp3.everybeam.create_telescope
    assert everybeam.Options is dp3.pydp3.everybeam.Options
    assert everybeam.StationNode is dp3.pydp3.everybeam.StationNode
    assert everybeam.Telescope is dp3.pydp3.everybeam.Telescope
    assert everybeam.TelescopeType is dp3.pydp3.everybeam.TelescopeType
    assert not hasattr(dp3, "create_telescope")
    assert not hasattr(dp3, "Options")

    docstring = everybeam.create_telescope.__doc__
    assert "telescope_type" in docstring
    assert "options" in docstring
    assert "station_tree" in docstring
    assert "preapplied_beam_mode" in docstring


def test_telescope_type_values_are_exported():
    everybeam = dp3.everybeam

    assert everybeam.TelescopeType.AARTFAAC.name == "AARTFAAC"
    assert everybeam.TelescopeType.LOFAR.name == "LOFAR"
    assert everybeam.TelescopeType.OSKAR.name == "OSKAR"
    assert everybeam.TelescopeType.SKA_MID.name == "SKA_MID"


def test_create_telescope_uses_dp3_everybeam_types():
    everybeam = dp3.everybeam

    options = everybeam.Options()
    options.element_response_model = (
        everybeam.ElementResponseModel.oskar_dipole
    )
    station_tree = everybeam.StationNode()

    telescope = everybeam.create_telescope(
        everybeam.TelescopeType.OSKAR,
        options,
        station_tree,
        delay_directions=[[0.0, 0.0]],
    )

    assert isinstance(telescope, everybeam.Telescope)

    info = dp3.DPInfo()
    assert not info.has_telescope

    info.set_telescope(telescope)
    assert info.has_telescope
    assert isinstance(info.get_telescope(), everybeam.Telescope)

    info.set_telescope(None)
    assert not info.has_telescope


def test_create_telescope_accepts_station_tree():
    everybeam = dp3.everybeam

    options = everybeam.Options()
    options.element_response_model = (
        everybeam.ElementResponseModel.oskar_dipole
    )

    axes = everybeam.StationCoordinateSystemAxes()
    coordinate_system = everybeam.StationCoordinateSystem(
        [0.0, 0.0, 0.0], axes
    )
    child = everybeam.StationNode()
    child.add_child_element([0.0, 0.0, 0.0])

    station_tree = everybeam.StationNode(coordinate_system, "station")
    station_tree.add_child_node(child, [1.0, 2.0, 3.0])

    telescope = everybeam.create_telescope(
        everybeam.TelescopeType.OSKAR,
        options,
        station_tree,
        delay_directions=[[0.0, 0.0]],
    )

    assert isinstance(telescope, everybeam.Telescope)


def test_create_telescope_requires_telescope_type():
    everybeam = dp3.everybeam

    options = everybeam.Options()
    station_tree = everybeam.StationNode()

    with pytest.raises(TypeError):
        everybeam.create_telescope(
            options, station_tree, delay_directions=[[0.0, 0.0]]
        )


def test_create_telescope_isolated_from_everybeam_python_module():
    pytest.importorskip("everybeam")

    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [tcf.PYTHONDIR, environment.get("PYTHONPATH", "")]
    )
    subprocess.run(
        [
            sys.executable,
            "-c",
            """
import everybeam
import dp3

options = dp3.everybeam.Options()
options.element_response_model = dp3.everybeam.ElementResponseModel.oskar_dipole
telescope = dp3.everybeam.create_telescope(
    dp3.everybeam.TelescopeType.OSKAR,
    options,
    dp3.everybeam.StationNode(),
    delay_directions=[[0.0, 0.0]],
)
info = dp3.DPInfo()
info.set_telescope(telescope)
assert type(telescope).__module__ == "dp3.pydp3.everybeam"
assert info.has_telescope
""",
        ],
        check=True,
        env=environment,
    )
