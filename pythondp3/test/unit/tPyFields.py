# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests are checking that the python bindings for the Fields class behave
correctly. The tests check the same functionalities checked in
DP3/common/test/unit/tFields.cc for the C++ version of the Fields class.

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/unit directory,
  using `pytest source/tPyFields.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/unit/CMakeLists.txt
"""

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass


def import_fields():
    """Validates the module can be imported.

    Importing a module may fail while running pytest --collect-only,
    but it should work when running the tests."""
    import dp3


def check_fields(fields, set_fields):
    """Check whether the fields set in `fields` (a `Fields` object) are those
    in `set_fields` (a list of strings)."""

    assert fields.data == ("data" in set_fields)
    assert fields.flags == ("flags" in set_fields)
    assert fields.weights == ("weights" in set_fields)
    assert fields.uvw == ("uvw" in set_fields)
    assert fields.fullresflags == ("fullResFlags" in set_fields)


def test_constructor():
    # Test constructor with empty field
    empty_fields = dp3.Fields()
    check_fields(empty_fields, [])

    # Test constructor with single field
    data_field = dp3.Fields(dp3.Fields.Single.DATA)
    check_fields(data_field, ["data"])

    flag_field = dp3.Fields(dp3.Fields.Single.FLAGS)
    check_fields(flag_field, ["flags"])

    weights_field = dp3.Fields(dp3.Fields.Single.WEIGHTS)
    check_fields(weights_field, ["weights"])

    uvw_field = dp3.Fields(dp3.Fields.Single.UVW)
    check_fields(uvw_field, ["uvw"])

    fullresflags_field = dp3.Fields(dp3.Fields.Single.FULLRESFLAGS)
    check_fields(fullresflags_field, ["fullResFlags"])

    # Test constructor with combined fields
    check_fields(flag_field | data_field, ["data", "flags"])


def test_or_operator():
    test_fields = dp3.Fields()
    test_fields |= dp3.Fields.FLAGS
    check_fields(test_fields, ["flags"])

    test_fields |= dp3.Fields.FULLRESFLAGS
    check_fields(test_fields, ["flags", "fullResFlags"])

    data_field = dp3.Fields.DATA
    uvw_field = dp3.Fields.UVW

    check_fields(data_field | uvw_field, ["data", "uvw"])


def test_equality_operator():
    left_empty = dp3.Fields()
    right_empty = dp3.Fields()
    left_three = dp3.Fields.DATA | dp3.Fields.FLAGS | dp3.Fields.UVW
    right_three = dp3.Fields.DATA | dp3.Fields.FLAGS | dp3.Fields.UVW

    # Check operator ==
    assert left_empty == left_empty
    assert left_empty == right_empty
    assert right_three == right_three
    assert left_three == right_three
    assert ~(left_empty == right_three)
    assert ~(right_three == left_empty)

    # Check operator!=
    assert ~(left_empty != left_empty)
    assert ~(left_empty != right_empty)
    assert ~(right_three != right_three)
    assert ~(left_three != right_three)
    assert left_empty != right_three
    assert right_three != left_empty


def test_string_operator():
    assert str(dp3.Fields()) == "[]"
    assert str(dp3.Fields.DATA) == "[data]"
    assert str(dp3.Fields.FLAGS) == "[flags]"
    assert (
        str(
            dp3.Fields.DATA
            | dp3.Fields.FLAGS
            | dp3.Fields.WEIGHTS
            | dp3.Fields.FULLRESFLAGS
            | dp3.Fields.UVW
        )
        == "[data, flags, weights, fullresflags, uvw]"
    )


def test_update_requirements():
    check_fields(
        dp3.Fields(dp3.Fields.DATA).update_requirements(
            dp3.Fields(), dp3.Fields.DATA
        ),
        [],
    )
    check_fields(
        dp3.Fields(dp3.Fields.DATA).update_requirements(
            dp3.Fields.DATA,
            dp3.Fields.DATA,
        ),
        ["data"],
    )

    check_fields(
        dp3.Fields(dp3.Fields.DATA).update_requirements(
            dp3.Fields.DATA, dp3.Fields()
        ),
        ["data"],
    )

    check_fields(
        dp3.Fields(dp3.Fields.DATA).update_requirements(
            dp3.Fields(), dp3.Fields()
        ),
        ["data"],
    )

    check_fields(
        dp3.Fields(dp3.Fields.DATA).update_requirements(
            dp3.Fields(), dp3.Fields.UVW
        ),
        ["data"],
    )
    check_fields(
        dp3.Fields(dp3.Fields.DATA).update_requirements(
            dp3.Fields.UVW,
            dp3.Fields.UVW,
        ),
        ["data", "uvw"],
    )
    check_fields(
        dp3.Fields(dp3.Fields.DATA).update_requirements(
            dp3.Fields.UVW, dp3.Fields()
        ),
        ["data", "uvw"],
    )
