# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests are checking that the python bindings for the ParameterSet class
behave correctly.

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/unit directory,
  using `pytest source/tPyParameterSet.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/unit/CMakeLists.txt
"""

import sys

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3.parameterset
except ImportError:
    pass


def test_parameter_set_magic():
    """
    Test the magic (double underscore) member functions  of ParamterSet and ParameterValue

    Tests the member functions __contains__ and __getitem__ of ParameterSet, and
    __bool__, __float__, __int__, __iter__ and __str__ of ParameterValue
    """
    ps = dp3.parameterset.ParameterSet()

    ps.add("prefix.somekey", "42")
    ps.add("prefix.vecint", "[1,2,3]")
    ps.add("prefix.vecfloat", "[2.71828, 3.14159]")
    ps.add("prefix.vecbool", "[T, F, True, False, 1, 0]")
    ps.add("prefix.vecstr", "[CS001,CS002,CS003]")
    ps.add("prefix.vecvecstr", "[[CS001,CS002],[CS003]]")

    assert "prefix.not_existent_key" not in ps
    assert "prefix.somekey" in ps

    assert [int(v) for v in ps["prefix.vecint"]] == [1, 2, 3]
    assert [str(v) for v in ps["prefix.vecstr"]] == ["CS001", "CS002", "CS003"]
    assert [float(v) for v in ps["prefix.vecfloat"]] == [2.71828, 3.14159]
    assert [bool(v) for v in ps["prefix.vecbool"]] == [
        True,
        False,
        True,
        False,
        True,
        False,
    ]
    assert [[str(vv) for vv in v] for v in ps["prefix.vecvecstr"]] == [
        ["CS001", "CS002"],
        ["CS003"],
    ]
