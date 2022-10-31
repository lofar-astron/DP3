# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

sys.path.append(f"{tcf.BINDIR}/pythondp3")

"""
Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tImportDP3.py` (extended with pytest options of your choice)
- using ctest, see DP3/pythondp3/test/integration/CMakeLists.txt
"""


def test_import_dp3():
    """
    Test import of dp3 module
    One reason this might fail is when the linkage of this module, a shared library,
    is incomplete. In that case, the run-time linker will report unresolved symbols.
    """
    import dp3
