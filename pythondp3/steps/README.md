[//]: # "Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy"
[//]: # "SPDX-License-Identifier: GPL-3.0-or-later"

This directory is there to collect all Steps implemented in Python.

There is no `test` subdirectory here. Tests should be placed in
either `pythondp3/test/unit` or `pythondp3/test/integration`.

The `unit` subdirectory contains both C++ tests, which use the Boost Test framework,
and Python test, which use the pytest test framework. The tests in the `integration`
subdirectory are written in Python and only use pytest as test framework.

It is possible to test Python steps from C++ by
instantiating them with the DP3::MakeSingleStep() function
as is done for example in `pythondp3/test/unit/tPyStep.cc`
But it probably more convenient to test them in Python 
by creating a pytest test in `pythondp3/test/integration`

When adding a Step in this directory, please update `pythondp3/CMakeLists.txt`
to make sure it is installed both in the pseudo install directory in the build
environment and in the true install directory.
