# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
from subprocess import check_call, check_output
from testconfig import TAQLEXE


def assert_taql(command, expected_rows=0):
    result = check_output([TAQLEXE, "-noph", command]).decode().strip()
    assert result == f"select result of {expected_rows} rows"


def untar_ms(source):
    if not os.path.isfile(source):
        raise IOError(f"Not able to find {source} containing the reference solutions.")
    check_call(["tar", "xf", source])
