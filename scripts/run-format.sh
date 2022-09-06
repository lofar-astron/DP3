#!/bin/bash

# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#Script configuration for this repo. Adjust it when copying to a different repo.

#The directory that contains the source files.
SOURCE_DIR=$(dirname "$0")/..

#Directories that must be excluded from formatting. These paths are
#relative to SOURCE_DIR.
EXCLUDE_DIRS=(external build CMake)

# clang-format version 12 and 14 produce slightly different results for
# pythondp3/parameterset.cc, so hard-code this to 14.
CLANG_FORMAT_BINARY=clang-format-14

#End script configuration.

#The common formatting script has further documentation.
source $(dirname "$0")/../external/aocommon/scripts/format.sh
