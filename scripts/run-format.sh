#!/bin/bash

# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#Script configuration for this repo. Adjust it when copying to a different repo.

cd $(dirname "$0")/..
DIRECTORY_LIST=`ls -1d */|grep -v 'build\|ci\|CMake\|CPack\|docker\|docs\|external\|pythondp3\|resources\|scripts'`

# Disable globbing. This is needed when defining patterns that have wildcards.
set -e -f

#The directory that contains the source files.
SOURCE_DIR=$(dirname "$0")/..

#Directories that must be excluded from formatting. These paths are
#relative to SOURCE_DIR.
EXCLUDE_DIRS=(external build CMake)

#The patterns of the C++ source files, which clang-format should format.
CXX_SOURCES=(*.cc *.h *.cu *.cuh)

# clang-format version 12 and 14 produce slightly different results for
# pythondp3/parameterset.cc, so hard-code this to 14.
CLANG_FORMAT_BINARY=clang-format-14

#End script configuration.

if grep 'cout\|cerr' -r ${DIRECTORY_LIST} ; then
    echo -e "\e[1m\e[31mAt least one file makes use of 'cout' or 'cerr'\e[0m"
    echo -e "\e[1m\e[31mUse aocommon::Logger instead of cout/cerr.\e[0m"
    exit 1
fi

#The common formatting script has further documentation.
source $(dirname "$0")/../external/aocommon/scripts/format.sh
