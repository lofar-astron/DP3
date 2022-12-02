#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install extra boost libraries from source

set -euo pipefail

# Boost is already in the container as leftover from casacore install
# AOFlagger v3.2.0 always requires the unit test framework. More recent versions
# no longer have this requirement.
pushd /build/boost_${BOOST_}
./bootstrap.sh --with-libraries=date_time,filesystem,math,program_options,system,test
./b2 -j${THREADS} cxxflags="-fPIC" link=static,shared install
popd
