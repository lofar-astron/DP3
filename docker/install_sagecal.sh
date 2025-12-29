#!/bin/bash

# Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Installs sagecal (libdirac) from source.

set -euo pipefail

# The version argument should not include the 'v' prefix, e.g., "0.8.5".
SAGECAL_VERSION=$1

pushd /tmp

echo "Downloading & unpacking SAGECal ${SAGECAL_VERSION}"
curl -fsSL -o sagecal-${SAGECAL_VERSION}.tgz https://github.com/nlesc-dirac/sagecal/archive/refs/tags/v${SAGECAL_VERSION}.tar.gz
tar xfz sagecal-${SAGECAL_VERSION}.tgz

pushd sagecal-${SAGECAL_VERSION}

echo "Configuring, building & installing SAGECal ${SAGECAL_VERSION}"
mkdir build
cd build
cmake -DLIB_ONLY=1 ..
make -j${THREADS} install

popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary SAGECal files"
rm -r sagecal-${SAGECAL_VERSION}.tgz sagecal-${SAGECAL_VERSION}

popd
