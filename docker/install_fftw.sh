#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install FFTW from source (both single and double precision).

set -euo pipefail

FFTW_VERSION=$1

pushd /tmp

echo "Downloading & unpacking FFTW ${FFTW_VERSION}"
curl https://fftw.org/pub/fftw/fftw-${FFTW_VERSION}.tar.gz --output fftw-${FFTW_VERSION}.tar.gz
tar xf fftw-${FFTW_VERSION}.tar.gz

pushd fftw-${FFTW_VERSION}
echo "Configuring, building & installing FFTW ${FFTW_VERSION}"

#Single precision
./configure --enable-threads --enable-shared --enable-float
make -j${THREADS}
make install

#Double precision (default)
./configure --enable-threads --enable-shared
make -j${THREADS}
make install
popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary files"
rm -r fftw-${FFTW_VERSION}
rm fftw-${FFTW_VERSION}.tar.gz

popd