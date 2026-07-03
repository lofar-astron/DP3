#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install EveryBeam from source

set -euo pipefail

EVERYBEAM_VERSION=$1

# On CentOS7, the 'blas' package does not provide cblas functions.
# -> Install OpenBLAS, which does provide them. CMake will prefer OpenBLAS.
# Install openblas-threads explicitly, otherwise erasing openblas-devel below
# will also erase openblas-threads.
yum -y install openblas-devel openblas-threads wget

pushd /tmp

echo "Cloning EveryBeam ${EVERYBEAM_VERSION}"

git clone --branch v${EVERYBEAM_VERSION} --depth 1 \
  https://git.astron.nl/RD/EveryBeam.git --recursive --shallow-submodules

echo "Configuring, building & installing EveryBeam ${EVERYBEAM_VERSION}"
pushd EveryBeam

mkdir build
cd build
# Ensure EveryBeam does not use python.
CMAKE_FLAGS="-DBUILD_WITH_PYTHON=OFF"
# Keep the vendored EveryBeam symbols in DP3 wheels separate from everybeam wheels.
CMAKE_FLAGS+=" -DCMAKE_CXX_FLAGS=-Deverybeam=dp3::everybeam"

# On AlmaLinux 8, the default OpenBLAS version is the serial version, which is
# incompatible with DP3. -> Use the version with *p*threads support.
CMAKE_FLAGS+=" -DBLAS_LIBRARIES=/usr/lib64/libopenblasp.so"
cmake ${CMAKE_FLAGS} ..
make -j${THREADS} install
popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary EveryBeam files"
rm -r EveryBeam

popd

echo "Cleaning up development libraries"
yum -y erase openblas-devel wget
