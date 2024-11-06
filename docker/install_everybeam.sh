#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install EveryBeam from source

set -euo pipefail

EVERYBEAM_VERSION=$1

# On CentOS7, the 'blas' package does not provide cblas functions.
# -> Install OpenBLAS, which does provide them. CMake will prefer OpenBLAS.
yum -y install openblas-devel wget

pushd /tmp

echo "Cloning EveryBeam ${EVERYBEAM_VERSION}"

git clone --branch v${EVERYBEAM_VERSION} --depth 1 \
  https://git.astron.nl/RD/EveryBeam.git --recursive --shallow-submodules

echo "Configuring, building & installing EveryBeam ${EVERYBEAM_VERSION}"
pushd EveryBeam

# Ensure EveryBeam does not use python.
sed -i '/find_package(PythonInterp .*)/d' CMakeLists.txt
sed -i '/pybind11/d' CMakeLists.txt

mkdir build
cd build
# On CentOS 7, the default OpenBLAS version is the serial version, which is
# incompatible with DP3. -> Use the version with *p*threads support.
cmake -DBLAS_LIBRARIES=/usr/lib64/libopenblasp.so ..
make -j${THREADS} install
popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary EveryBeam files"
rm -r EveryBeam

popd

echo "Cleaning up development libraries"
yum -y erase openblas-devel wget
