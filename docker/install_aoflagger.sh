#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install AOFlagger from source.

set -euo pipefail

AOFLAGGER_VERSION=$1

echo "Installing development libraries with yum"
# In the Almalinux 8.10 base image, pybind11-devel is tagged with a Python version.
# Since pybind11 is a header-only library, it works with any Python version.
yum -y install libpng-devel python3.11-pybind11-devel zlib-devel

pushd /tmp

echo "Downloading & unpacking AOFlagger ${AOFLAGGER_VERSION}"
# Download a custom source tarball, which includes external submodules like
# aocommon and schaapcommon.
if [ "${AOFLAGGER_VERSION}" != "3.5.1" ]; then
  echo $0: Please update the link for downloading AOFlagger ${AOFLAGGER_VERSION}.
  exit 1
fi
curl -fsSL -o aoflagger-v${AOFLAGGER_VERSION}.bz2 "https://gitlab.com/aroffringa/aoflagger/-/package_files/264635245/download"
tar -xjf aoflagger-v${AOFLAGGER_VERSION}.bz2

pushd aoflagger-v${AOFLAGGER_VERSION}

# Wheels should not actually link to libpython (see py_wheel.docker).
sed -i '/^ *\${Python_LIBRARIES/d' CMakeLists.txt

echo "Configuring, building & installing AOFlagger ${AOFLAGGER_VERSION}"
mkdir build
cd build
#-DFFTW3_LIB=${FFTW_DIR}/lib/libfftw3.so \
cmake \
  -DENABLE_GUI=False \
  -DPython_EXECUTABLE=/opt/python/${TARGET}/bin/python \
  ..
make -j${THREADS} aoflagger-lib
cmake --install . --component core_library
popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary AOFlagger files"
rm -r aoflagger-v${AOFLAGGER_VERSION}
rm aoflagger-v${AOFLAGGER_VERSION}.bz2

popd

echo "Cleaning up development libraries"
yum -y erase libpng-devel python3.11-pybind11-devel zlib-devel
