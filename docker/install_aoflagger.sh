#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install AOFlagger from source.

set -euo pipefail

AOFLAGGER_VERSION=$1

echo "Installing development libraries with yum"
yum -y install libpng-devel zlib-devel

pushd /tmp

echo "Downloading & unpacking AOFlagger ${AOFLAGGER_VERSION}"
if [ "${AOFLAGGER_VERSION}" != "3.4.0" ]; then
  echo $0: Please update the link for downloading AOFlagger ${AOFLAGGER_VERSION}.
  exit 1
fi
curl -fsSL -o aoflagger-v${AOFLAGGER_VERSION}.bz2 "https://gitlab.com/aroffringa/aoflagger/-/package_files/96704214/download"
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
yum -y erase libpng-devel zlib-devel
