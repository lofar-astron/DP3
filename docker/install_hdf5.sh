#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install HDF5 from source. CentOS 7 has a too old version.

set -euo pipefail

HDF5_VERSION=$1

echo "Installing zlib with yum"
yum -y install zlib-devel

pushd /tmp

echo "Downloading & unpacking HDF5 ${HDF5_VERSION}"
#                                   Remove trailing .*, to get e.g. '1.12' ↓
curl -fsSLO "https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_${HDF5_VERSION//./_}/source/hdf5-${HDF5_VERSION}.tar.gz"
tar -xzvf hdf5-${HDF5_VERSION}.tar.gz

pushd hdf5-${HDF5_VERSION}
echo "Configuring, building & installing HDF5 ${HDF5_VERSION}"
mkdir build
cd build
# Overriding CMAKE_INSTALL_PREFIX is necessary since HDF5 installs into
# /usr/local/HDF_Group/HDF5/${HDF5_VERSION} by default.
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=Off -DHDF5_BUILD_CPP_LIB=On ..
make -j${THREADS}
make install
popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary files"
rm -r hdf5-${HDF5_VERSION}
rm hdf5-${HDF5_VERSION}.tar.gz

popd

echo "Cleaning up development libraries"
yum -y erase zlib-devel
