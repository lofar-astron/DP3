#!/bin/bash

# Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install GSL from source. The manylinux wheel base image ships an
# older GSL that lacks gsl_multilarge_nlinear.h.

set -euo pipefail

GSL_VERSION=$1

pushd /tmp

echo "Downloading & unpacking GSL ${GSL_VERSION}"
curl -fsSLO "https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz"
tar -xzf gsl-${GSL_VERSION}.tar.gz

pushd gsl-${GSL_VERSION}
echo "Configuring, building & installing GSL ${GSL_VERSION}"
./configure --prefix=/usr/local
make -j${THREADS}
make install
ldconfig
popd

echo "Cleaning up unnecessary GSL files"
rm -r gsl-${GSL_VERSION}
rm gsl-${GSL_VERSION}.tar.gz

popd
