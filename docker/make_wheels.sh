#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script that makes python wheels for several versions.
# The working directory matters since the parent directory should contain the
# source code. The script should be run with
# ./make_wheels.sh <python versions>
# For example, when <python versions> equals 312 314, the script builds wheels
# for python 3.12 and 3.14.
# If <python versions> is empty, it becomes: 314 313 312 311 310 .

set -euo pipefail

if ! DOCKER=$(command -v docker || command -v podman); then
  echo "docker or podman is required for building the wheels"
  exit 1
fi

for py_version in ${@:-314 313 312 311 310}; do
    pushd $(dirname "$0")/..

    ## Build docker image from docker-file. The current wheel is created there.
    [ ${py_version:1} -le 7 ] && py_unicode="m" || py_unicode=
    $DOCKER build \
      -t dp3-py${py_version}-$USER \
      -f docker/py_wheel.docker \
      --progress=plain \
      --build-arg PYMAJOR=${py_version:0:1} \
      --build-arg PYMINOR=${py_version:1} \
      --build-arg PYUNICODE=${py_unicode} .
    ## Create a docker container from that image, and extract the wheel
    containerid=$($DOCKER create dp3-py${py_version}-$USER)
    echo "Docker container ID is: $containerid"
    $DOCKER cp $containerid:/output/ output-${py_version}  # Copies whole dir
    $DOCKER rm ${containerid}
    $DOCKER image rm dp3-py${py_version}-$USER

    popd
done
