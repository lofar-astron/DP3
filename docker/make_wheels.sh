#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script that makes python wheels for several versions.
# The working directory matters since the parent directory should contain the
# source code. The script should be run with
# ./make_wheels.sh <python versions>
# For example, when <python versions> equals 38 310, the script builds wheels
# for python 3.8 and 3.10.
# If <python versions> is empty, it becomes: 310 39 38 37 36 .

set -euo pipefail
for py_version in ${@:-313 312 311 310 39 38}; do
    pushd ..

    ## Build docker image from docker-file. The current wheel is created there.
    [ ${py_version:1} -le 7 ] && py_unicode="m" || py_unicode=
    docker build -t dp3-py${py_version}-$USER -f docker/py_wheel.docker \
      --build-arg PYMAJOR=${py_version:0:1} \
      --build-arg PYMINOR=${py_version:1} \
      --build-arg PYUNICODE=${py_unicode} .
    ## Create a docker container from that image, and extract the wheel
    containerid=$(docker create dp3-py${py_version}-$USER)
    echo "Docker container ID is: $containerid"
    docker cp $containerid:/output/ output-${py_version}  # Copies whole dir
    docker rm ${containerid}
    docker image rm dp3-py${py_version}-$USER

    popd
done
