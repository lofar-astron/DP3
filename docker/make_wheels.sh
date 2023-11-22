#!/usr/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to make python wheels for several versions
# The working directory matters, script should be run with ./make_wheels.sh
# The parent directory should contain the source code with the version to build.

set -euo pipefail
for py_version in 310 39 38 37 36; do
    pushd ..
    sed -i "s=\(master_wheel\)[0-9]*=\1${py_version}=" docker/py310_wheel.docker
    grep wheel3 docker/py310_wheel.docker

    ## Build docker image from docker-file. The current wheel is created there
    #docker build -t dp3-py${py_version}-$USER -f docker/py310_wheel.docker .
    ## Create a docker container from that image, and extract the wheel
    #containerid=$(docker create dp3-py${py_version}-$USER)
    #echo "Docker container ID is: $containerid"
    #docker cp $containerid:/output/ output-${py_version}  # Copies whole dir
    #docker rm ${containerid}
    #docker image rm dp3-py${py_version}-$USER
    popd
done
