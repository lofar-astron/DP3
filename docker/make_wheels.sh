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

    #docker build -t dp3-py${py_version} -f docker/py310_wheel.docker .
    #dockerid=$(docker create dp3-py${py_version})
    #docker cp $dockerid:/output/ output-${py_version}
    #docker rm ${dockerid}
    popd
done
