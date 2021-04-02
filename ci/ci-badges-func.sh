#!/bin/bash

# Retrieved from https://gitlab.com/ska-telescope/ci-metrics-utilities
# Used to create badges needed for Gitlab by SKA

mkdir -p build/badges
apt-get -y update
apt-get install -y python3-pip python3-setuptools python3-wheel --no-install-recommends
pip3 install anybadge

python3 $(dirname $0)/collect_metrics.py
python3 $(dirname $0)/create_badges.py

if [ $CI_COMMIT_REF_NAME = $CI_DEFAULT_BRANCH ]; then
    # Upload badges and metrics to master folder
    find build/badges -type f -exec curl --user $RAW_USER:$RAW_PASS --upload-file {} $RAW_HOST/repository/raw/gitlab-ci-metrics/$CI_PROJECT_PATH/$CI_COMMIT_REF_NAME/badges/ \;
    find build/reports -type f -exec curl --user $RAW_USER:$RAW_PASS --upload-file {} $RAW_HOST/repository/raw/gitlab-ci-metrics/$CI_PROJECT_PATH/$CI_COMMIT_REF_NAME/reports/ \;

    # Upload badges and metrics to sha folder
    find build/badges -type f -exec curl --user $RAW_USER:$RAW_PASS --upload-file {} $RAW_HOST/repository/raw/gitlab-ci-metrics/$CI_PROJECT_PATH/$CI_COMMIT_SHA/badges/ \;
    find build/reports -type f -exec curl --user $RAW_USER:$RAW_PASS --upload-file {} $RAW_HOST/repository/raw/gitlab-ci-metrics/$CI_PROJECT_PATH/$CI_COMMIT_SHA/reports/ \;
fi
