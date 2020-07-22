#!/bin/bash
#
# run-clang-format.sh: Formats source code in this repo in accordance with .clang-format file.
# This file is part of the DP3 software package.
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# License: GNU General Public License version 3 or any later version

set -e

# print in bold-face
echo -e "\e[1mRunning clang-format...\e[0m"

find $(dirname "$0")/.. -type f -name *.cc -o -name *.h \
  -exec clang-format -i -style=file \{\} +

if git diff --exit-code --quiet; then
    # print in bold-face green
    echo -e "\e[1m\e[32mGreat job, git shows no changed files!\e[0m"
    exit 0;
else
    # Print in bold-face red
    echo -e "\e[1m\e[31mGit shows at least one changed file now!\e[0m"
    exit 1;
fi
