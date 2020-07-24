#!/bin/bash
#
# run-clang-format.sh: Formats source code in this repo in accordance with .clang-format file.
# This file is part of the DP3 software package.
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# License: GNU General Public License version 3 or any later version
#
# To hook this script to pre-commit include the line
# "./scripts/run-clang-format.sh" to .git/hooks/pre-commit
# and make sure pre-commit is an executable shell script.

#Script configuration for this repo. Adjust it when copying to a different repo.

#The directory that contains the source files, which clang-format should format.
SOURCE_DIR=$(dirname "$0")/..

#Directories that must be excluded from formatting. These paths are
#relative to SOURCE_DIR.
EXCLUDE_DIRS=(external)

#The extensions of the source files, which clang-format should format.
SOURCE_EXT=(*.cc *.h)

#End script configuration.

set -e

# print in bold-face
echo -e "\e[1mRunning clang-format...\e[0m"

# Convert SOURCE_EXT into "-name ext1 -o -name ext2 -o name ext3 ..."
FIND_NAMES="-name ${SOURCE_EXT[0]}"
for i in `seq 1 $((${#SOURCE_EXT[*]} - 1))`; do
  FIND_NAMES+=" -o -name ${SOURCE_EXT[$i]}"
done

# Convert EXCLUDE_DIRS into "-path ./dir1 -prune -o -path ./dir2 -prune -o ..."
FIND_EXCLUDES=
for e in ${EXCLUDE_DIRS[*]}; do
  FIND_EXCLUDES+="-path ./$e -prune -o "
done

cd $SOURCE_DIR
find . $FIND_EXCLUDES -type f \( $FIND_NAMES \) \
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
