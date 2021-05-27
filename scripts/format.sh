# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# format.sh: Formats source code in a repository in accordance with
# .clang-format and .cmake-format.py files.
#
# This script uses the following variables:
# - SOURCE_DIR: The directory that contains the source files.
# - EXCLUDE_DIRS: (Optional) directories that must be excluded from formatting.
#                 These paths are relative to SOURCE_DIR.
# - CXX_SOURCES: Patterns of the C++ files, which clang-format should format.
# - CMAKE_SOURCES: Patterns of the CMake files, which cmake-format should format.
#
# A repository that uses format.sh should define its own run-format.sh script
# that defines these variables and then sources this script.
# If you want to automatically check formatting in each commit, include the line
# "./scripts/run-format.sh" to .git/hooks/pre-commit
# and make sure pre-commit is an executable shell script.

# Disable globbing
set -e -f

# Check arguments
if [ -z "$SOURCE_DIR" -o -z "$CXX_SOURCES" -o -z "$CMAKE_SOURCES" ]; then
  echo "Please define SOURCE_DIR, CXX_SOURCES and CMAKE_SOURCES when using $BASH_SOURCE"
  exit 1
fi

# Detect run environment.
if [ -n "$CI" ]; then
  DRYRUN=" (dry run on CI)"
elif [ -n "$GIT_AUTHOR_DATE" ]; then
  DRYRUN=" (dry run in git hook)"
fi

# print in bold-face
echo -e "\e[1mRunning formatters$DRYRUN...\e[0m"

# Convert SOURCES into "-name ext1 -o -name ext2 -o name ext3 ..."
CXX_FIND_NAMES="-name ${CXX_SOURCES[0]}"
for i in `seq 1 $((${#CXX_SOURCES[*]} - 1))`; do
  CXX_FIND_NAMES+=" -o -name ${CXX_SOURCES[$i]}"
done

CMAKE_FIND_NAMES="-name ${CMAKE_SOURCES[0]}"
for i in `seq 1 $((${#CMAKE_SOURCES[*]} - 1))`; do
  CMAKE_FIND_NAMES+=" -o -name ${CMAKE_SOURCES[$i]}"
done

# Convert EXCLUDE_DIRS into "-path ./dir1 -prune -o -path ./dir2 -prune -o ..."
FIND_EXCLUDES=
for e in ${EXCLUDE_DIRS[*]}; do
  FIND_EXCLUDES+="-path ./$e -prune -o "
done

cd $SOURCE_DIR
CXX_FILES=$(find . $FIND_EXCLUDES -type f \( $CXX_FIND_NAMES \) -print)
CMAKE_FILES=$(find . $FIND_EXCLUDES -type f \( $CMAKE_FIND_NAMES \) -print)

if [ -n "$DRYRUN" ]; then
  # If the clang-format xml has no replacement entries, all files are formatted.
  if !(clang-format -style=file --output-replacements-xml $CXX_FILES |
       grep -q "<replacement ") && cmake-format --check $CMAKE_FILES; then
    # print in bold-face green
    echo -e "\e[1m\e[32mGreat job, all files are properly formatted!\e[0m"
    exit 0;
  else
    # Print in bold-face red
    echo -e "\e[1m\e[31mAt least one file is not properly formatted!\e[0m"
    echo -e "\e[1m\e[31mRun $0 for formatting all files!\e[0m"
    exit 1;
  fi
else
  clang-format -i -style=file $CXX_FILES
  cmake-format -i $CMAKE_FILES
  # print in bold-face
  echo -e "\e[1mSuccessfully formatted all files.\e[0m"
fi
