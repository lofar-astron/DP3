#!/bin/bash
# Do check on clang-format, script is inspired on:
# - https://dev.to/10xlearner/formatting-cpp-c-javascript-and-other-stuff-2pof
# - https://gitlab.freedesktop.org/monado/monado/-/blob/master/scripts/format-and-spellcheck.sh
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl
#
# To hook this script to pre-commit include the line "./scripts/clang-format-check.sh" to .git/hooks/pre-commit
# and make sure pre-commit is an executable shell script.

programname=$0
set -e 

# Function to print usage of this script
usage () {
    s=$(printf "%-50s" "=")
    echo "${s// /=}"
    echo "Script to check whether run-clang-format.sh leads to any git diffs. If so, return non-zero exit status."
    echo ""
    echo "Usage: $programname [-i filename/pattern] [-s skip directory]"
    echo "  -i      include file(pattern) in diff, pass wildcard arguments as literals!"
    echo "  -s      skip directory in diff. Must be specified as a path, e.g. ./external"
    echo "  -h      display help"
    exit 1
}

# Function to join array to a string
join_by () { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

# Fill the include_files and exclude_directories
include_files=()
exclude_directories=()
while getopts ":hi:s:" opt
do
    case $opt in
        h) 
        usage
        exit;;
        i) 
        include_files=("${include_files[@]}" "-i '$OPTARG'")
        ;;
        s) 
        exclude_directories=("${exclude_directories[@]}" "-s $OPTARG")
        ;;
        \? ) 
        echo "Execute $programname -h to get usage of this script."
        exit;;
    esac
done

PATCH_NAME=clang-fixes.diff

echo -e "\e[1mRunning clang-format...\e[0m"
echo

# Include files?
if [ ${#include_files[@]} -eq 0 ]; then
    INC_FILE=""
else
    INC_FILE+=$(join_by " " "${include_files[@]}")
fi

# Exclude files or dirs?
if [ ${#exclude_directories[@]} -eq 0 ]; then
    EX_DIR=""
else
    EX_DIR+=$(join_by " " "${exclude_directories[@]}")
fi

# Run clang-format from run-clang-format.sh
SCRIPT_PATH=$(dirname "$0")
$SCRIPT_PATH/run-clang-format.sh $EX_DIR $INC_FILE

# Can't use tee because it hides the exit code
if git diff --patch --exit-code > $PATCH_NAME; then
    echo
    # print in bold-face green
    echo -e "\e[1m\e[32mGreat job, clang-format changed nothing!\e[0m"
    rm -f $PATCH_NAME
    exit 0;
else
    echo
    # Print in bold-face red
    echo -e "\e[1m\e[31mClang-format made at least one change. Run clang-format before pushing!\e[0m"
    rm -f $PATCH_NAME
    exit 1;
fi

