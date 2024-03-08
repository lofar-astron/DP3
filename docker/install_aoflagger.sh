#!/bin/bash

# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script to install AOFlagger from source.

set -euo pipefail

AOFLAGGER_VERSION=$1
PYTHON_VERSION=$2

if [ -z "$PYTHON_VERSION" ]; then
  echo "Usage: $0 <aoflagger version> <python version>"
  exit 1
fi

echo "Installing development libraries with yum"
yum -y install libpng-devel zlib-devel

pushd /tmp

echo "Downloading & unpacking AOFlagger ${AOFLAGGER_VERSION}"
if [ "${AOFLAGGER_VERSION}" != "3.2.0" ]; then
  echo $0: Please update the link for downloading AOFlagger ${AOFLAGGER_VERSION}.
  exit 1
fi
curl -fsSL -o aoflagger-v${AOFLAGGER_VERSION}.bz2 "https://gitlab.com/aroffringa/aoflagger/-/package_files/38600488/download"
tar -xjf aoflagger-v${AOFLAGGER_VERSION}.bz2

pushd aoflagger-v${AOFLAGGER_VERSION}

# Wheels should not actually link to libpython (see py310_wheel.docker).
sed -i '/find_package(PythonLibs 3 REQUIRED)/d' CMakeLists.txt

# Ensure AOFlagger uses the correct python version and finds it before pybind11 does.
sed -i "s=# Include aocommon/pybind11 headers=find_package(PythonInterp ${PYTHON_VERSION} EXACT REQUIRED)=" CMakeLists.txt

# The patches below are already in the AOFlagger master, and should be removed
# for AOFlagger > 3.2.0. See MR 204 and MR 205.
if [ "${AOFLAGGER_VERSION}" != "3.2.0" ]; then
  echo $0: Please remove the patches below and this check!
  exit 1
fi

# Using a component for the library allows installing it separately.
sed -i 's=install(TARGETS aoflagger-lib DESTINATION lib)=install(TARGETS aoflagger-lib DESTINATION lib COMPONENT core)=' CMakeLists.txt
sed -i 's=install(FILES interface/aoflagger.h DESTINATION include)=install(FILES interface/aoflagger.h DESTINATION include COMPONENT core)=' CMakeLists.txt
# This line updates two install() commands.
sed -i 's=DESTINATION share/aoflagger=COMPONENT core DESTINATION share/aoflagger=' CMakeLists.txt

# Remove H5P_DEFAULT arguments, which cause compilation errors.
sed -i -E -z 's=,\n? *H5P_DEFAULT\);=);=g' imagesets/sdhdfimageset.cpp

# Make AOFlagger find the strategies relative to the library directory.
patch -p0 <<EOF
--- lua/telescopefile.cpp	2022-11-21 16:49:58.806867642 +0100
+++ lua/telescopefile.cpp.new	2022-11-21 16:49:54.120925348 +0100
@@ -1,5 +1,7 @@
 #include "telescopefile.h"

+#include <dlfcn.h>
+
 #include "../imagesets/bhfitsimageset.h"
 #include "../imagesets/filterbankset.h"
 #include "../imagesets/fitsimageset.h"
@@ -116,9 +118,16 @@

   std::filesystem::path search;

-  search = std::filesystem::path(AOFLAGGER_INSTALL_PATH) /
-           "share/aoflagger/strategies" / filename;
-  if (std::filesystem::exists(search)) return search.string();
+  // Try using the directory of the AOFlagger library.
+  // When bundled as a python binary wheel, the strategies are there.
+  Dl_info dl_info;
+  if (dladdr(reinterpret_cast<const void*>(&TelescopeFile::FindStrategy),
+             &dl_info)) {
+    std::filesystem::path aoflagger_library_path(dl_info.dli_fname);
+    search = aoflagger_library_path.remove_filename() / "aoflagger/strategies" /
+             filename;
+    if (std::filesystem::exists(search)) return search.string();
+  }

   if (!argv0.empty()) {
     const std::filesystem::path root =
EOF

echo "Configuring, building & installing AOFlagger ${AOFLAGGER_VERSION}"
mkdir build
cd build
#-DFFTW3_LIB=${FFTW_DIR}/lib/libfftw3.so \
cmake \
  -DENABLE_GUI=False \
  ..
make -j${THREADS} aoflagger-lib
cmake --install . --component core
popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary AOFlagger files"
rm -r aoflagger-v${AOFLAGGER_VERSION}
rm aoflagger-v${AOFLAGGER_VERSION}.bz2

popd

echo "Cleaning up development libraries"
yum -y erase libpng-devel zlib-devel
