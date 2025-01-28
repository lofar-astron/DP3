#!/usr/bin/env python
# coding: utf-8

# Copyright (c) 2016 The Pybind Development Team, All rights reserved.

# Setup.py that calls CMake. Largely taken from https://github.com/pybind/cmake_example/

import os
import re
import shutil
import subprocess
import sys
from glob import glob

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name))
        )

        # required for auto-detection & inclusion of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        debug = (
            int(os.environ.get("DEBUG", 0))
            if self.debug is None
            else self.debug
        )
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        cmake_args = [
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            "-DCMAKE_CXX_FLAGS='-Wl,--unresolved-symbols=ignore-all'",
            # py310_wheel.docker creates a fake libpython.so so the BLAS check
            # binary fails due to missing python symbols.
            "-DCHECK_BLAS=Off",
        ]
        build_args = []

        # Using Ninja-build since it a) is available as a wheel and b)
        # multithreads automatically. MSVC would require all variables be
        # exported for Ninja to pick it up, which is a little tricky to do.
        # Users can override the generator with CMAKE_GENERATOR in CMake
        # 3.15+.
        if not cmake_generator or cmake_generator == "Ninja":
            try:
                import ninja  # noqa: F401

                ninja_executable_path = os.path.join(ninja.BIN_DIR, "ninja")
                cmake_args += [
                    "-GNinja",
                    f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_executable_path}",
                ]
            except ImportError:
                pass

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += [
                    "-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))
                ]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        build_temp = os.path.join(self.build_temp, ext.name)
        if not os.path.exists(build_temp):
            os.makedirs(build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=build_temp
        )

        # Move the core wheel files to 'extdir'
        # 'auditwheel' will add necessary libraries, including libDP3.so.
        shutil.move(os.path.join(build_temp, "pythondp3", "dp3"), extdir)


# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="DP3",
    version="6.3.0",  # Please keep in sync with the version in CMakeLists.txt.
    author="ASTRON",
    author_email="dijkema@astron.nl",
    description="DP3",
    long_description="LOFAR preprocessing software, including averaging, flagging, various kinds of calibration and more.",
    ext_modules=[CMakeExtension("dp3")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    entry_points={"console_scripts": ["DP3.py = dp3:entrypoint"]},
    python_requires=">=3.6",
    data_files=[
        (
            "../../DP3.libs/aoflagger/strategies",
            glob("/usr/local/share/aoflagger/strategies/*.lua"),
        ),
    ],
    url="https://dp3.readthedocs.io/",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
)
