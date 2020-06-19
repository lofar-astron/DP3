# DP3
LOFAR preprocessing software, including averaging, flagging, various kinds of calibration and more.

The DPPP documentation can be found at: https://www.astron.nl/lofarwiki/doku.php?id=public:user_software:documentation:ndppp

This repository is a continuation of the one at svn.astron.nl/LOFAR. In particular, it has branched off at LOFAR Release 3.2 (Sept 2018). The version of DP3 that is in the ASTRON repository is no longer maintained.

## Installation
Some non-standard dependencies of this project are: armadillo, boost, boost-python, casacore, hdf5, aoflagger, and the LOFAR beam model libstationresponse. These can be installed using KERN, see the file .travis/Dockerfile in this repository for an example installation.

Typical installation commands:
```
mkdir build
cd build
cmake ..
make -j4
make install
```
