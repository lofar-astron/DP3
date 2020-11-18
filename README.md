# DP3
LOFAR preprocessing software, including averaging, flagging, various kinds of calibration and more.

The DPPP documentation can be found at: https://www.astron.nl/citt/DP3

This repository is a continuation of the one at svn.astron.nl/LOFAR. In particular, it has branched off at LOFAR Release 3.2 (Sept 2018). The version of DP3 that is in the ASTRON repository is no longer maintained.

## Installation
Some non-standard dependencies of this project are: armadillo, boost, boost-python, casacore, hdf5, aoflagger, and EveryBeam. See the Dockerfiles [`docker/ubuntu_20_04_base`](docker/ubuntu_20_04_base) and [`docker/ubuntu_20_04_lofar`](docker/ubuntu_20_04_lofar) as examples.

Typical installation commands:
```
mkdir build
cd build
cmake ..
make -j4
make install
```
