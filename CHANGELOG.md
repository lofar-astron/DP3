# DP3 Changelog

## [6.0] - 2023-08-11

### New features

- Optional support for Sagecal (libdirac), if that is installed.
- Re-using model visibilities in multiple steps (e.g. ddecal steps) is now possible through the 'extradata' feature.
- Python bindings are thoroughly refactored; entire DP3 pipelines can now be constructed from Python.
- DP3 now has some steps with increased performance specifically for the AARTFAAC instrument.

### Improvements

- Under the hood some casacore structures were replaced by xtensor structures, which may reduce memory usage.
- Deprecate the use of `FULL_RES_FLAGS`. These are not written anymore.
- DP3 now requires EveryBeam v0.5.1
- Experimental use of AVX and XSIMD features to speed up various operations.

### Bug fixes

- Fix predicting sky components with negative n-value (further than 90 degrees from phase centre) - note that this bug did not affect the demix step.
