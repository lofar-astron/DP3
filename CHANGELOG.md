# DP3 Changelog

## [6.1] - 2024-06-18

- Support EveryBeam 0.6 (in addition to EveryBeam 0.5.8)
- Faster beam predict by using cached time directions
- Add support for dish telescopes like SKA-mid
- Add preliminary clipper task for VLBI processing
- Reduce memory in predict step
- Read sky model only once to speed up processing
- Various AARTFAAC processing improvements
- More use of XTensor
- Some use of AVX(2) instructions to speed up tasks
- Improve threading and use of threadpool
- Fix wrong h5 time axis when solint larger than number
- Fix prediction for source models with only negative polarized sources
- Fix compilation on platforms where size_t != 64 bit
- Fix multiple baseline selection (!1255)
- Fix ApplyCal with full-Jones correction
- Various small bugs and code improvements

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
