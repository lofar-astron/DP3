# DP3 Changelog

## Upcoming release

- Add flagtransfer step, which transfers flags from a low-resolution MS.

## [6.2.2] - 2024-11-08

### Improvements

- Remove the executable `__DP3_from_pip__` from python binary wheel, replace it
  with a python command-line interface DP3.py which behaves like DP3.

## [6.2.1] - 2024-11-06

### New features

- Allow wildcards in DDECal directions.

### Improvements

- Many internal quality improvements to the input step.
- Optimize threading in solvers.
- Reduce beam evaluations.

### Bug fixes

- Fix the `__DP3_from_pip__` application in binary wheels.

## [6.2] - 2024-08-29

### New features

- Support reading extra data columns from the input Measurement Set(s).
- Support extra (model) data columns in applybeam step.
- Add wgridderpredict step, which uses ducc wgridder for image based prediction.
- Add rotationdiagonalmode setting to ddecal step.
- Add smoothnessspectralexponent setting to ddecal step.
- Add usedualvisibilites setting to ddecal step, for using XX/YY visibilities only.
- Support a restricted LBFGS solution range in ddecal and demixer.
- Add skipstations setting to applybeam step.
- Support DD intervals in bda ddecal step.
- Add wscleanwriter step which writes data in WSClean reordered format.

### Improvements

- Interface change: Remove start channel from DPInfo constructor, also in Python.
- Format python files using isort.
- Use Logger for all output. Do not allow std::cout / std::cerr anymore.

### Bug fixes

- Fix using msin.startchan / filter.startchan when updating a Measurement Set.

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
