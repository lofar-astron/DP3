# DP3 Changelog

## [(Upcoming release)] - 202?-??-??

### New features

- Add `predict.usefastpredict` setting to enable the new fast predict at runtime, deprecating the `USE_FAST_PREDICT` CMake option.

### Improvements

### Bug fixes

## [6.5.1] - 2025-10-09

### Bug fixes

- Fixed compilation and updated package for BDA support.
- Fixed argument handling in Python API for `execute` and `execute_from_command_line`

## [6.5] - 2025-10-01

### New features

- Add `ddecal.antenna_averaging_factors` setting for specifying different solutions intervals per antenna.
- Add `ddecal.antenna_smoothness_factors` setting for specifying smoothness factors per antenna.
- Support `ddecal.keepmodel` and `ddecal.reusemodel` settings for `ddecal` steps that process BDA data.
- Add `predict.coefficients_path` setting, e.g., for specifying MWA beam model coefficients.
- Add `predict.outputbuffername` setting, for selecting named data buffers.
- Add `transfer` step for transfering data and flags from low- to high-resolution MSs.
- Add `combine` step for adding or subtracting named data buffers.
- Support new simulated signal compression (Sisco) in CasaCore.
- Use new fast `predict` step instead of the normal `predict` step when USE_FAST_PREDICT is passed to CMake.
- Support streaming cobalt (LOFAR) data in addition to ALMA data.

### Improvements

- Remove dependency on C compiler.
- Support numpy 2.x.
- Use 50% of memory instead of 10% for wgridder predict buffer.
- Show number of solver iterations in BDA `ddecal` step.
- `ddecal.h5parm` is now a mandatory setting. The default `instrument.h5` filename is now deprecated.
- Improve `preflagger` performance.
- Remove `ddecal.initialsolutions.gaintype` setting. Always deduce the gain type from the initial solutions.
- When using `onlypredict=True` and `keepmodel=True` in a `ddecal` step, it no longer overwrites the main data buffer. This improvement allows predicting model data using a regular `ddecal` step and reusing it in a BDA `ddecal` step.
- Show total/combined averaging factor for all baselines in `bdaaverager` step.
- In the `predict` step, replace the upsampling method for time-smearing correction by an approximation using the sinc function.
- Require EveryBeam 0.7.4, since DP3 now uses the new multi-frequency interface of EveryBeam.

### Bug fixes

- Fix and update default `msout.tilenchan` setting. The default value is now 64.
- Fix updating phase centre when writing BDA MSs.
- Fix a bug in antenna-uvw calculations causing incorrect data prediction. This bug did not manifest itself for LOFAR observations, it was observed in specific MWA tests.
- In `ddecal`, avoid printing empty lines with sub step timings.
- In the BDA implementation of `ddecal`, fix using multiple solutions per direction.
- In `bdaexpander`, avoid floating point rounding errors for the averaging factor.


## [6.4.1] - 2025-05-22

### New features
- Add a new rotationconstraint in DDECal for fitting Faraday rotation.
- Support multiple directions in the rotation and rotation-and-diagonal constraints.

### Improvements
- Add missing beam keywords to BDA data.


## [6.4] - 2025-04-14

This release completes the meta-data compression feature, which makes LOFAR MSs ~14% smaller. Compressing meta-data is off by default but can be enabled in CMake and/or in the parset. This release also adds various calibration options.

### New features
- Add `msout.uvwcompression` and `msout.antennacompression` options, for compressing
UVW and antenna values using new CasaCore storage managers.
- Add `clipper.flagallcorrelations` option, for flagging all correlations instead of a single correlation only.
- Add `ddecal.initialsolutions` option, which allows starting from existing solutions when running the solver.
- Add `ddecal.smoothness_dd_factors` option, for direction-dependent scaling of the smoothing kernel.
- Add `ddecal.smoothness_kernel_truncation` option, for disabling truncating the smoothing kernel.

### Improvements
- Update handling `msin.starttime` / `msin.starttimeslots` options. DP3 now uses a time slot from the input MS as start time, instead of using the first time slot in the MS and increasing it by an integer times the MS interval.
- Require CasaCore 3.7.1 when building DP3.
- Average extra (model) data buffers too in the BDA averager step.
- Add documentation for the `null` step.
- Build binary wheel for Python 3.13; Drop binary wheel for Python 3.7.

### Bug fixes
- Fix angular proximity clustering for models with sources near RA=0.
- Fix writing the h5 time axis in the ddecal step when the solution interval is larger than the number of time slots.
- Allow updating the input MS after running a 'predict' step.
- Fix showing ddecal settings when using the 'hybrid' solver.

## [6.3] - 2025-01-28

### New features
- Add flagtransfer step, which transfers flags from a low-resolution MS.
- Support new Casacore Stokes I storage manager.
- Add `msout.scalarflags` option, for compressing flags.

### Improvements
- The `elementmodel` parset key of the ApplyBeam and Predict step is now parsed
by EveryBeam making all element models in that library available.
The default value is changed from "hamaker" to "default" meaning that
EveryBeam selects the default element model for the telescope in
the measurement set. For a LOFAR MS that will still be "hamaker".
- Support EveryBeam 0.7.x.
- Apply per-direction weights in constraints.
- Use C++20. DP3 now requires at least GCC-10.

### Bug fixes
- Fix flag counting in UVWFlagger when using BDA.

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
