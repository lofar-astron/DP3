# Baseline Dependent Averaging : Design

## Introduction
Baseline Dependent Averaging(BDA) aims at reducing the amount of data in a
Measurement Set(MS) by averaging measurements with low significance.
The following criteria determine the significance of a measurement:
- The length of the baseline. Measurements with short baselines are less
  significant than measurements with long baselines.
- The frequency. Low frequencies are more significant than high frequencies.
  (Is this statement correct???)

## Current situation

- The [DPStep::process()](@ref DP3::DPPP::DPStep::process) function has a
  single [DPBuffer](@ref DP3::DPPP::DPBuffer) argument.
- A [DPBuffer](@ref DP3::DPPP::DPBuffer) has a regular multidimensional
  cube structure:
  - The lenght of a time step is equal for all baselines.
  - The amount of frequency channels is equal for all baselines.
    All channels have the same size.

## Requirements

DPPP should support the following:
1. Converting a regular buffer to BDA data, using a BDA algorithm.
2. Converting BDA data to a regular buffer, using an expansion algorithm.
3. Read /write BDA data from/to an MS using the Casacore library.
   The MS format should be aligned with other tools.
   (Which tools???)
4. Perform the following processing steps on BDA data:
   - Calibrate
   - Predict
   - ApplyBeam
   - ApplyCal
   - DDECal
   - Preflagger
   - Phaseshift
   - Filter
   - Scaledata
   - Preflagger
   - UVWFlagger
5. The DPPP input should allow specifying the BDA operations above.

### Non-requirements

The following processing steps explicitly do have to support BDA data,
since it is impossible:
- Average
- AOFlagger
- Interpolate
- StationAdder

## Design

The current [DPStep](@ref DP3::DPPP::DPStep) class can be reused for supporting
BDA operations, by extending the process() method with support for BDA data.
This approach supports all requirements: A step can perform a BDA process()
call sending BDA data to its next step and implement the BDA process()
method for receiving BDA data from its previous step.

With this approach, the DPPP input format can remain the same: The steps
will receive and send data in the requested format.

### BDA data structure

A crucial point in this design is the data structure for holding BDA data.
There are various alternatives for this data structure:
- foo
- bar
