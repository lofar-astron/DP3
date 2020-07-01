# Baseline Dependent Averaging : Design

## Introduction
Baseline Dependent Averaging(BDA) aims at reducing the amount of data in a
Measurement Set(MS) by averaging measurements that provide similar information.
The following criteria determine the similarity of a measurement:
- The length of the baseline. Measurements with short baselines are more
  similar than measurements with long baselines.
- The frequency. Measurements for low frequencies are less similar than 
  for high
  frequencies. (Is this statement correct???)

## Current situation

- A [DPBuffer](@ref DP3::DPPP::DPBuffer) has a regular multidimensional
  cube structure:
  - Main layout: 
  - The length of a time step is equal for all baselines.
  - A DPBuffer contains the measurements for a single time step.
  - The amount of frequency channels is equal for all baselines.
    All channels have the same size.
- The [DPStep::process()](@ref DP3::DPPP::DPStep::process) function has a
  single [DPBuffer](@ref DP3::DPPP::DPBuffer) argument. It processes
  the measurements for a single time step.

## Requirements

DP3 should support the following:
1. Converting a regular buffer to BDA data, using a BDA algorithm.
2. Converting BDA data to a regular buffer, using an expansion algorithm.
3. Read /write BDA data from/to an MS using the Casacore library.
   The MS format should be aligned with other tools.
   (Which tools???)
4. Perform processing steps that only use the data for the current iteration:
   - Predict
   - ApplyBeam
   - ApplyCal
   - Preflagger
   - Phaseshift
   - Filter
   - Scaledata
   - Preflagger
   - UVWFlagger
5. Perform processing steps that use data from multiple iterations.
   The steps store data items of the current iteration for
   use in later iterations. These steps include:
   - Calibrate
   - DDECal
6. The DPPP input should allow specifying the BDA operations above.

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

### DPStep structure

DPSteps can support BDA data using a new process() function for BDA data.
The default implementatation for this function will throw an exception.
The existing process() function for regular data, which is now abstract,
should also get a default implementation that throws an exception. In code:

    class DPStep {
    public:
      virtual bool process(DPBuffer) { throw exception; } // Process regular data.
      virtual bool process(BDABuffer) { throw exception; } // Process BDA data.
    };

With these two process functions, a DPStep can support
regular and/or BDA data: Existing DPSteps, which only support regular
data will remain working. Adding support for BDA data can be done
incrementally, one DPStep at a time.

A disadvantage of this approach is that DPSteps will become larger and
more complex: They may get two modes of operation, which
violates the principle of having a single purpose for each class.
However, if a DPStep becomes too complex, this design allows
splitting it into a regular variant and a BDA variant for reducing complexity

### BDA data structure

A crucial point in this design is the BDABuffer which holds BDA data.
The smallest data element in this structure is a measurement row, which
contains the measurements for all polarisations in all channels for a single
baseline and a single time step. Since BDA averages data for multiple channels,
different rows may have different amounts of channels. In code:

    struct Row {
      double time;
      double exposure;
      Matrix<complex> data(nChannels, nPolarisations);
      Matrix<bool> flags(nChannels, nPolarisations);
      Matrix<float> weights(nChannels, nPolarisations);
      Matrix<bool> fullResFlags(nChannels, nPolarisations);
    };

For reducing overhead, a BDABuffer should not contain a single Row. Using
a variable-sized vector of rows provides a flexible structure: The BDABuffer
size can be adjusted to the number of channels in the Rows.
Combining and splitting buffers is also possible. In code:

    class BDABuffer
    {
      vector<Row> itsRows; // Buffer contains a variable number of rows.
    };

This approach still has a major drawback: Each Matrix object in each row
performs its own memory allocation, which may yield high overhead.
A BDABuffer can mitigate this overhead using a shared memory pool for all Rows.

### Steps that use data from multiple iterations.

Steps that use data from multiple iterations have an internal storage
that holds the data from previous iterations. For regular data, this internal
storage holds a fixed number of time steps.

For BDA data, the exposure duration of a time step is different for each Row.
For gathering BDA data for a time interval, an additional BDAIntervalBuffer
class will be necessary. It should provide the following functionality
to a step:
- Set an initial time interval as the current time interval.
- Store data from BDABuffers.
- Detect if all relevant data for the current time interval is present.
- Provide the data for the current time interval. The result should take
  the BDABuffer start
  time and exposure duration into account, since they may not match the start
  time and duration of the requested time interval.
- Advance the current time interval. The BDAIntervalBuffer should then
  discard stored data that is no longer necessary.
