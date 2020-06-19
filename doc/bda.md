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

DPSteps can support BDA data in many ways:

- Keep the existing interface, but adjust DPBuffer, so it can also
  contain BDA data.
  - Advantages:
    - There is no need for a new process() method in the DPStep interface.
  - Disadvantages:
    - All existing code that uses DPBuffers has to check if a DPBuffer
      has regular data or BDA data.
    - The process() method will probably call new processRegular() or
      processBda() methods, depending on the DPBuffer, so new methods will
      be introduced anyway.
  - Cost: High.
    - A lot of existing code must be modified. 
  - Risk: High.
    - Forgetting updating some existing code is easy.

- Introduce a new DPBDABuffer for BDA data, and augment the DPStep interface
  with a new process() function for this new buffer type.
  - Advantages:
    - Existing code that uses DPBuffers requires no changes.
    - DPStep classes can be updated incrementally, one at a time.
    - Regular steps and BDA steps can be easily mixed.
  - Disadvantages:
    - The DPStep interface requires a new process() method for 
      DPBDABuffers.
    - Existing DPStep classes will become larger.
    - A single DPStep class may get two modes of operation, which
      violates the principle of having a single purpose for each class.
  - Cost: Medium.
    - Most DPStep classes will require a new process() for DPBDABuffers.
  - Risk: Low:
    - Changes do not touch existing code much.
    - If a class becomes too complex, splitting it into a
      regular variant and a BDA variant is possible.

- Introduce a new DPBDABuffer and use it in a new DPBDAStep interface.
  - Advantages:
    - Existing code will remain working.
  - Disadvantages:
    - Mixing DPStep and DPBDAStep classes, in the conversion steps, will 
      be non-trivial.
    - Reusing code between DPSteps and DPBDASteps that perform similar
      operations will be non-trivial.
    - The step set-up code needs to know upfront if it should create DPSteps
      or DPBDASteps. Towards the user, regular steps and bda steps may have a
      different name, which is not needed with the other approaches.
  - Cost: High.
    - Even existing code may have to be modified, for avoiding duplicate code.
    - Maintenance cost will be high, as similar operations are provided by
      different classes.
  - Risk: High.
    - Non-trivial solutions are required.

Based on this overview the best solution is augmenting the DPStep interface.

### BDA data structure

A crucial point in this design is the DPBDABuffer which holds BDA data.
The smallest data element in this structure is a Measurement Line, which
contains the measurements for all polarisations in all channels for a single
baseline and a single time step. Since BDA averages data for multiple channels,
different MLs may have different amounts of channels. In code:

    using Line = struct {
      Time startTime;
      Duration duration;
      Matrix<complex> data(nChannels, nPolarisations);
    };

There are various alternatives for the DPBDABuffer:

- 

      class DPBDABuffer
      {
        Line itsData; // Buffer contains a single Line
      };

  - Disadvantages:
    - Since the buffer has little data, the relative overhead of DPStep function
      calls, where are executed for each buffer, will be high.
  - Cost: Low
    - The structure is simple.
  - Risk: High
    - Great risk that the overhead is too high.

-

    class DPBDABuffer
    {
      Line itsData[100]; // Buffer contains fixed number of lines.
    };

  - Disadvantages:
    - The last buffer will not be full and requires special treatment.
    - Since the number of channels may vary between Lines, the actual
      buffer size will still vary.
  - Cost: Low
    - The structure is simple.
  - Risk: Medium
    - There are special border cases.
    - The structure is not flexible: The buffer size can not be adjusted
      to the number of channels in the Lines.

-

    class DPBDABuffer
    {
      vector<Line> itsData; // Buffer contains a variable number of lines.
    };

  - Advantages:
    - The structure is flexible. Combining and splitting buffers is even possible.
    - The buffer size can be adjusted to the number of channels in the Lines.
  - Cost: Low
    - The structure is only slightly more difficult than when using a fixed
      number of lines.
    - Special border cases for the last buffer are no longer present.
      The last buffer is simply smaller.
  - Risk: Low
    - The flexibility of the structure allows using it in unforeseen cases, too.

Based on this overview the best solution is the buffer with a variable number of lines.

### Steps that use data from multiple iterations.

Steps that use data from multiple iterations have an internal storage
that holds the data from previous iterations. For regular data, this internal
storage holds a fixed number of time steps.

For BDA data, the duration of a time step is different for each Line.
For gathering BDA data for a time interval, an additional BDAIntervalBuffer
class will be necessary. It should provide the following functionality:
- Add new DPBDABuffers to the BDAIntervalBuffers. The BDAIntervalBuffer
  should ideally only store references to data from DPBDABuffers instead of
  copying the data.
- Query if the BDAIntervalBuffer has all relevant data for a time interval.
- Query data for a time interval. The result should take the DPBDABuffer start
  time and duration into account, since they may not match the start time
  and duration of the query interval.
- Discard data for previous time intervals.

If the BDAIntervalBuffer needs to have step-specific routines, using a
(derived) class for a specific step always remains an option.

