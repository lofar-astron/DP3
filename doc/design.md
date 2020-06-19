# Design

## DPPP: a framework for streaming radio interferometric pipelines
DPPP (DP3, short for Default PreProcessing Pipeline) is a framework
for efficiently processing time-ordered visibilities,
e.g. freshly correlated data. It was developed for the LOFAR
telescope, and is used for preprocessing its imaging pipelines.
DPPP user documentation is found at the [LOFAR Wiki](https://goo.gl/2uAaN9).

When DPPP reads a chunk of data, it calls all the steps to do
all processing possible on this data. The last step will then
write out the processed data. In this way, the data is read
and written only once.

![High level overview of DPPP](doc/images/diagram.png)

DPPP runs a pipeline in clearly defined stages, to make spotting
errors as easy as possible.
* **Initialise** This connects all the processing steps together. This step is data independent.
* **Tune At** this stage only metadata, such as the number of channels, is read, and steps are adjusted to this.
* **Process** Data is piped through the steps in buffers of time. The buffer size is defined by the steps.

## Measurement Set standard
DPPP adheres to the Measurement Set (MS) standard, so it
can be used with data from all radio telescopes. DPPP can
also create a new MS. When an appropriate input step is in
place (e.g. a streaming connection to a correlator), MSs can
be written from scratch.

DPPP supports only regularly shaped MSs, meaning that every
time slot must have the same number of baselines and antennas.
Only one spectral window can be used.

## Built-in steps
DPPP contains built-in steps that are linked together and process one time slot at a time and is therefore suitable for large measurement sets.
All the built-in steps derive from the [DPStep](@ref DP3::DPPP::DPStep) class.

[DPRun](@ref DP3::DPPP::DPRun) initialises and orchestrates the configured steps.
It begins with setting the required info for all the first step and its next step(s).
Next, each buffer (containing a single time slot) is processed by the first step.
When processed, it invokes the process function of the next step.
Once all buffers are processed it finishes the processing of the first step and subsequent steps.
It wraps it with adding some data to the MeasurementSet written/updated which is started at the last step.
The figure below gives a graphical overview of the [DPStep](@ref DP3::DPPP::DPStep) flow.

![Process flow of DPPP](doc/images/flow.png)

### Calibration
The DPPP step [GainCal](@ref DP3::DPPP::GainCal) implements many variants of direction
independent calibration. It uses StefCal, with many
extra features, such as fitting a function over frequency.
Also full-Jones calibration is supported. Results are written
to ParmDB; H5Parm support is in development.

### AOFlagger
A shallow wrapper exists that calls AOFlagger from DPPP.
This makes it possible to flag with AOFlagger and immediately
afterwards average the data, without writing the
data to disk in between.

### User-defined steps
Apart from the built-in steps, also custom steps can
be defined, following the existing interface description.
Recompiling DPPP is not necessary: steps are found as shared
libraries.

A special custom step is the [Python step](@ref DP3::DPPP::PythonStep), which finds a
python module and applies that. This makes it very easy to
implement custom steps even for non C++-wizards. An
example Python step is bundled with the LOFAR software.

## Example reduction
Arguments are given to DPPP as a 'parameter set', a simple
key-value format. Optionally, parameters can be overridden
on the command line. The example also shows that you can 
use the out step to store intermediate results.

    msin = [bla.MS, bla2.MS]
    steps = [aoflag, avg1, out1, cal, out2, avg2 ]
    avg1.type = average
    avg1.timestep = 4 # Average down every 4 timeslots
    out1.type = out
    out1.name = myfile_averaged.MS
    cal.type = calibrate
    cal.skymodel = 3C196.skymodel
    cal.applysolution = true
    cal.parmdb = instrument.parmdb
    out2.type = out
    out2.name = .
    out2.datacolumn = CORRECTED_DATA
    avg2.type = average
    avg2.timestep = 2 # average down some more
    avg2.freqstep = 16
    msout = bla_combined.MS
