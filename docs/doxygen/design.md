# Design

## DP3: a framework for streaming radio interferometric pipelines
DP3 (short for Default PreProcessing Pipeline) is a framework
for efficiently processing time-ordered visibilities,
e.g. freshly correlated data. It was developed for the LOFAR
telescope, and is used for preprocessing its imaging pipelines.
DP3 user documentation is found at the [Read The Docs page](../index.html).

For avoiding high memory usage, DP3 reads the input data in chunks instead of reading all input data at once.
For each chunk, it calls all the steps to do the processing on this data. The last step will then
write out the processed data.

![High level overview of DP3](docs/doxygen/images/diagram.png)

DP3 runs a pipeline in clearly defined stages, to make spotting
errors as easy as possible.
* **Initialise** This connects all the processing steps together in a linked list data structure. This stage is data independent.
* **Get fields to read** At this stage the fields required by each step (data, flags, weights, uvw) are combined, to get the information on what to read from the measurement set.
* **Tune At** This stage reads metadata, such as the number of channels, and adjusts the steps accordingly.
* **Process** This stage pipes the data through the steps, one buffer at a time. The steps define the buffer size.

## Reading input
Internally, DP3 passes all bulk data (visibilities, flags, weights, uvw) as so-called [DPBuffer](@ref dp3::base::DPBuffer)s or [BdaBuffer](@ref dp3::base::BdaBuffer)s. As output, a new MS can be written (respectively with [MSWriter](@ref dp3::steps::MSWriter) and [MSBDAWriter](@ref dp3::steps::MSBDAWriter) or the input MS can be updated with [MSUpdater](@ref dp3::steps::MSUpdater)). Updating is currently not implemented for BDA MSs, though. If an MS is updated, and some data from the input MS is not used by any step in the pipeline, the data is not read into the DPBuffer. DP3 implements this mechanism using the [getRequiredFields](@ref dp3::steps::Step::getRequiredFields) and [getProvidedFields](@ref dp3::steps::Step::getProvidedFields) functions of each Step.
If the first step is an MsReader and a field is not required, it will not read it into the buffer.
The RequiredFields can also hold metadata about which things need to be read into the metadata, which is stored in the [DPInfo](@ref dp3::base::DPInfo) object. The same logic holds: if a field is not required by the pipeline, the MsReader will not read it into the DPInfo object. An example of an optional metadata object is the beam information.

## Measurement Set standard
DP3 adheres to the Measurement Set (MS) standard, which is the common format for radio telescope data. DP3 can
also create a new MS. When an appropriate input step is in
place (e.g. a streaming connection to a correlator), MSs can
be written from scratch.

DP3 supports baseline-dependent averaged (BDA) MSs, though not all steps support BDA data.
For regular MSs the following conditions hold:
- Every time slot has the same number of antennas.
- Time slots have the same duration and are equally spaced (missing time slots are allowed).
- Every integration has the same number of correlations and channels.
- Only one spectral window can be handled.

## Built-in steps
DP3 contains built-in steps that are linked together in a linked list, where each step processes one time slot at a time. This way, DP3 is suitable for large measurement sets.
All the built-in steps derive from the [Step](@ref dp3::steps::Step) class.
Its inheritance diagram shows all step types, including the built-in steps.

[DPRun](@ref dp3::base::DPRun) initialises and orchestrates the configured steps.
It begins with setting the required info for the first step and its next step(s).
Next, each buffer (containing a single time slot) is processed by the first step.
When processed, it invokes the process function of the next step.
Once all buffers are processed, DP3 calls the `finish` function of the first step and subsequent steps.
This function flushes all remaining data from the steps to their output channel(s).
Finally, DP3 calls the `addToMS` function on each step for updating the metadata of the MS.
The figure below gives a graphical overview of the [Step](@ref dp3::steps::Step) flow.

![Process flow of DP3](docs/doxygen/images/flow.png)

### Calibration
The DP3 steps [GainCal](@ref dp3::steps::GainCal) and [DDECal](@ref dp3::steps::DDECal) implement many variants of direction
independent and direction dependent calibration. GainCal predates DDECal and can only be used for direction-independent calibration. DDECal is intended for direction-dependent calibration, but it can also be used for direction-independent, because it contains some features that have not been backported to GainCal.
Calibration results are stored as [H5Parm](https://github.com/revoltek/losoto/wiki/H5parm-specifications) files. GainCal also supports writing to the deprecated [ParmDB](https://www.astron.nl/lofarwiki/doku.php?id=public:user_software:documentation:makesourcedb) format.

### User-defined steps
Apart from the built-in steps, also custom steps can
be defined, following the existing interface description.
Recompiling DP3 is not necessary: steps are found as shared
libraries.

A special custom step is the [Python step](@ref dp3::pythondp3::PyStep), which uses a step defined in a Python module.
This should make it easy to implement custom steps even for non C++-wizards. An
example Python step is included in DP3.

## Example reduction
Arguments are given to DP3 as a 'parameter set' or 'parset', a simple
key-value format. Optionally, parameters can be overridden
on the command line. The example also shows that you can
use the out step to store intermediate results.

    msin = [bla.MS, bla2.MS]

    steps = [aoflag, avg1, out1, cal, out2, avg2 ]
    avg1.type = average
    avg1.timestep = 4 # Average down every 4 time slots

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
