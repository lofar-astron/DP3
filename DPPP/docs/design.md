# Design

## DPPP: a framework for streaming radio interferometric pipelines
DPPP (DP3, short for Default PreProcessing Pipeline) is a framework
for efficiently processing time-ordered visibilities,
e.g. freshly correlated data. It was developed for the LOFAR
telescope, and is used for preprocessing its imaging pipelines.
DPPP documentation is found at the [LOFAR Wiki](https://goo.gl/2uAaN9).

## Optimised for IO-heavy tasks
When DPPP reads a chunk of data, it calls all the steps to do
all processing possible on this data. The last step will then
write out the processed data. In this way, the data is read
and written only once, and memory consumption is minimal.

This is specially useful for IO-limited operations, where lots of
data is involved and not too much computations are done.

## Measurement Set standard
DPPP adheres to the Measurement Set (MS) standard, so it
can be used with data from all radio telescopes. DPPP can
also create a new MS. When an appropriate input step is in
place (e.g. a streaming connection to a correlator), MSs can
be written from scratch.

DPPP supports only regularly shaped MSs, meaning that every
time slot must have the same number of baselines and antennas.
Only one spectral window can be used.

## Fail early
DPPP runs a pipeline in clearly defined stages, to make spotting
errors as easy as possible.
* **Initialise** This connects all the processing steps together. This step is data independent.
* **Tune At** this stage only metadata, such as the number of channels, is read, and steps are adjusted to this.
* **Process** Data is piped through the steps in buffers of time. The buffer size is defined by the steps.

## Calibration
The DPPP step [GainCal](@ref DP3::DPPP::GainCal) implements many variants of direction
independent calibration. It uses StefCal, with many
extra features, such as fitting a function over frequency.
Also full-Jones calibration is supported. Results are written
to ParmDB; H5Parm support is in development.

## AOFlagger
A shallow wrapper exists that calls AOFlagger from DPPP.
This makes it possible to flag with AOFlagger and immediately
afterwards average the data, without writing the
data to disk in between.


## Built-in steps
All the built-in steps derive from the [DPStep](@ref DP3::DPPP::DPStep) class. Alternatively, an overview can be found at the [LOFAR Wiki](https://goo.gl/2uAaN9). Some of the core built-in steps are given in the table below.
|              |                                                |
|--------------|------------------------------------------------|
| [aoflag](@ref DP3::DPPP::AOFlaggerStep)       | Automatic flagging in time/freq windows        |
| [applybeam](@ref DP3::DPPP::DemixInfo)    | Apply the LOFAR beam model                     |
| [applycal](@ref DP3::DPPP::ApplyCal)     | Apply instrument corrections                   |
| [average](@ref DP3::DPPP::Averager)      | Average data in time and/or frequency          |
| [demix](@ref DP3::DPPP::Demixer)        | Calibrate and subtract strong sources          |
| [filter](@ref DP3::DPPP::Filter)       | Leave out specified data in further processing |
| [gaincal](@ref DP3::DPPP::GainCal)      | Calibrate against a given skymodel             |
| [phaseshift](@ref DP3::DPPP::PhaseShift)   | Shift data to a different phase center         |
| [scaledata](@ref DP3::DPPP::ScaleData)    | Scale data with a polynomial in frequency      |
| [stationadder](@ref DP3::DPPP::StationAdder) | Add signals from stations (e.g. superterp)     |
| [predict](@ref DP3::DPPP::Predict)      | Predict visibilities of a given skymodel       |
| [preflag](@ref DP3::DPPP::PreFlagger)      | Flag given baselines, antennas, channels, etc. |

## User-defined steps
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
key-value format. Optionally, some parameters can be overridden
on the command line.

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
