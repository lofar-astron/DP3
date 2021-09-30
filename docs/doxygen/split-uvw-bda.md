@section uvw_splitting_bda UVW splitting of BDA expanded data

## Problem description

To save computations the [OnePredict](@ref dp3::steps::OnePredict) step computes station based uvw-coordinates. The phasors can then be computed per station instead of per baseline.
The computation a phasor involves a sincos evaluation, which is expensive.

To compute station based uvw-coordinates only part of the uvw-coordinates are needed.
That is because there is redundant information in the uvw-coordinates.
The baseline A-B can be composed out of other baselines if these from a route from A to B.
A minimal set of baselines that reach all station form a spanning tree.

It can occur that a set of baselines connects disjoint groups of stations.
In that case there multiple disjoint spanning trees can be formed.

Although for regular data uvw-coordinates per baseline is truly redundant, is this
not entirely true for BDA-ed data. Consider a simple linear array

![Linear array](docs/doxygen/images/linear_array.svg)

The long baseline 0-3 can be composed of the short baselines 0-1, 1-2 and 2-3.
The short baselines have been averaged to a lower resolution than the long baseline.
When the spanning tree is formed by the short baselines, the long baseline is 
reconstructed from the short ones. This leads to a loss of information.

## Solution

To make sure no information is lost when computing station based uvw-coordinates from
baseline based uvw-coordinates, one needs to make sure the uvw-coordinates of the
longest baseline that reaches a station is used for the computation of that station's
uvw-coordinates.

The function [nsetupSplitUVW](@ref dp3::base::nsetupSplitUVW) returns a set of indices that
the function [nsplitUVW](@ref dp3::base::nsplitUVW) uses to compute station based
uvw-coordinates from baseline based uvw-coordinates. It starts at a reference station
and then extends the set of known stations by adding baselines connecting a known and
an unknown station. That way the baselines are found in the correct order to 
compute the uvw-coordinates.

The function [nsetupSplitUVW](@ref dp3::base::nsetupSplitUVW) does not take the
baseline lenght into account.
To support BDA data another [nsetupSplitUVW](@ref dp3::base::nsetupSplitUVW) has been added
that accepts station positions as extra input parameters, and uses that information to
select the long baselines for station based uvw-coordinate computation.
To do that it needs to build up the set in a different way, greedily selecting the 
longest baselines, even if they connect two unkown stations. The order the baselines are 
found is not suitable to compute the uvw-coordinates. 
Therefore the regular [nsetupSplitUVW](@ref dp3::base::nsetupSplitUVW) is
called with the selected baselines to get the indices in the right order.

## Algorithm

The algorithm iterates over the baselines in order of ascending length.
It keeps track of the selected baselines and selected stations, and in what group the stations are.

When a baseline is considered for selection, four different cases can be distinguished:

Case 1.
![Baseline connects two unkown stations](docs/doxygen/images/add_baseline0.svg)
The baseline is selected and a new group is formed containing the two stations.

Case 2.
![Baseline connects a known and an unknown station](docs/doxygen/images/add_baseline1.svg)
The baseline is selected and the unkown station is added to the group of the known station.

Case 3.
![Baseline connects two known stations in the same group](docs/doxygen/images/add_baseline2.svg)
The baseline adds no new information, and is thus not selected.

![Baseline connects two known stations from different groups](docs/doxygen/images/add_baseline3.svg)
The baseline is selected and the two groups are joined

## Results

Below an example is shown where it matters whether nsetupSplitUVW uses baseline length information,
or not. The graphs compare visibilities that were obtained by either first predicting and then
bda-averaging, or the other way around. Ideally they produce the same result. For short baselines
there is little difference, but for longer baselines a much better agreement is obtained
when using the uvw-coordinates of the long baselines.

Short baseline before the fix | Long baseline before the fix
:----------------------------:|:-------------------------:
![ ](docs/doxygen/images/bdapredict-bl-3-20-fail.svg) | ![ ](docs/doxygen/images/bdapredict-bl-3-45-fail.svg)

Short baseline after the fix | Long baseline after the fix
:-------------------------:|:-------------------------:
![ ](docs/doxygen/images/bdapredict-bl-3-20-success.svg) | ![Long baseline after the fix](docs/doxygen/images/bdapredict-bl-3-45-success.svg)

