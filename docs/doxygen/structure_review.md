# DP3 Structure

The goal of this document is to review the current structure of DP3 and identify potential impediments for adding more functionality in the future and identify possible improvements to its structure to mitigate or resolve those impediments.

We focus on the impediments for allowing multiple directions to be streamed across the application. This currently cannot be done in the DDECal step. To enable that, we need multiple predict types, e.g. a combination of Predict, IDGPredict and model column(s). The wanted data flow in DDECal is given in the image below where the current DDECal implementation is highlighted in yellow. In this example there are two streams for two directions (Predict --> Applycal and IDGPredict --> UVWFlag). The current implementation only allows a single data stream for all directions.

![Detailed process of a DDECal pipeline](docs/doxygen/images/ddecal.png)

## Identified problems
We have identified the following impediments in the structure of DP3 that prevent the DDECal step to enable multiple predict types.
1. A DPStep can only have one input and one output. For example, DDECal requires inputs from three previous step, but can only receive a Buffer from one previous step. This is currently done by having the entire data stream implemented and orchestrated in the DDECal step (except the MSReader and MSWriter).
2. A DPStep requires a Buffer as input, but it does not always need the (entire) Buffer. For example, the IDGPredict step only needs the UVW and not the entire buffer. Additionally, the second ApplyCal step needs an hdf5 file which cannot be given, but since all processing is done in the DDECal step, the hdf5 file that is needed in that step is already in available.

## Improvements
1. To solve problem 1, we need to implement internal streams. The process function will no longer have a DPBuffer as input but rather a vector of DPBuffer. An entity orchestrates the process and makes sure that the data flows as described.
2. DPBuffer will be split into multiple types. For example, UVW will be a buffer on its own. Using a vector as input for the process function, the DPSteps have access to the information they need and nothing more or less. It allows any kind of data type to be streamed in DP3.

## Implementation plan
During the review sessions we concluded that to enable multiple predict types in DDECal not much needs to be done. Therefore, we have decided to not refactor the DP3 structure significantly because we can implement the required feature without changing the core structure of DP3.

We have identified 2 key refactors that fulfil the required use cases:
1. Refactor the creation of ApplyCal steps in Prediction steps. DDECal should be responsible for the creation of ApplyCal, not Predict. This should therefore be done in the DDECal step and note in the Predict step.
2. Factor out IDGPredict as separate step so that it can be used outside of DDECal.
3. Refactor the nested DDECal model pipeline so that it can be defined separately for every direction instead of one to rule them all. The metadata must not change in any of the nested model steps of DDECal. This requires a refactors in the process function of DDECal.cc.

From a user perspective we propose the following additional parset definitions for DDECal:
- `ddecal.modelnextsteps`: Defines the default steps for the models that need to be applied to every direction before applying DDECal. Each step can have its own parset entries.
- `ddecal.modelnextsteps.<direction>`: Allows steps to be defined for a specific direction. Overrides the default steps for this direction.

An example of a parset for DDECal using different models for different directions is as follows:
```
ddecal.images = [image1.fits, image2.fits]  
ddecal.regions = ds9.reg  
ddecal.sourcedb = ATeam.sourcedb  
ddecal.directions = [CasA, CygA] # subselection of patches from sourcedb

ddecal.modelnextsteps = [applycal3]  
ddecal.modelnextsteps.CasA = [applycal4, uvwflag3, applybeam3]  
ddecal.modelnextsteps.CygA = [applycal5, applybeam4]  

applycal5.h5parm = known_corruptions.h5 # Will automagically get right direction
```
