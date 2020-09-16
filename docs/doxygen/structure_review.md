# DP3 Structure

The goals of this document are as follows:
- Review the current structure of DP3.
- Identify potential impediments for adding more functionality in the future.
- Identify possible improvements to its structure to mitigate or resolve those impediments.

We focus on the impediments related to handling multiple directions. The DDECal step currently does not suppport different input types for different directions. To enable that, we need multiple predict types, e.g. a combination of Predict, IDGPredict and model column(s). The image below describes the desired data flow in the DDECal step. The current DDECal implementation is highlighted in yellow. In this example there are two streams for two directions (Predict --> Applycal and IDGPredict --> UVWFlag). The current implementation only allows a single data stream for all directions.

![Detailed process of a DDECal pipeline](docs/doxygen/images/ddecal.png)

## Identified problems
We have identified the following impediments in the structure of DP3 that prevent the DDECal step to enable multiple predict types:
1. A DPStep can only have one input and one output buffer. For example, DDECal requires inputs for the model(s) besides the regular input, but can only receive one buffer from its previous step. This is currently done by using substeps inside DDECal for the models. On each invocation, DDECal invokes those substeps and uses their results as extra inputs.
2. A DPStep like the second ApplyCal step needs an hdf5 file which cannot be given, but since all processing is done in the DDECal step, the hdf5 file that is needed in that step is already in available.
3. Different directions cannot be specified in DDECal. Instead, the same processing is applied for every direction.

## Improvements
1. To solve problem 1, the process function should no longer have a single DPBuffer as input but rather a vector of DPBuffers. DP3 should orchestrate the process and make sure that the data flows as described.
2. For solving problem 2, we could introduce different DPBuffer types. For example, a buffer then only holds UVW data and nothing else. Using a vector as input for the process function, the DPSteps have access to the information they need and nothing more or less. DP3 then supports streaming any kind of data type.
3. Problem 3 can be solved by providing the user with a syntax to define with model is to be applied per direction. DDECal will then initialize the models and applies it to said directions.

## Implementation plan
During the review sessions we concluded that to enable multiple predict types in DDECal not much needs to be done. Therefore, we have decided to not refactor the DP3 structure significantly because we can implement the required feature without changing the core structure of DP3.

We have identified four key refactors that fulfil the required use cases:
1. Refactor the creation of ApplyCal steps in Predict steps. Currently, the DDECal step creates Predict substeps which in turn create ApplyCal subsubsteps. DDECal should be responsible for the creation of ApplyCal, not Predict, so the ApplyCal steps become substeps of DDECal instead of being subsubsteps inside Predict.
2. Currently, IDGPredict is implemented using a custom Interface. We should factor out this code and convert it to the generic DPStep interface, so it can also be used outside of DDECal.
3. Refactor the nested DDECal model pipeline so that it can be defined separately for every direction instead of one to rule them all. The metadata must not change in any of the nested model steps of DDECal. This requires refactoring the process function of DDECal.cc.
4. Remove the model column from the buffer and create a ColumnReader that can read model columns that can be used in DDECal.

From a user perspective we propose the following additional parset definitions for DDECal:
- `ddecal.modelnextsteps`: Defines the default steps for the models that need to be applied to every direction before using the output for calibration. Each step can have its own parset entries.
- `ddecal.modelnextsteps.<direction>`: Allows steps to be defined for a specific direction. Overrides the default steps for this direction.
- `ddecal.modeldatacolumns`: (optional) Names of the columns that are to be used in DDECal as models.

An example of a parset for DDECal using different models for different directions is as follows:
```
ddecal.images = [image1.fits, image2.fits]  
ddecal.regions = ds9.reg  
ddecal.sourcedb = ATeam.sourcedb  
ddecal.modeldatacolumns = [sol]
ddecal.directions = [CasA, CygA] # subselection of patches from sourcedb

ddecal.modelnextsteps = [applycal3]  
ddecal.modelnextsteps.CasA = [applycal4, uvwflag3, applybeam3]  
ddecal.modelnextsteps.CygA = [applycal5, applybeam4]  

applycal5.h5parm = known_corruptions.h5 # Will automagically get right direction
```
