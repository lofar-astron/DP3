# Full resolution flags

## Introduction
After averaging a dataset, one can have the necessity to go back and check the
flag's value at full resolution. That's the reason why we have the fullResFlags
field in the [RegularBuffer](@ref dp3::base::DPBuffer).

## When are they updated?
The full resolution flags are only updated with the flag's values when the 
Averager is called. Other steps which read/use them will not update the full
resolution flags.

NB. This means that the flags and the full resolution flags are not always
synchronized with each other, but the full resolution flags will be synchronized
with the flag values in the Averager step.

## Example 1

Let's assume the following settings: 
msin=in.ms 
msout=out.ms 
steps=[preflagger, averager]

1. Since the averager's and writer's required fields contain fullResFlags, the 
MSReader will read/generate them (copying from the flags if not available in 
the MS).

2. Let's assume that the preflagger will flag the 4th and 5th data indexes 
(fullResFlags are not updated).

3. The averager averages the five times to 1, resulting in a discrepancy between 
the size of flags and fullresflags. Before compressing the flags, they are 
OR-ed with the fullResFlags.

NB. The averager will flag the averaged data if more than a certain percentage 
of the full resolution data was flagged. In this case it's not flagging the 
averaged data.

4. The full resolution flags are saved by default in the MS (that can be 
disabled via the msout.writefullresflags parameter).

<table>
<tr><th>Step <th>Flags <th>FullResFlags
<tr><td rowspan="1">MSReader <td> 0 0 0 0 0 <td> 0 0 0 0 0
<tr><td rowspan="1">PreFlagger <td> 0 0 0 1 1 <td> 0 0 0 0 0
<tr><td rowspan="1">Averager <td> 0 <td> 0 0 0 1 1
<tr><td rowspan="1">MSWriter <td> 0 <td> 0 0 0 1 1
</table>

## Example 2

Let's assume the following settings: 
msin=in.ms 
msout=out.ms 
steps=[averager, aoflagger]

1. Since the writer's and averager's required fields contain fullResFlags, the 
MSReader will read/generate them (copying from the flags if not available in 
the MS).

2. The averager averages the five times to one, resulting in a discrepancy between 
the size of flags and fullresflags.

3. Let's assume that the aoflagger will flag the data (fullResFlags are not 
updated):

4. The full resolution flags are saved by default in the MS (that can be 
disabled via the msout.writefullresflags parameter).

NB. The flags and the full resolution flags may not be consistent with each 
other due to the design, and that's the intended behavior.

<table>
<tr><th>Step <th>Flags <th>FullResFlags
<tr><td rowspan="1">MSReader <td> 0 0 0 0 0 <td> 0 0 0 0 0
<tr><td rowspan="1">Averager <td> 0 <td> 0 0 0 0 0
<tr><td rowspan="1">AOFlagger <td> 1 <td> 0 0 0 0 0
<tr><td rowspan="1">MSWriter <td> 1 <td> 0 0 0 0 0
</table>


## Example 3

Let's assume the following settings: 
msin=in.ms 
msout=out.ms 
steps=[aoflagger] 
msout.writefullresflag=false

1. Since the writer's and aoflagger's required fields do not contain 
fullResFlags, the MSReader will not generate them.

2. Let's assume that the aoflagger will flag the data.

3. The full resolution flags are not saved in the MS, since the 
msout.writefullresflag parameter is set to false.

<table>
<tr><th>Step <th>Flags <th>FullResFlags
<tr><td rowspan="1">MSReader <td> 0 0 0 0 0 <td> n.a.
<tr><td rowspan="1">AOFlagger <td> 1 1 1 1 1 <td> n.a.
<tr><td rowspan="1">MSWriter <td> 1 1 1 1 1 <td> n.a.
</table>

## Example 4

Let's assume the following settings: 
msin=in.ms 
msout=out.ms 
steps=[aoflagger]

Behaviour before AST-1047

1. The reader does not read fullResFlags as they are not needed by the reader 
step.

2. Let's assume that the aoflagger will flag the data.

3. The MSWriter will read/generate the full resolution flags (copying from the 
flags if not available in the MS) and save them to the MS.

<table>
<tr><th>Step <th>Flags <th>FullResFlags
<tr><td rowspan="1">MSReader <td> 0 0 0 0 0 <td> n.a.
<tr><td rowspan="1">AOFlagger <td> 1 1 1 1 1 <td> n.a.
<tr><td rowspan="1">MSWriter <td> 1 1 1 1 1 <td> 1 1 1 1 1
</table>

Behaviour after AST-1047

1. Since the writer's required fields contain fullResFlags, the MSReader will 
read/generate them (copying from the flags if not available in 
the MS).

2. Let's assume that the aoflagger will flag the data.

3. The full resolution flags are saved in the MS.

<table>
<tr><th>Step <th>Flags <th>FullResFlags
<tr><td rowspan="1">MSReader <td> 0 0 0 0 0 <td> 0 0 0 0 0
<tr><td rowspan="1">AOFlagger <td> 1 1 1 1 1 <td> 0 0 0 0 0
<tr><td rowspan="1">MSWriter <td> 1 1 1 1 1 <td> 0 0 0 0 0
</table>