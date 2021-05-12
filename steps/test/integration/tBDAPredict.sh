#!/bin/bash

# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DPPP/test.
  exit 1;
fi
source $INIT

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

cat > test.skymodel <<EOF
FORMAT = Name, Type, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation
point-0, POINT, 0.4362457236387493, 0.5287469737178224, 1.0, 0, 0, 0, , ,
EOF

$makesourcedbexe in=test.skymodel out=test.sourcedb append=false

echo; echo "Test BDAPredict"; echo

bda_args="bdaaverager.timebase=20000 bdaaverager.frequencybase=20000 bdaaverager.maxinterval=61"

cmd="$dp3exe msin=tNDPPP-generic.MS msout=bdapredict0.MS msout.overwrite=T steps=[predict,bdaaverager] $bda_args predict.sourcedb=test.sourcedb predict.usebeammodel=F"
echo $cmd
$cmd

cmd="$dp3exe msin=tNDPPP-generic.MS msout=bdapredict1.MS msout.overwrite=T steps=[bdaaverager,predict] $bda_args predict.sourcedb=test.sourcedb predict.usebeammodel=F"
echo $cmd
$cmd

# Compare the DATA columns of the output MSs.
# Because of flagging and weighting the difference can be as large
# the difference between visibilities at the edge and centre of the interval.
# For a scenario without flagging and with uniform weighting the tolerance can be lower.
taqlcmd='select ANTENNA1, ANTENNA2 from bdapredict0.MS t1, bdapredict1.MS t2 where not all(abs(t1.DATA-t2.DATA)<15e-2 || t1.FLAG || t1.WEIGHT_SPECTRUM==0)'
echo $taqlcmd
$taqlexe -noph $taqlcmd > taql.out
diff taql.out taql.ref  ||  exit 1

