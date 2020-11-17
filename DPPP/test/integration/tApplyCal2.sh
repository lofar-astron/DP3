#!/bin/bash

set -e # Stop on any error

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DPPP/test.
  exit 1;
fi
source $INIT

tar zxf ${srcdir}/tApplyCal2.parmdb.tgz

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

echo; echo "Testing without updateweights"
cmd='$dpppexe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 msout.weightcolumn=WEIGHTS_NEW steps=[applycal] applycal.parmdb=tApplyCal.parmdb showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=9*DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1
$taqlexe 'select from tNDPPP-generic.MS where not(all(WEIGHTS_NEW~=WEIGHT_SPECTRUM))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing with updateweights"
cmd='$dpppexe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 msout.weightcolumn=WEIGHTS_NEW steps=[applycal] applycal.parmdb=tApplyCal.parmdb showcounts=false applycal.updateweights=true'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(WEIGHTS_NEW~=81*WEIGHT_SPECTRUM))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing CommonScalarPhase"
cmd='$dpppexe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=commonscalarphase showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing ScalarAmplitude values"
cmd='$dpppexe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=scalaramplitude showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=9*DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing RotationAngle:*:phase_center values"
cmd='$dpppexe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=rotationmeasure showcounts=false'
echo $cmd
eval $cmd
