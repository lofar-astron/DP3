#!/bin/bash

set -e # Stop on any error

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DPPP/test.
  exit 1;
fi
source $INIT

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

echo; echo "Creating parmdb with defvalues 3"
parmdbm <<EOL
open table="tApplyCal.parmdb"
adddef Gain:0:0:Real values=3.
adddef Gain:1:1:Real values=3.
EOL

echo; echo "Testing without updateweights"
cmd='NDPPP msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 msout.weightcolumn=WEIGHTS_NEW steps=[applycal] applycal.parmdb=tApplyCal.parmdb showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=9*DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1
$taqlexe 'select from tNDPPP-generic.MS where not(all(WEIGHTS_NEW~=WEIGHT_SPECTRUM))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing with updateweights"
cmd='NDPPP msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 msout.weightcolumn=WEIGHTS_NEW steps=[applycal] applycal.parmdb=tApplyCal.parmdb showcounts=false applycal.updateweights=true'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(WEIGHTS_NEW~=81*WEIGHT_SPECTRUM))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing CommonScalarPhase"
rm -rf tApplyCal.parmdb
parmdbm <<EOL
open table="tApplyCal.parmdb"
adddef CommonScalarPhase values=0
EOL
cmd='NDPPP msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=commonscalarphase showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing ScalarPhase"
rm -rf tApplyCal.parmdb
parmdbm <<EOL
open table="tApplyCal.parmdb"
adddef ScalarPhase values=0
EOL
cmd='NDPPP msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=commonscalarphase showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing ScalarPhase values"
rm -rf tApplyCal.parmdb
parmdbm <<EOL
open table="tApplyCal.parmdb"
add ScalarPhase:CS001HBA0 values=[0.,0.,0.,0.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarPhase:CS002HBA0 values=[0.,0.,0.,0.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarPhase:CS002HBA1 values=[0.,0.,0.,0.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarPhase:CS004HBA1 values=[0.,0.,0.,0.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarPhase:RS106HBA values=[0.,0.,0.,0.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarPhase:RS208HBA values=[0.,0.,0.,0.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarPhase:RS305HBA values=[0.,0.,0.,0.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarPhase:RS307HBA values=[0.,0.,0.,0.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
EOL
cmd='NDPPP msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=commonscalarphase showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing ScalarAmplitude values"
rm -rf tApplyCal.parmdb
parmdbm <<EOL
open table="tApplyCal.parmdb"
add ScalarAmplitude:CS001HBA0 values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarAmplitude:CS002HBA0 values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarAmplitude:CS002HBA1 values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarAmplitude:CS004HBA1 values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarAmplitude:RS106HBA values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarAmplitude:RS208HBA values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarAmplitude:RS305HBA values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add ScalarAmplitude:RS307HBA values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
EOL
cmd='NDPPP msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=scalaramplitude showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=9*DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Testing RotationAngle:*:phase_center values"
rm -rf tApplyCal.parmdb
parmdbm <<EOL
open table="tApplyCal.parmdb"
add RotationMeasure:CS001HBA0:phase_center values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:CS002HBA0:phase_center values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:CS002HBA1:phase_center values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:CS004HBA1:phase_center values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:RS106HBA:phase_center values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:RS208HBA:phase_center values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:RS305HBA:phase_center values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:RS307HBA:phase_center values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
EOL
cmd='NDPPP msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=rotationmeasure showcounts=false'
echo $cmd
eval $cmd

echo; echo "Testing RotationAngle:* values"
rm -rf tApplyCal.parmdb
parmdbm <<EOL
open table="tApplyCal.parmdb"
add RotationMeasure:CS001HBA0 values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:CS002HBA0 values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:CS002HBA1 values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:CS004HBA1 values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:RS106HBA values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:RS208HBA values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:RS305HBA values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
add RotationMeasure:RS307HBA values=[3.,3.,3.,3.], domain=[10e6, 200e6, 4472025735, 4972025795], shape=[2,2], shape=[2,2]
EOL
cmd='NDPPP msin=tNDPPP-generic.MS msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.parmdb=tApplyCal.parmdb applycal.correction=rotationmeasure showcounts=false'
echo $cmd
eval $cmd
