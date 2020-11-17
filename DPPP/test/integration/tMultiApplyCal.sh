#!/bin/bash

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DPPP/test.
  exit 1;
fi
source $INIT

set -e # Stop on any error

tar zxf ${srcdir}/tApplyCal2.parmdb.tgz

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

cmd='$dpppexe msin=tNDPPP-generic.MS checkparset=1 msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.steps="[gain,csp]" applycal.gain.parmdb=tApplyCal.parmdb applycal.gain.correction=gain applycal.csp.parmdb=tApplyCal.parmdb applycal.csp.correction=commonscalarphase showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=9*DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1

cmd='$dpppexe msin=tNDPPP-generic.MS checkparset=1 msout=. msout.datacolumn=DATA3 steps=[applycal] applycal.steps="[gain,csp]" applycal.parmdb=tApplyCal.parmdb applycal.gain.correction=gain applycal.csp.correction=commonscalarphase showcounts=false'
echo $cmd
eval $cmd
$taqlexe 'select from tNDPPP-generic.MS where not(all(DATA~=9*DATA3))' > taql.out
diff taql.out taql.ref  ||  exit 1

