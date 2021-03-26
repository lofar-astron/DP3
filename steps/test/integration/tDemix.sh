#!/bin/bash

# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

set -e # Stop on any error

INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DPPP/test.
  exit 1;
fi
source $INIT

# Unpack the MS and other files and do the DPPP run.
tar zxf ${srcdir}/tDemix.in_MS.tgz
cd tDemix_tmp

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

echo; echo "### Test without target."; echo
rm -rf instrument
$dp3exe tDemix.parset demix.ignoretarget=true
# Compare some columns of the output MS with the reference output.
$taqlexe 'select from tDemix_out.MS t1, tDemix_ref1.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  not all(t1.LOFAR_FULL_RES_FLAG = t2.LOFAR_FULL_RES_FLAG)  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME' > taql.out
diff taql.out taql.ref  ||  exit 1
# Compare the instrument table.
##$taqlexe 'select from instrument t1, instrument_ref1 t2 where t1.NAMEID!=t1.NAMEID  ||  t1.STARTX!~=t2.STARTX  ||  t1.ENDX!~=t2.ENDX  ||  t1.STARTY!~=t2.STARTY  ||  t1.ENDY!~=t2.ENDY  ||  not all(near(t1.VALUES, t2.VALUES))' > taql.out
##diff taql.out taql.ref  ||  exit 1

echo; echo "### Test with target projected away."; echo
rm -rf instrument
$dp3exe tDemix.parset demix.ignoretarget=false
# Compare some columns of the output MS with the reference output.
$taqlexe 'select from tDemix_out.MS t1, tDemix_ref2.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  not all(t1.LOFAR_FULL_RES_FLAG = t2.LOFAR_FULL_RES_FLAG)  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME' > taql.out
diff taql.out taql.ref  ||  exit 1
# Compare the instrument table.
##$taqlexe 'select from instrument t1, instrument_ref2 t2 where t1.NAMEID!=t1.NAMEID  ||  t1.STARTX!~=t2.STARTX  ||  t1.ENDX!~=t2.ENDX  ||  t1.STARTY!~=t2.STARTY  ||  t1.ENDY!~=t2.ENDY  ||  not all(near(t1.VALUES, t2.VALUES))' > taql.out
##diff taql.out taql.ref  ||  exit 1

echo; echo "### Test with target."; echo
rm -rf instrument
$dp3exe tDemix.parset demix.target=CIZA.SP1A.FITS.pbcor_patch_s537 demix.freqstep=32 demix.timestep=5
# Compare some columns of the output MS with the reference output.
# $taqlexe 'select from tDemix_out.MS t1, tDemix_ref3.MS t2 where not all(near(t1.DATA,t2.DATA,1e-3) || (isnan(t1.DATA) && isnan(t2.DATA)))  ||  not all(t1.FLAG = t2.FLAG)  ||  not all(near(t1.WEIGHT_SPECTRUM, t2.WEIGHT_SPECTRUM))  ||  not all(t1.LOFAR_FULL_RES_FLAG = t2.LOFAR_FULL_RES_FLAG)  ||  t1.ANTENNA1 != t2.ANTENNA1  ||  t1.ANTENNA2 != t2.ANTENNA2  ||  t1.TIME !~= t2.TIME' > taql.out
# diff taql.out taql.ref  ||  exit 1
# Compare the instrument table.
##$taqlexe 'select from instrument t1, instrument_ref3 t2 where t1.NAMEID!=t1.NAMEID  ||  t1.STARTX!~=t2.STARTX  ||  t1.ENDX!~=t2.ENDX  ||  t1.STARTY!~=t2.STARTY  ||  t1.ENDY!~=t2.ENDY  ||  not all(near(t1.VALUES, t2.VALUES))' > taql.out
##diff taql.out taql.ref  ||  exit 1

exit 0
