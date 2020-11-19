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

tar zxf ${srcdir}/tGainCal.tab.tgz

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

echo "Creating MODEL_DATA so that residual can be computed"
cmd="$dpppexe checkparset=1 showprogress=false msin=tNDPPP-generic.MS msout=. msout.datacolumn=MODEL_DATA steps=[predict] predict.sourcedb=tNDPPP-generic.MS/sky predict.usebeammodel=false"
echo $cmd
$cmd

echo; echo "Test caltype=diagonal"; echo
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout= steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-diagonal.h5 gaincal.usebeammodel=false gaincal.caltype=diagonal gaincal.propagatesolutions=true gaincal.solint=1"
echo $cmd
$cmd

cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL steps=[applycal] applycal.parmdb=tNDPPP-generic.MS/inst-diagonal.h5 applycal.steps=[amplitude,phase] applycal.phase.correction=phase000 applycal.amplitude.correction=amplitude000 applycal.amplitude.correction=amplitude000"
echo $cmd
$cmd

echo "Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima."
$taqlexe 'select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL-t1.MODEL_DATA)))) as bbsres from tNDPPP-generic.MS t1, tGainCal.tab t2) where dpppres>bbsres*1.02' > taql.out
diff taql.out taql.ref  ||  exit 1
echo "Checking that not everything was flagged"
$taqlexe 'select from tNDPPP-generic.MS where all(FLAG) groupby true having gcount()>100' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=fulljones"; echo
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_FULLJONES_GAINCAL steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-fulljones.h5 gaincal.usebeammodel=false gaincal.caltype=fulljones gaincal.solint=1 gaincal.applysolution=true"
echo $cmd
$cmd

echo; echo "Test caltype=diagonal, nchan=2"; echo
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_NCHAN_GAINCAL steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-diagonal-nchan.h5 gaincal.usebeammodel=false gaincal.caltype=diagonal gaincal.solint=4 gaincal.nchan=2 gaincal.applysolution=true"
echo $cmd
$cmd
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_NCHAN steps=[applycal] applycal.parmdb=tNDPPP-generic.MS/inst-diagonal-nchan.h5 applycal.steps=[phase,amplitude] applycal.phase.correction=phase000 applycal.amplitude.correction=amplitude000"
echo $cmd
$cmd

echo "Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima."
$taqlexe 'select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL_NCHAN-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL_NCHAN-t1.MODEL_DATA)))) as bbsres from tNDPPP-generic.MS t1, tGainCal.tab t2) where dpppres>bbsres*1.02' > taql.out
diff taql.out taql.ref  ||  exit 1

echo "Comparing the solutions from gaincal + applycal with gaincal directly"
$taqlexe 'select from tNDPPP-generic.MS where not(all(DPPP_DIAGONAL_NCHAN_GAINCAL ~= DPPP_DIAGONAL_NCHAN))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo "Checking that not everything was flagged"
$taqlexe 'select from tNDPPP-generic.MS where all(FLAG) groupby true having gcount()>100' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=diagonal, nchan=2, solint=7"; echo
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_NCHAN_7_GAINCAL steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-diagonal-nchan7.h5 gaincal.usebeammodel=false gaincal.caltype=diagonal gaincal.solint=4 gaincal.nchan=2 gaincal.applysolution=true"
echo $cmd
$cmd
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_NCHAN_7 steps=[applycal] applycal.parmdb=tNDPPP-generic.MS/inst-diagonal-nchan7.h5 applycal.steps=[amplitude,phase] applycal.amplitude.correction=amplitude000 applycal.phase.correction=phase000"
echo $cmd
$cmd

echo "Comparing the solutions from gaincal + applycal with gaincal directly"
$taqlexe 'select from tNDPPP-generic.MS where not(all(DPPP_DIAGONAL_NCHAN_7_GAINCAL ~= DPPP_DIAGONAL_NCHAN_7))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=tec"; echo
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_TEC steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-tec.h5 gaincal.caltype=tec gaincal.solint=2"
echo $cmd
$cmd

echo; echo "Test caltype=tecandphase"; echo
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_TEC steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-tecandphase.h5 gaincal.caltype=tecandphase gaincal.solint=2"
echo $cmd
$cmd

echo; echo "Test filter"; echo
cmd="$dpppexe checkparset=1 msin=tNDPPP-generic.MS msout=tNDPPP-filtered.MS steps=[filter,gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-filter.h5 filter.baseline='!CS001HBA0&&*' gaincal.baseline='!CS002HBA1,RS305HBA&&*' gaincal.caltype=diagonal"
echo $cmd
$cmd
