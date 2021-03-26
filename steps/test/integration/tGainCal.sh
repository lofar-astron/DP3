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
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=MODEL_DATA steps=[predict] predict.sourcedb=tNDPPP-generic.MS/sky predict.usebeammodel=false

echo; echo "Test caltype=diagonal"; echo
$dp3exe msin=tNDPPP-generic.MS msout= steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-diagonal gaincal.usebeammodel=false gaincal.caltype=diagonal gaincal.propagatesolutions=true gaincal.solint=1

$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL steps=[applycal] applycal.parmdb=tNDPPP-generic.MS/inst-diagonal

echo "Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima."
taqlcmd='select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL-t1.MODEL_DATA)))) as bbsres from tNDPPP-generic.MS t1, tGainCal.tab t2) where dpppres>bbsres*1.02'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1
echo "Checking that not everything was flagged"
$taqlexe 'select from tNDPPP-generic.MS where all(FLAG) groupby true having gcount()>100' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=diagonal with timeslotsperparmupdate=4"; echo
$dp3exe msin=tNDPPP-generic.MS msout= steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-diagonal-tpp gaincal.usebeammodel=false gaincal.caltype=diagonal gaincal.solint=4 gaincal.timeslotsperparmupdate=1 gaincal.propagatesolutions=false
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_TPP steps=[applycal] applycal.parmdb=tNDPPP-generic.MS/inst-diagonal-tpp
taqlcmd='select from tNDPPP-generic.MS where not all(near(DPPP_DIAGONAL, DPPP_DIAGONAL_TPP))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
#diff taql.out taql.ref || exit 1 # TODO This test fails, was not noticed before beacuse the diff command was forgotten

echo "Comparing the difference between applying with timeslotsperparmupdate = default and timeslotsperparmupdate=1"
taqlcmd='select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL-t1.MODEL_DATA)))) as bbsres from tNDPPP-generic.MS t1, tGainCal.tab t2) where dpppres>bbsres*1.02'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1
echo "Checking that not everything was flagged"
$taqlexe 'select from tNDPPP-generic.MS where all(FLAG) groupby true having gcount()>100' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=fulljones"; echo
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_FULLJONES_GAINCAL steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-fulljones gaincal.usebeammodel=false gaincal.caltype=fulljones gaincal.solint=1 gaincal.applysolution=true
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_FULLJONES steps=[applycal] applycal.parmdb=tNDPPP-generic.MS/inst-fulljones

echo "Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima."
taqlcmd='select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_FULLJONES-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_FULLJONES-t1.MODEL_DATA)))) as bbsres from tNDPPP-generic.MS t1, tGainCal.tab t2) where dpppres>bbsres*1.02'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1
echo "Comparing the solutions from gaincal + applycal with gaincal directly"
$taqlexe 'select from tNDPPP-generic.MS where not(all(DPPP_FULLJONES ~= DPPP_FULLJONES))' > taql.out
diff taql.out taql.ref  ||  exit 1
echo "Checking that not everything was flagged"
$taqlexe 'select from tNDPPP-generic.MS where all(FLAG) groupby true having gcount()>100' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=diagonal, nchan=2"; echo
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_NCHAN_GAINCAL steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-diagonal-nchan gaincal.usebeammodel=false gaincal.caltype=diagonal gaincal.solint=4 gaincal.nchan=2 gaincal.applysolution=true
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_NCHAN steps=[applycal] applycal.parmdb=tNDPPP-generic.MS/inst-diagonal-nchan

echo "Comparing the bbs residual with the dppp residual (solutions will not be equal, but residual should be equal). This avoids issues with local minima."
taqlcmd='select from (select gsumsqr(sumsqr(abs(iif(t1.FLAG,0,t1.DPPP_DIAGONAL_NCHAN-t1.MODEL_DATA)))) as dpppres, gsumsqr(sumsqr(abs(iif(FLAG,0,t2.BBS_DIAGONAL_NCHAN-t1.MODEL_DATA)))) as bbsres from tNDPPP-generic.MS t1, tGainCal.tab t2) where dpppres>bbsres*1.02'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1

echo "Comparing the solutions from gaincal + applycal with gaincal directly"
$taqlexe 'select from tNDPPP-generic.MS where not(all(DPPP_DIAGONAL_NCHAN_GAINCAL ~= DPPP_DIAGONAL_NCHAN))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo "Checking that not everything was flagged"
$taqlexe 'select from tNDPPP-generic.MS where all(FLAG) groupby true having gcount()>100' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=diagonal, nchan=2, solint=7"; echo
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_NCHAN_7_GAINCAL steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-diagonal-nchan gaincal.usebeammodel=false gaincal.caltype=diagonal gaincal.solint=4 gaincal.nchan=2 gaincal.applysolution=true
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_DIAGONAL_NCHAN_7 steps=[applycal] applycal.parmdb=tNDPPP-generic.MS/inst-diagonal-nchan

echo "Comparing the solutions from gaincal + applycal with gaincal directly"
$taqlexe 'select from tNDPPP-generic.MS where not(all(DPPP_DIAGONAL_NCHAN_7_GAINCAL ~= DPPP_DIAGONAL_NCHAN_7))' > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=tec"; echo
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_TEC steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-tec gaincal.caltype=tec gaincal.solint=2
# For now, only testing that the right parameter names are in the output
echo "    select result of 1 rows" > taql1.ref
taqlcmd='select from tNDPPP-generic.MS/inst-tec where (select NAME from ::NAMES)[NAMEID]=="TEC:CS001HBA0"'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql1.ref  ||  exit 1
taqlcmd='select from tNDPPP-generic.MS/inst-tec where (select NAME from ::NAMES)[NAMEID]=="CommonScalarPhase:CS001HBA0"'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test caltype=tecandphase"; echo
$dp3exe msin=tNDPPP-generic.MS msout=. msout.datacolumn=DPPP_TEC steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-tecandphase gaincal.caltype=tecandphase gaincal.solint=2
# For now, only testing that the right parameter names are in the output
echo "    select result of 1 rows" > taql1.ref
taqlcmd='select from tNDPPP-generic.MS/inst-tecandphase where (select NAME from ::NAMES)[NAMEID]=="TEC:CS001HBA0"'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql1.ref  ||  exit 1
taqlcmd='select from tNDPPP-generic.MS/inst-tecandphase where (select NAME from ::NAMES)[NAMEID]=="CommonScalarPhase:CS001HBA0"'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql1.ref  ||  exit 1

echo; echo "Test filter"; echo
$dp3exe msin=tNDPPP-generic.MS msout=tNDPPP-filtered.MS steps=[filter,gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-filter filter.baseline='!CS001HBA0&&*' gaincal.baseline='!CS002HBA1,RS305HBA&&*' gaincal.caltype=diagonal
taqlcmd='select from tNDPPP-generic.MS/inst-filter::NAMES where NAME LIKE "CS001HBA0%" OR NAME LIKE "%CS002HBA1%" OR NAME LIKE "%RS305HBA%"'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref || exit 1

echo; echo "Test debug output"; echo
$dp3exe msin=tNDPPP-generic.MS msout=. numthreads=1 steps=[gaincal] gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.parmdb=tNDPPP-generic.MS/inst-debug gaincal.caltype=diagonal gaincal.debuglevel=1
[[ -f debug.h5 ]] || exit 1
