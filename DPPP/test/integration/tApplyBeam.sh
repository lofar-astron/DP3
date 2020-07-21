#!/bin/bash

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DPPP/test.
  exit 1;
fi
source $INIT

tar zxf ${srcdir}/tApplyBeam.tab.tgz

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

echo; echo "Test with invert=true and usechannelfreq=false"; echo
cmd="$dpppexe msin=tNDPPP-generic.MS msout=outinv.ms steps=[applybeam] applybeam.usechannelfreq=false applybeam.invert=true"
echo $cmd
$cmd
# Compare the DATA column of the output MS with the BBS reference output.
taqlcmd='select from outinv.ms t1, tApplyBeam.tab t2 where not all(near(t1.DATA,t2.DATA_noucf,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA_noucf)))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "### Test with invert=false on the output of the previous step"; echo
cmd="$dpppexe msin=outinv.ms msout=out.ms steps=[applybeam] applybeam.usechannelfreq=false applybeam.invert=false"
echo $cmd
$cmd
# Compare the DATA column of the output MS with the original MS.
taqlcmd='select from out.ms t1, tNDPPP-generic.MS t2 where not all(near(t1.DATA,t2.DATA,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA)))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test with invert=true and usechannelfreq=true"; echo
cmd="$dpppexe msin=tNDPPP-generic.MS msout=outinv.ms msout.overwrite=true steps=[applybeam] applybeam.usechannelfreq=true applybeam.invert=true"
echo $cmd
$cmd
# Compare the DATA column of the output MS with the BBS reference output.
taqlcmd='select from outinv.ms t1, tApplyBeam.tab t2 where not all(near(t1.DATA,t2.DATA_ucf,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA_ucf)))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test with invert=false on the output of the previous step"; echo
cmd="$dpppexe msin=outinv.ms msout=out.ms msout.overwrite=true steps=[applybeam] applybeam.usechannelfreq=true applybeam.invert=false"
echo $cmd
$cmd
# Compare the DATA column of the output MS with the original MS.
taqlcmd='select from out.ms t1, tNDPPP-generic.MS t2 where not all(near(t1.DATA,t2.DATA,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA)))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test with beammode=ARRAY_FACTOR"; echo
cmd="$dpppexe msin=tNDPPP-generic.MS msout=outinv.ms msout.overwrite=true steps=[applybeam] applybeam.usechannelfreq=true applybeam.invert=true applybeam.beammode=ARRAY_FACTOR"
echo $cmd
$cmd
# Compare the DATA column of the output MS with the BBS reference output.
taqlcmd='select from outinv.ms t1, tApplyBeam.tab t2 where not all(near(t1.DATA,t2.DATA_ARRAY_FACTOR,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA_ARRAY_FACTOR)))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1

echo; echo "Test with beammode=ELEMENT"; echo
echo "!!!!!!!!!!! DISABLED, SEE https://github.com/lofar-astron/DP3/issues/273 !!!!!!!!!!!!"
# cmd="$dpppexe msin=tNDPPP-generic.MS msout=outinv.ms msout.overwrite=true steps=[applybeam] applybeam.usechannelfreq=true applybeam.invert=true applybeam.beammode=ELEMENT"
# echo $cmd
# $cmd
# # Compare the DATA column of the output MS with the BBS reference output.
# taqlcmd='select from outinv.ms t1, tApplyBeam.tab t2 where not all(near(t1.DATA,t2.DATA_ELEMENT,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA_ELEMENT)))'
# echo $taqlcmd
# $taqlexe "$taqlcmd" > taql.out
# diff taql.out taql.ref  ||  exit 1

echo; echo "Test with updateweights=true"; echo
cmd="$dpppexe msin=tNDPPP-generic.MS msout=. steps=[applybeam] applybeam.updateweights=truue msout.weightcolumn=NEW_WEIGHT_SPECTRUM"
echo $cmd
$cmd
# Check that the weights have changed
taqlcmd='select from tNDPPP-generic.MS where all(near(WEIGHT_SPECTRUM, NEW_WEIGHT_SPECTRUM))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1
