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

tar zxf ${srcdir}/tApplyBeam.tab.tgz

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

cmd="$dp3exe msin=tNDPPP-generic.MS steps=[split] split.steps=[applybeam,out] split.replaceparms=[out.name,applybeam.usechannelfreq] out.name=[splitout1.ms,splitout2.ms] applybeam.usechannelfreq=[false,true] applybeam.invert=true"
echo $cmd
$cmd

# Compare the DATA column of the output MS for usechannelfreq=false with the BBS reference output.
taqlcmd='select from splitout1.ms t1, tApplyBeam.tab t2 where not all(near(t1.DATA,t2.DATA_noucf,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA_noucf)))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1

# Compare the DATA column of the output MS for usechannelfreq=true with the BBS reference output.
taqlcmd='select from splitout2.ms t1, tApplyBeam.tab t2 where not all(near(t1.DATA,t2.DATA_ucf,8e-5) || (isnan(t1.DATA) && isnan(t2.DATA_ucf)))'
echo $taqlcmd
$taqlexe "$taqlcmd" > taql.out
diff taql.out taql.ref  ||  exit 1
