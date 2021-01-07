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

chmod a-w -R tNDPPP-generic.MS
chmod u+w -R tNDPPP-generic.MS/sky
cmd="$dpppexe msin=tNDPPP-generic.MS steps=[gaincal] gaincal.parmdb=tNDPPP-generic.MS/inst-diagonal gaincal.sourcedb=tNDPPP-generic.MS/sky gaincal.caltype=diagonal gaincal.parmdb=gaincal.h5 msout="
echo $cmd
$cmd || exit 1

# ddecal.sourcedb=tNDPPP-generic.MS/sky ddecal.mode=diagonal ddecal.h5parm=dde.h5
