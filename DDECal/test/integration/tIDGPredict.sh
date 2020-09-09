#!/bin/bash

set -e

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DDECal/test/integration.
  exit 1;
fi
source $INIT
gunzip -c -f $srcdir/foursources.fits.gz > foursources.fits
cp $srcdir/foursources.reg .

echo "IDG calibration using four sources"
cmd="$dpppexe checkparset=1 msin=tDDECal.MS msout=.\
  steps=[ddecal] ddecal.useidg=True ddecal.idg.regions=foursources.reg\
  ddecal.idg.images=[foursources.fits]"
echo $cmd
$cmd >& /dev/null
