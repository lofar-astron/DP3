#!/bin/bash

set -e

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DDECal/test/integration.
  exit 1;
fi
source $INIT
gunzip -c -f $srcdir/foursources.fits.gz > foursources-model.fits
cp $srcdir/foursources.reg .

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

$taqlexe 'update tDDECal.MS set WEIGHT_SPECTRUM=1, FLAG=False'

#Use wsclean for generating a reference prediction, in the MODEL_DATA column.
wsclean -predict -name foursources tDDECal.MS

echo "Predict four point sources using IDG"
cmd="$dpppexe checkparset=1 msin=tDDECal.MS msout=.\
  steps=[ddecal] ddecal.useidg=True ddecal.idg.regions=foursources.reg\
  ddecal.idg.images=[foursources-model.fits]\
  ddecal.onlypredict=True msout.datacolumn=IDG_DATA"
echo $cmd
$cmd

cmd="$taqlexe 'select from tDDECal.MS where not all(near(MODEL_DATA,IDG_DATA,1e-3))' > taql.out"
echo $cmd
eval $cmd
diff -q taql.out taql.ref || exit 1
