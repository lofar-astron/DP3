#!/bin/bash

# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

set -e

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DPPP/test.
  exit 1;
fi
source $INIT

$taqlexe 'update tDDECal.MS set WEIGHT_SPECTRUM=1, FLAG=False'

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

echo "Create corrupted model"
cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. steps=[ddecal]\
  ddecal.sourcedb=tDDECal.MS/sky ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]\
  ddecal.h5parm=instrumentcorrupted.h5 ddecal.mode=complexgain"
echo $cmd
$cmd >& /dev/null

# Modify h5 file for multiple solution intervals
python3 << EOF
import h5py
import numpy as np
h5file = h5py.File("instrumentcorrupted.h5", "r+")
sol = h5file['sol000/amplitude000/val']
sol[:4,:,:,0,:] = np.sqrt(5)
sol[4:,:,:,0,:] = np.sqrt(5+2)
sol[:4,:,:,1,:] = np.sqrt(9)
sol[4:,:,:,1,:] = np.sqrt(9+2)
sol[:4,:,:,2,:] = np.sqrt(13)
sol[4:,:,:,2,:] = np.sqrt(13+2)
h5file.close()
EOF

echo "Predict corrupted visibilities"
cmd="$dp3exe checkparset=1 msin=tDDECal.MS/ msout=. msout.datacolumn=DATA steps=[h5parmpredict]\
  h5parmpredict.sourcedb=tDDECal.MS/sky h5parmpredict.applycal.parmdb=instrumentcorrupted.h5\
  h5parmpredict.applycal.correction=amplitude000"
echo $cmd
$cmd >& /dev/null

echo "Predict model data column"
cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. msout.datacolumn=MODEL_DATA\
  steps=[]"
echo $cmd
$cmd >& /dev/null

for caltype in complexgain scalarcomplexgain amplitudeonly scalaramplitude
do
  for solint in 0 1 2 4
  do
    [[ $solint = 0 ]] && tolerance="0.15" || tolerance="0.015"
    for nchan in 1 2 5
    do
          echo "Calibrate on the original sources, caltype=$caltype"
          cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. steps=[ddecal]\
             ddecal.sourcedb=tDDECal.MS/sky ddecal.solint=$solint ddecal.nchan=$nchan \
             ddecal.directions=[[center,dec_off],[ra_off],[radec_off]] \
             ddecal.h5parm=instrument.h5 ddecal.mode=$caltype"
          echo $cmd
          $cmd

          echo "Apply solutions with multiple predict steps, caltype=$caltype"
          cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. msout.datacolumn=SUBTRACTED_DATA\
            steps=[predict1,predict2,predict3]\
              predict1.sourcedb=tDDECal.MS/sky\
              predict1.applycal.parmdb=instrument.h5 predict1.sources=[center,dec_off]\
              predict1.operation=subtract predict1.applycal.correction=amplitude000 \
              predict2.sourcedb=tDDECal.MS/sky\
              predict2.applycal.parmdb=instrument.h5 predict2.sources=[radec_off]\
              predict2.operation=subtract predict2.applycal.correction=amplitude000 \
              predict3.sourcedb=tDDECal.MS/sky\
              predict3.applycal.parmdb=instrument.h5 predict3.sources=[ra_off]\
              predict3.operation=subtract predict3.applycal.correction=amplitude000"
          echo $cmd
          $cmd

          #h5sols.py instrument.h5

          echo "Check that residual is small, caltype=$caltype, nchan=$nchan, solint=$solint"
          cmd="$taqlexe 'select norm_residual/norm_data FROM (select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*DATA))) as norm_data, sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*SUBTRACTED_DATA))) as norm_residual from tDDECal.MS)'"
          echo $cmd
          eval $cmd

          cmd="$taqlexe 'select FROM (select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*DATA))) as norm_data, sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*SUBTRACTED_DATA))) as norm_residual from tDDECal.MS) where norm_residual/norm_data > $tolerance or isinf(norm_residual/norm_data) or isnan(norm_residual/norm_data)' > taql.out"
          echo $cmd
          eval $cmd

          diff taql.out taql.ref || exit 1
    done # Loop over nchan
  done # Loop over solint
done # Loop over caltype

echo "Apply solutions with h5parmpredict"
cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. msout.datacolumn=SUBTRACTED_DATA_H5PARM\
  steps=[predict]\
    predict.type=h5parmpredict\
    predict.sourcedb=tDDECal.MS/sky\
    predict.applycal.parmdb=instrument.h5\
    predict.operation=subtract predict.applycal.correction=amplitude000"
echo $cmd
$cmd

echo "Check that h5parmpredict creates the same output as multiple predict steps"
cmd="$taqlexe 'select abs(gsumsqr(SUBTRACTED_DATA-SUBTRACTED_DATA)) as diff from tDDECal.MS'"
echo $cmd
eval $cmd
cmd="$taqlexe 'select from (select abs(gsumsqr(SUBTRACTED_DATA-SUBTRACTED_DATA)) as diff from tDDECal.MS) where diff>1.e-6' > taql.out"
echo $cmd
eval $cmd
diff taql.out taql.ref || exit 1

echo "Check pre-apply"
cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. numthreads=4\
  steps=[ddecal]\
    ddecal.sourcedb=tDDECal.MS/sky\
    ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]\
    ddecal.applycal.parmdb=instrument.h5 ddecal.applycal.steps=applyampl\
    ddecal.applycal.applyampl.correction=amplitude000\
    ddecal.h5parm=instrument2.h5 ddecal.mode=scalarcomplexgain"
echo $cmd
$cmd

echo "Check tec"
cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. numthreads=4\
  steps=[ddecal]\
    ddecal.sourcedb=tDDECal.MS/sky\
    ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]\
    ddecal.h5parm=instrument-tec.h5 ddecal.mode=tec"
echo $cmd
$cmd

echo "Check tec and phase"
cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. numthreads=4\
  msin.baseline='!CS001HBA0'\
  steps=[ddecal]\
    ddecal.sourcedb=tDDECal.MS/sky\
    ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]\
    ddecal.h5parm=instrument-tecandphase.h5 ddecal.mode=tecandphase"
echo $cmd
$cmd


echo "Create MODEL_DATA"
cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. msout.datacolumn=MODEL_DATA\
  steps=[]"
echo $cmd
$cmd >& /dev/null

echo "Multiply MODEL_DATA by 42"
$taqlexe "update tDDECal.MS set MODEL_DATA=DATA*42"

cmd="$dp3exe checkparset=1 msin=tDDECal.MS msout=. msout.datacolumn=SUBTRACTED_DATA\
  steps=[ddecal] ddecal.modeldatacolumns=[MODEL_DATA] ddecal.modelnextsteps.MODEL_DATA=[scaledata]\
  ddecal.h5parm=instrument-modeldata ddecal.solint=2 ddecal.nchan=3\
  scaledata.stations='*' scaledata.scalesize=False scaledata.coeffs=1\
  ddecal.subtract=True"
echo $cmd
$cmd

$taqlexe 'select from (select abs(sumsqr(SUBTRACTED_DATA)/sumsqr(DATA)) as diff from tDDECal.MS) where diff>1.e-6' > taql.out
diff taql.out taql.ref || exit 1
