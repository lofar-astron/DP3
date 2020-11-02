#!/bin/bash

set -e

# Locate the executables and srcdir (script created by cmake's configure_file).
INIT=testInit.sh
if [ ! -f $INIT ]; then
  echo $INIT not found. Please run this script from build/DDECal/test/integration.
  exit 1;
fi
source $INIT

tar xfj $srcdir/resources/idg-fits-sources.tbz2
cp $srcdir/resources/*.reg .

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref

passed=0
failed=0
skipped=0

# Prevent and ignore tput errors an CI.
if [ -n "$CI" ]; then export TERM=dumb; set +e; fi
text_boldred="$(tput bold)$(tput setaf 1)"
text_boldgreen="$(tput bold)$(tput setaf 2)"
text_boldcyan="$(tput bold)$(tput setaf 6)"
text_normal=$(tput sgr0)
if [ -n "$CI" ]; then set -e; fi

compare_results() {
  # Ignore baselines with antennna 6 for now, since not all predictors
  # generate visibilities for those baselines.
  # TODO (AST-223): Investigate what the expected behavior is.
  $taqlexe "select from tDDECal.MS where ANTENNA2 != 6 and not all(near($1_DATA,MODEL_DATA,1e-3))" > taql.out
  echo -n "Predict source: $1 offset: $2 result: "
  if diff -q taql.out taql.ref; then
    echo "${text_boldgreen}Test passed.${text_normal}"
    passed=$((passed + 1))
  else
    echo "${text_boldred}Test failed.${text_normal}"
    failed=$((failed + 1))
  fi
}

# Test an input with four sources.

# Since wsclean on CI does not support IDG, tDDECal.MS has a foursources_DATA
# column with the reference output/visibilities for foursources-model.fits.
# These commands generated the column:
#wsclean -use-idg -predict -name resources/foursources tDDECal.MS
#taql "alter table tDDECal.MS rename column MODEL_DATA to foursources_DATA"

echo "Predict four sources using IDG"
cmd="$dpppexe checkparset=1 msin=tDDECal.MS msout=.\
  steps=[ddecal] ddecal.idg.regions=foursources.reg\
  ddecal.idg.images=[foursources-model.fits]\
  ddecal.onlypredict=True msout.datacolumn=MODEL_DATA"
echo $cmd
$cmd
compare_results foursources "(multiple)"

# Test if IDGPredict step will have the same results as DDECal
cmd="$dpppexe checkparset=1 msin=tDDECal.MS msout=idgout.MS\
  steps=[idgpredict] idgpredict.regions=foursources.reg\
  idgpredict.images=[foursources-model.fits]"
echo $cmd
$cmd

# Compare the MODEL_DATA column of the output MS with the original data minus the BBS reference output.
taqlcmd='select from idgout.MS t1, tDDECal.MS t2 where not all(near(t1.DATA,t2.MODEL_DATA,5e-2) || (isnan(t1.DATA) && isnan(t2.MODEL_DATA)))'
echo $taqlcmd
$taqlexe $taqlcmd > taql.out
diff -q taql.out taql.ref  ||  exit 1

# Test multiple data sources for DDECal
echo "Create model data column with 3 sources"
cmd="$dpppexe checkparset=1 msin=tDDECal.MS msout=.\
  steps=[ddecal] ddecal.idg.regions=threesources.reg\
  ddecal.idg.images=[foursources-model.fits]\
  ddecal.onlypredict=True msout.datacolumn=MODEL_DATA"
echo $cmd
$cmd >& /dev/null

echo "Run DDECal with 3 directions in the MODEL_DATA column and 1 direction using IDG"
cmd="$dpppexe checkparset=1 msin=tDDECal.MS msout=.\
  steps=[ddecal] ddecal.idg.regions=onesource.reg\
  ddecal.idg.images=[foursources-model.fits]\
  ddecal.onlypredict=True\
  ddecal.usemodelcolumn=true\
  msout.datacolumn=MULTIPLE_SOURCES"
echo $cmd
$cmd >& /dev/null

# Results of a the DDECal above (3 directions modeldata and 1 direction IDG) should be equal to a run with foursources (4 direction IDG).
taqlcmd='select from idgout.MS t1, tDDECal.MS t2 where not all(near(t1.DATA,t2.MULTIPLE_SOURCES,5e-2) || (isnan(t1.DATA) && isnan(t2.MULTIPLE_SOURCES)))'
echo $taqlcmd
$taqlexe $taqlcmd > taql.out
diff taql.out taql.ref || exit 1

# Test inputs that contain a single source.
# Since these tests take quite some time, they only run locally, and only
# if the foursources test fails or if it is commented out.
# CI runs don't work since wsclean does not support IDG on CI.
if [ $passed != 1 -a -z "$CI" ]; then
  SOURCES="center ra dec radec"
  OFFSETS="center dl dm dldm"
  for source in $SOURCES; do
    #Generate reference predictions in the ${source}_DATA column.
    wsclean -use-idg -predict -name $source tDDECal.MS
    taql "alter table tDDECal.MS rename column MODEL_DATA to ${source}_DATA"

    for offset in $OFFSETS; do
      echo "Predict source: $source offset: $offset using IDG"
      if grep -q "^polygon" $source-$offset.reg; then
        cmd="$dpppexe checkparset=1 msin=tDDECal.MS msout=.\
          steps=[ddecal] ddecal.idg.regions=$source-$offset.reg\
          ddecal.idg.images=[$source-model.fits]\
          ddecal.onlypredict=True msout.datacolumn=MODEL_DATA"
        echo $cmd
        $cmd
        compare_results $source $offset
      else
        echo "${text_boldcyan}Test skipped: No polygon(s) defined!${text_normal}"
        skipped=$((skipped + 1))
      fi
    done
  done
fi

echo Test polynomial frequency term corrections...
cmd="$dpppexe checkparset=1 msin=tDDECal.MS msout=.\
  steps=[ddecal] ddecal.idg.regions=center-center.reg\
  ddecal.idg.images=[term0-model.fits,term1-model.fits,term2-model.fits]\
  ddecal.onlypredict=True msout.datacolumn=TERMS_DATA"
echo $cmd
$cmd

#Since the test involves the center pixel only, taql can check the values.
ch_count=`taql -nop "select count(CHAN_FREQ) from tDDECal.MS::SPECTRAL_WINDOW"`
term_failed=0
for ch in `seq 0 $((ch_count - 1))`; do
  # The factors 10, 20000 and 30000 match those in tIDGPredict_ref.py.
  taql "select from tDDECal.MS
        where not(TERMS_DATA[$ch,0]=0 or near(TERMS_DATA[$ch,0],
        (select 10+20000*(CHAN_FREQ[$ch]/CHAN_FREQ[0]-1)
                  +30000*(CHAN_FREQ[$ch]/CHAN_FREQ[0]-1)**2
         from ::SPECTRAL_WINDOW)[0], 1e-3))" > taql.out
  if ! diff -q taql.out taql.ref; then
    term_failed=$((term_failed + 1))
  fi
done

echo -n "Frequency term test result: "
if [ $term_failed -eq 0 ]; then
  echo "${text_boldgreen}Test passed.${text_normal}"
  passed=$((passed + 1))
else
  echo "${text_boldred}Test has $term_failed failures.${text_normal}"
  failed=$((failed + 1))
fi

echo Tests passed: $passed skipped: $skipped failed: $failed
exit $failed
