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

text_boldred="$(tput bold)$(tput setaf 1)"
text_boldgreen="$(tput bold)$(tput setaf 2)"
text_boldcyan="$(tput bold)$(tput setaf 6)"
text_normal=$(tput sgr0)

compare_results() {
  # Ignore baselines with antennna 6 for now, since not all predictors
  # generate visibilities for those baselines.
  # TODO (AST-223): Investigate what the expected behavior is.
  $taqlexe "select from tDDECal.MS where ANTENNA2 != 6 and not all(near($1_DATA,IDG_DATA,1e-3))" > taql.out
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
  steps=[ddecal] ddecal.useidg=True ddecal.idg.regions=foursources.reg\
  ddecal.idg.images=[foursources-model.fits]\
  ddecal.onlypredict=True msout.datacolumn=IDG_DATA"
echo $cmd
$cmd
compare_results foursources "(multiple)"

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
          steps=[ddecal] ddecal.useidg=True ddecal.idg.regions=$source-$offset.reg\
          ddecal.idg.images=[$source-model.fits]\
          ddecal.onlypredict=True msout.datacolumn=IDG_DATA"
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

echo Tests passed: $passed skipped: $skipped failed: $failed
exit $failed
