# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Since there are many DDECAL tests, defining a CTest for each subtest allows
# runnning those subtests in parallel, which speeds up testing.
add_python_tests(ONE_CTEST_PER_SUBTEST tDDECal)

# Running all subtests together takes at most 11 seconds for these tests
# on DAS-6. A single CTest, running all subtests, suffices here.
add_python_tests(tBdaDdeCal tFaradayConstraint tIDGPredict
                 tWGridderPredictDdeCal)
