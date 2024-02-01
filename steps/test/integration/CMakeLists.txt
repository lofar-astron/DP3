# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

add_python_tests(
  tApplyBeam
  tApplyCal2
  tBdaExpander
  tBdaPredict
  tClipper
  tColumnReader
  tDemix
  tGainCal
  tGainCalH5Parm
  tIDGImager
  tMsIn
  tMsOut
  tMultiApplyCal
  tPhaseshiftPredict
  tPredict
  tReadOnly
  tNullStokes
  tSplit)

if(LIBDIRAC_FOUND)
  add_python_tests(tSagecalPredict)
endif()
