# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

add_python_tests(
  tApplyBeam
  tApplyCal2
  tBdaExpander
  tBdaPredict
  tDemix
  tDemixerNew
  tGainCal
  tGainCalH5Parm
  tMsIn
  tMsOut
  tMultiApplyCal
  tPhaseshiftPredict
  tPredict
  tReadOnly
  tSplit)

if(LIBDIRAC_FOUND)
  add_python_tests(tSagecalPredict)
endif()
