# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

configure_file(testInit.sh.in testInit.sh)

add_python_integrationtests(
  tApplyBeam
  tApplyCal2
  tBdaExpander
  tBdaPredict
  tDemix
  tGainCal
  tGainCalH5Parm
  tMsIn
  tMultiApplyCal
  tPredict
  tReadOnly
  tSplit)
