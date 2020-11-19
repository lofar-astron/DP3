// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Adriaan Renting

#include "../Common/SystemUtil.h"

#include "../DPPP/Exceptions.h"

#include <PLC/ACCmain.h>

#include "CombinerProcessControl.h"

#include <iostream>

using namespace LOFAR;

// Use a terminate handler that can produce a backtrace.
// Exception::TerminateHandler t(Exception::terminate);

int main(int argc, char* argv[]) {
  try {
    DP3CS1::CombinerProcessControl myProcess;
    return DP3ACC::PLC::ACCmain(argc, argv, &myProcess);
  } catch (Exception& ex) {
    std::cerr << ex << std::endl;
    return 1;
  }
  return 0;
}
