/***************************************************************************
 *   Copyright (C) 2007-8 by ASTRON, Adriaan Renting                       *
 *   renting@astron.nl                                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <lofar_config.h>
#include <iostream>

#include <CS1_pp_lib/RunDetails.h>

using namespace LOFAR::CS1;

//===============>>>  RunDetails::RunDetails  <<<===============

RunDetails::RunDetails()
{
}

//===============>>>  RunDetails::~RunDetails  <<<===============

RunDetails::~RunDetails()
{
}

//===============>>> RunDetails::PrintInfo  <<<===============

bool RunDetails::CheckValues(void)
{
  bool result = false;
  if (!((FreqWindow %2) && (TimeWindow % 2)))
  { std::cout << "The window sizes need to be uneven" << std::endl;
    result = true;
  }
  return result;
}
//===============>>> RunDetails::PrintInfo  <<<===============

void RunDetails::PrintInfo(void)
{
  std::cout << "Fixed:             " << Fixed << std::endl;        // BandpassCorrector
  std::cout << "FreqWindow:        " << FreqWindow << std::endl;   // FrequencyFlagger, MADFlagger
  std::cout << "TimeWindow:        " << TimeWindow << std::endl;   // ComplexMedianFlagger, MADFlagger
  std::cout << "Threshold:         " << Threshold << std::endl;     // FrequencyFlagger
  std::cout << "MinThreshold:      " << MinThreshold << std::endl; // ComplexMedianFlagger
  std::cout << "MaxThreshold:      " << MaxThreshold << std::endl; // ComplexMedianFlagger
  std::cout << "Algorithm:         " << Algorithm << std::endl;    // FrequencyFlagger
  std::cout << "Existing:          " << Existing << std::endl;     // all flaggers
  std::cout << "NChan:             " << NChan << std::endl;        // DataSquasher
  std::cout << "Start:             " << Start << std::endl;        // DataSquasher
  std::cout << "Step:              " << Step << std::endl;         // DataSquasher
  std::cout << "Skip:              " << Skip << std::endl;         // DataSquasher
  std::cout << "Columns:           " << Columns << std::endl;      // DataSquasher
  std::cout << "TimeStep:          " << TimeStep << std::endl;     // DataSquasher
}

//===============>>> RunDetails  <<<===============

