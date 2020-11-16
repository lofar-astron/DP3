// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @author Adriaan Renting

#ifndef LOFARSPWCOMBINE_H
#define LOFARSPWCOMBINE_H

#include <iostream>
#include <cstdlib>
#include <string>
#include <ms/MeasurementSets.h>
#include <casa/Arrays.h>
#include <tables/Tables.h>
#include <tables/Tables/TableIter.h>

namespace LOFAR {
namespace CS1 {
using namespace casa;

class SPWCombine {
 private:
  int itsNumAntennae;
  int itsNumPairs;

  TableIterator CreateDataIterator(MeasurementSet& myMS);

 public:
  SPWCombine(void);
  ~SPWCombine(void);
  int itsNumPolarizations;
  int itsNumChannels;
  int itsNumBands;

  void TableResize(TableDesc tdesc, IPosition ipos, std::string name,
                   Table& table);
  void GetMSInfo(MeasurementSet& myMS);
  void Combine(vector<MeasurementSet*> inMS, MeasurementSet& myMS,
               std::string Data);
};
}  // namespace CS1
};  // namespace LOFAR
#endif
