/***************************************************************************
 *   Copyright (C) 2007 by Adriaan Renting, ASTRON                         *
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
#ifndef LOFARSPWCOMBINE_H
#define LOFARSPWCOMBINE_H

/**
@author Adriaan Renting
*/
#include <iostream>
#include <cstdlib>
#include <string>
#include <ms/MeasurementSets.h>
#include <casa/Arrays.h>
#include <tables/Tables.h>
#include <tables/Tables/TableIter.h>

namespace LOFAR
{
  namespace CS1
  {
    using namespace casa;

    class SPWCombine
    {
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

      void TableResize(TableDesc tdesc, IPosition ipos, std::string name, Table& table);
      void GetMSInfo(MeasurementSet& myMS);
      void Combine(vector<MeasurementSet*> inMS, MeasurementSet& myMS, std::string Data);
    };
  } //namespace CS1
}; //namespace LOFAR
#endif
