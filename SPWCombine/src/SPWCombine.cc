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
#include <lofar_config.h>
#include "SPWCombine.h"
namespace LOFAR
{
  namespace CS1
  {
    using namespace casa;

    //===============>>>  SPWCombine::SPWCombine  <<<===============

    SPWCombine::SPWCombine(void)
    {
    }

    //===============>>>  SPWCombine::~SPWCombine  <<<===============

    SPWCombine::~SPWCombine(void)
    {
    }

    //===============>>>  SPWCombine::GetMSInfo  <<<===============

    void SPWCombine::GetMSInfo(MeasurementSet& myMS)
    {
      //Number of channels in the Band
      MSSpectralWindow spectral_window = myMS.spectralWindow();
      ROScalarColumn<Int>              NUM_CHAN_col(spectral_window, "NUM_CHAN");
      itsNumChannels                   = NUM_CHAN_col(0);

      //Number of polarizations
      MSPolarization polarization      = myMS.polarization();
      ROScalarColumn<Int>              NUM_CORR_col(polarization, "NUM_CORR");
      itsNumPolarizations              = NUM_CORR_col(0);

      //Number of Bands
      itsNumBands                      = NUM_CHAN_col.nrow();

      //Number of Antennae
      MSAntenna antennae               = myMS.antenna();
      itsNumAntennae                   = antennae.nrow();

      //calculate number of baselines.
      itsNumPairs = (itsNumAntennae) * (itsNumAntennae + 1) / 2; //Triangular numbers formula
    }

    //===============>>>  SPWCombine::CreateDataIterator  <<<===============

    TableIterator SPWCombine::CreateDataIterator(MeasurementSet& myMS)
    {
      Block<String> ms_iteration_variables(4);
      ms_iteration_variables[0] = "TIME_CENTROID";
      ms_iteration_variables[1] = "DATA_DESC_ID";
      ms_iteration_variables[2] = "ANTENNA1";
      ms_iteration_variables[3] = "ANTENNA2";

      return TableIterator(myMS, ms_iteration_variables);
    }

    //===============>>>  SPWCombine::TableResize  <<<===============

    void SPWCombine::TableResize(TableDesc tdesc, IPosition ipos, string name, Table& table)
    {
      ColumnDesc desc = tdesc.rwColumnDesc(name);
      desc.setOptions(0);
      desc.setShape(ipos);
      desc.setOptions(4);
      if (table.tableDesc().isColumn(name))
      { table.removeColumn(name);
      }
      table.addColumn(desc);
    }

    //===============>>>  SPWCombine::Combine  <<<===============

    void SPWCombine::Combine(vector<MeasurementSet*> inMS, MeasurementSet& myMS, std::string Data)
    {
      int nrMS            = inMS.size();
      vector<int>           nrBands(nrMS);
      vector<int>           nrChannels(nrMS);
      vector<TableIterator> myIters(nrMS);

      for (int i = 0; i < nrMS; i++)
      {
        GetMSInfo(*(inMS[i]));
        nrBands[i]    = itsNumBands;
        nrChannels[i] = itsNumChannels;
        myIters[i]    = CreateDataIterator(*(inMS[i]));
      }

      TableIterator iter = CreateDataIterator(myMS);
      GetMSInfo(myMS);
      int step = myMS.nrow() / 10 + 1; //not exact but it'll do
      int row  = 0;
      while (!iter.pastEnd())
      {
        Cube<Bool>           myFlags(itsNumPolarizations, itsNumChannels, itsNumPairs);
        Cube<Complex>        myData(itsNumPolarizations, itsNumChannels, itsNumPairs);

        int bandMScounter = 0;
        for (int i = 0; i < nrMS; i++)
        {
          for (int j = 0; j < nrBands[i]; j++)
          {
            for (int k = 0; k < itsNumPairs; k++)
            {
              Table         oldTable = myIters[i].table();
              ROArrayColumn<Complex> Old(oldTable, Data);
              Matrix<Complex>        oldData(itsNumPolarizations, nrChannels[i]);
              Old.get(0, oldData);
              ROArrayColumn<Bool> OldFlags(oldTable, "FLAG");
              Matrix<Bool>        oldFlags(itsNumPolarizations, nrChannels[i]);
              OldFlags.get(0, oldFlags);
              (myIters[i])++;
              for (int m = 0; m < itsNumPolarizations; m++)
              {
                for (int n = 0; n < nrChannels[i]; n++)
                {
                  myData(m, bandMScounter + n, k)  = oldData(m, n);
                  myFlags(m, bandMScounter + n, k) = oldFlags(m, n);
                }
              }
            }
            bandMScounter += nrChannels[i];
          }
        }

        for (int i = 0; i < itsNumPairs; i ++)
        {
          if (row++ % step == 0) // to tell the user how much % we have processed,
          { std::cout << 10*(row/step) << "%" << std::endl; //not very accurate for low numbers of timeslots
          }

          Table                DataTable = iter.table();
          ArrayColumn<Bool>    Flags(DataTable, "FLAG");
          ArrayColumn<Complex> New(DataTable, Data);
          Flags.put(0, myFlags.xyPlane(i));
          New.put(0, myData.xyPlane(i));
          iter++;
        }
      }
    }
  } //namespace CS1
}; //namespace LOFAR
