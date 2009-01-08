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
#include <iostream>
#include <cstdlib>
#include <string>
#include <limits>
#include <tables/Tables.h>
#include <tables/Tables/TableParse.h>
#include <ms/MeasurementSets.h>
#include <casa/Exceptions.h>

#include <CS1_SPWCombine/CombinerProcessControl.h>
#include "SPWCombine.h"

#define COMBINER_VERSION "0.22"
// 0.10 Initial version based on DataSquasher
// 0.20 Ported additions and updates from DataSquasher
// 0.21 Fixed calculation of REF_FREQUENCY
// 0.22 Added handing of Measurementsets with different numbers of timesamples

namespace LOFAR
{
  namespace CS1
  {
    using namespace casa;
    //===============>>> CombinerProcessControl::CombinerProcessControl  <<<===============
    CombinerProcessControl::CombinerProcessControl()
    : ProcessControl()
    {
      itsCombiner = NULL;
    }

    //===============>>> CombinerProcessControl::~CombinerProcessControl  <<<==============
    CombinerProcessControl::~CombinerProcessControl()
    {
    }

    //===============>>> CombinerProcessControl::define  <<<==============================
    tribool CombinerProcessControl::define()
    {
      ParameterSet* ParamSet = globalParameterSet();
      itsInMS  = ParamSet->getStringVector("inms");
      itsOutMS = ParamSet->getString("outms");
      return true;
    }

    //===============>>> CombinerProcessControl::run  <<<=================================
    tribool CombinerProcessControl::run()
    {
      try{
        std::cout << "Creating " << itsOutMS << ", please wait..." << std::endl;
        unsigned int min_nrow = std::numeric_limits<unsigned int>::max();
        int sourceMS = 0;
        for (unsigned int i = 0; i < itsInMS.size(); i++) // to search the shortest MS
        {
          if (inMS[i]->nrow() < min_nrow)
          {
            min_nrow = inMS[i]->nrow();
            sourceMS = i;
          }
        }
        if (sourceMS > 0)
        { cout << "Not all sources are the same lenght, using te shortest one." << endl;
        }

        Table TempTable = tableCommand(string("SELECT UVW,FLAG_CATEGORY,WEIGHT,SIGMA,ANTENNA1,ANTENNA2,ARRAY_ID,DATA_DESC_ID,") +
                                       string("EXPOSURE,FEED1,FEED2,FIELD_ID,FLAG_ROW,INTERVAL,OBSERVATION_ID,PROCESSOR_ID,") +
                                       string("SCAN_NUMBER,STATE_ID,TIME,TIME_CENTROID,WEIGHT_SPECTRUM,FLAG FROM ")
                                       + itsInMS[sourceMS] + string(" WHERE DATA_DESC_ID = 0"));
        // Need FLAG to make it a valid MS
        TempTable.deepCopy(itsOutMS, Table::NewNoReplace, true);
        tableCommand(string("DELETE FROM ") + itsOutMS + string("/DATA_DESCRIPTION WHERE rownumber() > 1"));
        tableCommand(string("DELETE FROM ") + itsOutMS + string("/SPECTRAL_WINDOW WHERE rownumber() > 1"));

        MeasurementSet outMS = MeasurementSet(itsOutMS, Table::Update);
        int nchan = 0;
        for (unsigned int i = 0; i < itsInMS.size(); i++)
        {
          itsCombiner->GetMSInfo(*(inMS[i]));
          for (int j = 0; j < itsCombiner->itsNumBands; j++)
          { nchan += itsCombiner->itsNumChannels;
          }
        }
        TableDesc tdesc    = inMS[sourceMS]->tableDesc();
        Vector<Int> temp(2);
        temp(0)            = itsCombiner->itsNumPolarizations;
        temp(1)            = nchan;
        std::cout << "New number of channels: " << nchan << std::endl;
        IPosition data_ipos(temp);

        itsCombiner->TableResize(tdesc, data_ipos, "DATA", outMS);

        //fix the FLAGS column
        itsCombiner->TableResize(tdesc, data_ipos, "FLAG", outMS);

        //Fix the SpectralWindow values
        IPosition spw_ipos(1, nchan);
        //ugly workaround MSSpectral window does no allow deleting and then recreating columns
        Table outSPW = Table(itsOutMS + "/SPECTRAL_WINDOW", Table::Update);
        ScalarColumn<Int> channum(outSPW, "NUM_CHAN");
        channum.fillColumn(nchan);

        TableDesc SPWtdesc = outSPW.tableDesc();
        itsCombiner->TableResize(SPWtdesc, spw_ipos, "CHAN_FREQ", outSPW);
        itsCombiner->TableResize(SPWtdesc, spw_ipos, "CHAN_WIDTH", outSPW);
        itsCombiner->TableResize(SPWtdesc, spw_ipos, "EFFECTIVE_BW", outSPW);
        itsCombiner->TableResize(SPWtdesc, spw_ipos, "RESOLUTION", outSPW);

        ArrayColumn<Double> outFREQ(outSPW, "CHAN_FREQ");
        ArrayColumn<Double> outWIDTH(outSPW, "CHAN_WIDTH");
        ArrayColumn<Double> outBW(outSPW, "EFFECTIVE_BW");
        ArrayColumn<Double> outRESOLUTION(outSPW, "RESOLUTION");
        ScalarColumn<Double> outREF_FREQUENCY(outSPW, "REF_FREQUENCY");

        Vector<Double> new_FREQ(nchan, 0.0);
        Vector<Double> new_WIDTH(nchan, 0.0);
        Vector<Double> new_BW(nchan, 0.0);
        Vector<Double> new_RESOLUTION(nchan, 0.0);
        int total_channels   = 0;
        int total_bands      = 0;
        double ref_frequency = 0.0;

        for (unsigned int i = 0; i < itsInMS.size(); i++)
        {
          itsCombiner->GetMSInfo(*(inMS[i]));
          int old_nchan = itsCombiner->itsNumChannels;
          Vector<Double> old_temp(old_nchan, 0.0);

          MSSpectralWindow inSPW = inMS[i]->spectralWindow();

          ROArrayColumn<Double> inFREQ(inSPW, "CHAN_FREQ");
          ROArrayColumn<Double> inWIDTH(inSPW, "CHAN_WIDTH");
          ROArrayColumn<Double> inBW(inSPW, "EFFECTIVE_BW");
          ROArrayColumn<Double> inRESOLUTION(inSPW, "RESOLUTION");
          ROScalarColumn<Double> inREF_FREQUENCY(inSPW, "REF_FREQUENCY");

          for (unsigned int n = 0; n < inSPW.nrow(); n++)
          {
            for (int m = 0; m < old_nchan; m++)
            {
            inFREQ.get(n, old_temp); // could be outsid this loop
            new_FREQ(total_channels + m) = old_temp(m);

            inWIDTH.get(n, old_temp);
            new_WIDTH(total_channels + m) = old_temp(m);

            inBW.get(n, old_temp);
            new_WIDTH(total_channels + m) = old_temp(m);

            inRESOLUTION.get(n, old_temp);
            new_RESOLUTION(total_channels + m) = old_temp(m);
            }
            total_channels += old_nchan;

            double temp_freq;
            inREF_FREQUENCY.get(n, temp_freq);
            ref_frequency += temp_freq;
            total_bands++;
          }
          outFREQ.put(0, new_FREQ);
          outWIDTH.put(0, new_WIDTH);
          outBW.put(0, new_BW);
          outRESOLUTION.put(0, new_RESOLUTION);
          outREF_FREQUENCY.put(0, ref_frequency/total_bands);
        }

        //Do the real stuff
        itsCombiner->Combine(inMS, outMS, "DATA");
      }
      catch(casa::AipsError& err)
      {
        std::cerr << "Aips++ error detected: " << err.getMesg() << std::endl;
        return false;
      }
      return true;
    }

    //===============>>> CombinerProcessControl::init  <<<================================
    tribool CombinerProcessControl::init()
    {
      try {
      using std::cout;
      using std::cerr;
      using std::endl;

      cout  << string(COMBINER_VERSION) + string(" spw combine by Adriaan Renting for LOFAR CS1\n") +
              string("This is experimental software, please report errors or requests to renting@astron.nl\n") +
              string("Documentation can be found at: www.lofar.org/operations/doku.php?id=engineering:software:postprocessing_software\n");
      cout << string("Combining ");
      for (unsigned int i = 0; i < itsInMS.size(); i++)
      {
        cout << itsInMS[i] << ", ";
      }
      cout << string(" into ") << itsOutMS << endl;
      if (itsInMS.size() == 0 || itsOutMS == "")
      {
        cerr  << " Error missing input" << endl;
        return false;
      }
      inMS.resize(itsInMS.size());
      for (unsigned int i = 0; i < itsInMS.size(); i++)
      {
        inMS[i] = new MeasurementSet(itsInMS[i]);
      }
      itsCombiner = new SPWCombine();
      }
      catch(casa::AipsError& err)
      {
        std::cerr << "Aips++ error detected: " << err.getMesg() << std::endl;
        return false;
      }
      return true;
    }

    //===============>>> CombinerProcessControl::pause  <<<===============================
    tribool CombinerProcessControl::pause(const std::string&)
    { return false;
    }

    //===============>>> CombinerProcessControl::quit  <<<================================
    tribool CombinerProcessControl::quit()
    {
      return true;
    }

    //===============>>> CombinerProcessControl::release  <<<=============================
    tribool CombinerProcessControl::release()
    { return false;
    }
    //===============>>> CombinerProcessControl::recover  <<<=============================
    tribool CombinerProcessControl::recover(const std::string&)
    { return false;
    }

    //===============>>> CombinerProcessControl::reinit  <<<==============================
    tribool CombinerProcessControl::reinit(const  std::string&)
    { return false;
    }

    //===============>>> CombinerProcessControl::askInfo  <<<=============================
    std::string CombinerProcessControl::askInfo(const std::string&)
    { return std::string("");
    }

    //===============>>> CombinerProcessControl::snapshot  <<<============================
    tribool CombinerProcessControl::snapshot(const std::string&)
    { return false;
    }
  } //namespace CS1
}; //namespace LOFAR
