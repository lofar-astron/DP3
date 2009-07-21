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
#include <casa/BasicMath/Math.h>
#include <casa/Arrays.h>
#include <casa/Quanta/MVEpoch.h>
#include <casa/Utilities/Assert.h>

#include <iostream>

#include <DPPP/MsFile.h>
#include <DPPP/MsInfo.h>
#include <DPPP/RunDetails.h>
#include <DPPP/TimeBuffer.h>

using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  MsFile::MsFile  <<<===============

MsFile::MsFile(const std::string& msin, const std::string& msout):
  SELECTblock(20),
  InMS(NULL),
  OutMS(NULL)
{
  InName  = msin;
  OutName = msout;
/* tableCommand(string("SELECT UVW,FLAG_CATEGORY,WEIGHT,SIGMA,ANTENNA1,ANTENNA2,ARRAY_ID,DATA_DESC_ID,") +
                string("EXPOSURE,FEED1,FEED2,FIELD_ID,FLAG_ROW,INTERVAL,OBSERVATION_ID,PROCESSOR_ID,") +
                string("SCAN_NUMBER,STATE_ID,TIME,TIME_CENTROID FROM $1"),
                Data_iter.table());*/

  SELECTblock[0]  = "UVW";
  SELECTblock[1]  = "FLAG_CATEGORY";
  SELECTblock[2]  = "WEIGHT";
  SELECTblock[3]  = "SIGMA";
  SELECTblock[4]  = "ANTENNA1";
  SELECTblock[5]  = "ANTENNA2";
  SELECTblock[6]  = "ARRAY_ID";
  SELECTblock[7]  = "DATA_DESC_ID";
  SELECTblock[8]  = "EXPOSURE";
  SELECTblock[9]  = "FEED1";
  SELECTblock[10] = "FEED2";
  SELECTblock[11] = "FIELD_ID";
  SELECTblock[12] = "FLAG_ROW";
  SELECTblock[13] = "INTERVAL";
  SELECTblock[14] = "OBSERVATION_ID";
  SELECTblock[15] = "PROCESSOR_ID";
  SELECTblock[16] = "SCAN_NUMBER";
  SELECTblock[17] = "STATE_ID";
  SELECTblock[18] = "TIME";
  SELECTblock[19] = "TIME_CENTROID";
}

//===============>>>  MsFile::~MsFile  <<<===============

MsFile::~MsFile()
{
  delete InMS;
  delete OutMS;
}

//===============>>>  DataSquasher::TableResize  <<<===============

void MsFile::TableResize(TableDesc tdesc, IPosition ipos, string name, Table& table)
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

//===============>>> MsFile::Init  <<<===============
bool MsFile::Init(MsInfo& Info, RunDetails& Details, int Squashing)
{
  cout << "Please wait, preparing output MS" << endl;
  Block<String> tempblock(SELECTblock);
  tempblock.resize(21);
  tempblock[20] = "FLAG";
  Table temptable = Table(InName).project(tempblock);
  // NOT copying WEIGHT_SPECTRUM as it only contains dummy data anyway
  // We do need FLAG to make it a valid MS
  Table outtable = TableCopy::makeEmptyTable(OutName, Record(), temptable,
                                             Table::NewNoReplace, Table::AipsrcEndian, true, true);
  TableCopy::copyInfo(outtable, temptable);
  TableCopy::copySubTables(outtable, temptable);

  InMS    = new MeasurementSet(InName);
  OutMS   = new MeasurementSet(OutName, Table::Update);
  //some magic to create a new DATA column
  TableDesc tdesc = InMS->tableDesc();
  ColumnDesc desc = tdesc.rwColumnDesc("DATA");
  IPosition ipos  = desc.shape();
  Vector<Int> temp_pos = ipos.asVector();
  std::cout << "Old shape: " << temp_pos(0) << ":" <<  temp_pos(1) << std::endl;
  int old_nchan = temp_pos(1);
  int new_nchan;
  if (Squashing)
  { new_nchan     = Details.NChan/Details.Step;
    temp_pos(1)   = new_nchan;
  }
  else
  { new_nchan     = old_nchan;
    Details.NChan = old_nchan;
    Details.Step  = 1;
    Details.Start = 0;
  }
  std::cout << "New shape: " << temp_pos(0) << ":" <<  temp_pos(1) << std::endl;
  IPosition data_ipos(temp_pos);

  if (tdesc.isColumn("WEIGHT_SPECTRUM"))
  { tdesc.removeColumn("WEIGHT_SPECTRUM");
  }
  tdesc.addColumn(ArrayColumnDesc<Float>("WEIGHT_SPECTRUM", "Added by datasquasher",
                                          data_ipos, ColumnDesc::FixedShape));

  TableResize(tdesc, data_ipos, "DATA", *OutMS);
  TableResize(tdesc, data_ipos, "WEIGHT_SPECTRUM", *OutMS);

  //if present handle the CORRECTED_DATA column
  if (tdesc.isColumn("CORRECTED_DATA"))
  {
    if (Details.Columns)
    {
      cout << "CORRECTED_DATA detected for processing" << endl;
      TableResize(tdesc, data_ipos, "CORRECTED_DATA", *OutMS);
    }
  }

  //if present handle the MODEL_DATA column
  if (tdesc.isColumn("MODEL_DATA"))
  {
    if (Details.Columns)
    {
      cout << "MODEL_DATA detected for processing" << endl;
      desc = tdesc.rwColumnDesc("MODEL_DATA");
      desc.setOptions(0);
      desc.setShape(data_ipos);
      desc.setOptions(4);
      desc.rwKeywordSet().removeField("CHANNEL_SELECTION"); //messes with the Imager if it's there but has wrong values
      Matrix<Int> selection;
      selection.resize(2, Info.NumBands);
      selection.row(0) = 0; //start in Imager, will therefore only work if imaging whole SPW
      selection.row(1) = new_nchan;
      desc.rwKeywordSet().define("CHANNEL_SELECTION", selection); // #spw x [startChan, NumberChan] for the VisBuf in the Imager
      // see code/msvis/implement/MSVis/VisSet.cc
      OutMS->addColumn(desc);
      OutMS->addColumn(ArrayColumnDesc<Float>("IMAGING_WEIGHT","imaging weight", 1));
    }
  }

  //fix the FLAGS column
  TableResize(tdesc, data_ipos, "FLAG", *OutMS);

  //Fix the SpectralWindow values
  IPosition spw_ipos(1,new_nchan);
  MSSpectralWindow inSPW = InMS->spectralWindow();
  //ugly workaround MSSpectral window does no allow deleting and then recreating columns
  Table outSPW = Table(OutName + "/SPECTRAL_WINDOW", Table::Update);
  ScalarColumn<Int> channum(outSPW, "NUM_CHAN");
  channum.fillColumn(new_nchan);

  TableDesc SPWtdesc = inSPW.tableDesc();
  TableResize(SPWtdesc, spw_ipos, "CHAN_FREQ", outSPW);

  TableResize(SPWtdesc, spw_ipos, "CHAN_WIDTH", outSPW);

  TableResize(SPWtdesc, spw_ipos, "EFFECTIVE_BW", outSPW);

  TableResize(SPWtdesc, spw_ipos, "RESOLUTION", outSPW);

  ROArrayColumn<Double> inFREQ(inSPW, "CHAN_FREQ");
  ROArrayColumn<Double> inWIDTH(inSPW, "CHAN_WIDTH");
  ROArrayColumn<Double> inBW(inSPW, "EFFECTIVE_BW");
  ROArrayColumn<Double> inRESOLUTION(inSPW, "RESOLUTION");

  ArrayColumn<Double> outFREQ(outSPW, "CHAN_FREQ");
  ArrayColumn<Double> outWIDTH(outSPW, "CHAN_WIDTH");
  ArrayColumn<Double> outBW(outSPW, "EFFECTIVE_BW");
  ArrayColumn<Double> outRESOLUTION(outSPW, "RESOLUTION");

  Vector<Double> old_temp(old_nchan, 0.0);
  Vector<Double> new_temp(new_nchan, 0.0);

  for (unsigned int i = 0; i < inSPW.nrow(); i++)
  {
    for (int j = 0; j < new_nchan; j++)
    { inFREQ.get(i, old_temp);
      if (Details.Step % 2) //odd number of channels in step
      { new_temp(j) = old_temp(Details.Start + j*Details.Step + (Details.Step + 1)/2);
      }
      else //even number of channels in step
      { new_temp(j) = 0.5 * (old_temp(Details.Start + j*Details.Step + Details.Step/2 -1)
                              + old_temp(Details.Start + j*Details.Step + Details.Step/2));
      }
      outFREQ.put(i, new_temp);
    }
    for (int j = 0; j < new_nchan; j++)
    { inWIDTH.get(i, old_temp);
      new_temp(j) = old_temp(0) * Details.Step;
      outWIDTH.put(i, new_temp);
    }
    for (int j = 0; j < new_nchan; j++)
    { inBW.get(i, old_temp);
      new_temp(j) = old_temp(0) * Details.Step;
      outBW.put(i, new_temp);
    }
    for (int j = 0; j < new_nchan; j++)
    { inRESOLUTION.get(i, old_temp);
      new_temp(j) = old_temp(0) * Details.Step;
      outRESOLUTION.put(i, new_temp);
    }
  }
  OutMS->flush(true);
  cout << "Finished preparing output MS" << endl;
  return (tdesc.isColumn("CORRECTED_DATA") && tdesc.isColumn("MODEL_DATA"));
}

//===============>>> MsFile::PrintInfo  <<<===============
void MsFile::PrintInfo(void)
{
  std::cout << "In  MeasurementSet:   " << InName << std::endl;
  std::cout << "Out MeasurementSet:   " << OutName << std::endl;
}

//===============>>> MsFile::TimeIterator  <<<===============

TableIterator MsFile::TimeIterator()
{
  Block<String> ms_iteration_variables(1);
  ms_iteration_variables[0] = "TIME";

  return TableIterator((*InMS), ms_iteration_variables, TableIterator::Ascending);
}

//===============>>> MsFile::UpdateTimeslotData  <<<===============
void MsFile::UpdateTimeslotData(casa::TableIterator& Data_iter,
                                MsInfo& Info,
                                DataBuffer& Buffer,
                                TimeBuffer& TimeData)
{
  Table         TimeslotTable = Data_iter.table();
  int           rowcount      = TimeslotTable.nrow();
  bool          columns       = Buffer.ModelData.size() > 0;
  ROTableVector<Int>            antenna1      (TimeslotTable, "ANTENNA1");
  ROTableVector<Int>            antenna2      (TimeslotTable, "ANTENNA2");
  ROTableVector<Int>            bandnr        (TimeslotTable, "DATA_DESC_ID");
  ROArrayColumn<Complex>        data          (TimeslotTable, "DATA");
  ROArrayColumn<Double>         uvw           (TimeslotTable, "UVW");
  ROTableVector<Double>         time_centroid (TimeslotTable, "TIME_CENTROID");
  ROTableVector<Double>         time          (TimeslotTable, "TIME");
  ROTableVector<Double>         interval      (TimeslotTable, "INTERVAL");
  ROTableVector<Double>         exposure      (TimeslotTable, "EXPOSURE");
  ROArrayColumn<Bool>           flags         (TimeslotTable, "FLAG");
  ROArrayColumn<Complex>        modeldata;
  ROArrayColumn<Complex>        correcteddata;
  if (columns)
  { modeldata.attach(TimeslotTable, "MODEL_DATA");
    correcteddata.attach(TimeslotTable, "CORRECTED_DATA");
  }
  Cube<Complex>                 tempData(Info.NumPolarizations, Info.NumChannels, rowcount);
  Cube<Complex>                 tempModelData(Info.NumPolarizations, Info.NumChannels, rowcount);
  Cube<Complex>                 tempCorrectedData(Info.NumPolarizations, Info.NumChannels, rowcount);
  Cube<Bool>                    tempFlags(Info.NumPolarizations, Info.NumChannels, rowcount);

  data.getColumn(tempData); //We're not checking Data.nrow() Data.ncolumn(), assuming all data is the same size.
  if (columns)
  { modeldata.getColumn(tempModelData);
    correcteddata.getColumn(tempCorrectedData);
  }
  flags.getColumn(tempFlags);
  TimeData.Time.push_back(time(0));
  TimeData.TimeCentroid.push_back(time_centroid(0));
  TimeData.Interval.push_back(interval(0));
  TimeData.Exposure.push_back(exposure(0));
  TimeData.Uvw.push_back(uvw(0));

  for (int i = 0; i < rowcount; i++)
  {
    int bi    = Info.BaselineIndex[baseline_t(antenna1(i), antenna2(i))];
    int band  = bandnr(i);
    int index = (band % Info.NumBands) * Info.NumPairs + bi;
    Buffer.Data[index].xyPlane(Buffer.Position)  = tempData.xyPlane(i);
    Buffer.Flags[index].xyPlane(Buffer.Position) = tempFlags.xyPlane(i);
    if (columns)
    { Buffer.ModelData[index].xyPlane(Buffer.Position)     = tempData.xyPlane(i);
      Buffer.CorrectedData[index].xyPlane(Buffer.Position) = tempData.xyPlane(i);
    }
  }
}

//===============>>> MsFile::WriteFlags  <<<===============

void MsFile::WriteData(casa::TableIterator& Data_iter,
                       MsInfo& Info,
                       DataBuffer& Buffer,
                       TimeBuffer& TimeData)
{
  Table DataTable = *OutMS;
  int   rowcount  = Data_iter.table().nrow();
  int   nrows     = DataTable.nrow();
  int   pos       = (Buffer.Position+1) % Buffer.WindowSize;
  bool  columns   = Buffer.ModelData.size() > 0;
  Table temptable = Data_iter.table().project(SELECTblock);

  DataTable.addRow(rowcount);
  Table dummy = DataTable.project(SELECTblock);
  TableCopy::copyRows(dummy, temptable, nrows, 0, rowcount);
  ROTableVector<Int>        antenna1     (DataTable, "ANTENNA1");
  ROTableVector<Int>        antenna2     (DataTable, "ANTENNA2");
  ROTableVector<Int>        bandnr       (DataTable, "DATA_DESC_ID");
  TableVector<Double>       time         (DataTable, "TIME");
  TableVector<Double>       time_centroid(DataTable, "TIME_CENTROID");
  TableVector<Double>       exposure     (DataTable, "EXPOSURE");
  TableVector<Double>       interval     (DataTable, "INTERVAL");
  ArrayColumn  <Double>     uvw          (DataTable, "UVW");
  ArrayColumn  <Complex>    data         (DataTable, "DATA");
  ArrayColumn  <Bool>       flags        (DataTable, "FLAG");
  ArrayColumn  <Float>      weights      (DataTable, "WEIGHT_SPECTRUM");
  ArrayColumn  <Complex>    modeldata;
  ArrayColumn  <Complex>    correcteddata;
  if (columns)
  { modeldata.attach(DataTable, "MODEL_DATA");
    correcteddata.attach(DataTable, "CORRECTED_DATA");
  }
  //cout << "Processing: " << MVTime(temp(0)/(24*3600)).string(MVTime::YMD) << endl; //for testing purposes

  TimeData.Squash();

  for (int i = 0; i < rowcount; i++)
  {
    int bi    = Info.BaselineIndex[baseline_t(antenna1(i), antenna2(i))];
    int band  = bandnr(i);
    int index = (band % Info.NumBands) * Info.NumPairs + bi;

    data.put(nrows + i, Buffer.Data[index].xyPlane(pos));
    flags.put(nrows + i, Buffer.Flags[index].xyPlane(pos));
    weights.put(nrows + i, Buffer.Weights[index].xyPlane(pos));
    time.set(nrows + i, TimeData.Time[0]);
    time_centroid.set(nrows + i, TimeData.TimeCentroid[0]);
    exposure.set(nrows + i, TimeData.Exposure[0]);
    interval.set(nrows + i, TimeData.Interval[0]);
    if (columns)
    {
      modeldata.put(nrows + i, Buffer.ModelData[index].xyPlane(pos));
      correcteddata.put(nrows + i, Buffer.CorrectedData[index].xyPlane(pos));
    }
  }
  TimeData.Clear();
}

//===============>>> MsFile  <<<===============
