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

void MsFile::TableResize(ColumnDesc desc, const IPosition& ipos,
                         TiledColumnStMan* tsm, Table& table)
{
  desc.setOptions(0);
  desc.setShape(ipos);
  desc.setOptions(ColumnDesc::FixedShape);
  if (table.tableDesc().isColumn(desc.name())) {
    table.removeColumn(desc.name());
  }
  // Use tiled storage manager if given.
  if (tsm == 0) {
    table.addColumn (desc);
  } else {
    table.addColumn (desc, *tsm);
  }
}

//===============>>> MsFile::DetermineDATAshape  <<<===============

IPosition MsFile::DetermineDATAshape(const Table& MS)
{
  ROArrayColumn<Complex> temp_column(MS, "DATA");
  // First try to get it as a fixed column shape.
  IPosition shp = temp_column.shapeColumn();
  if (shp.empty()  &&  MS.nrow() > 0) {
    // Not fixed shape, so get shape from first row (if present).
    shp = temp_column.shape(0);
  }
  // Still unknown, so use default [4,1].
  if (shp.empty()) {
    THROW(PipelineException, "Error, can't figure out the shape of the DATA column");
  }
  return shp;
}

//===============>>> MsFile::Init  <<<===============
void MsFile::Init(MsInfo& Info, RunDetails& Details, int Squashing)
{
  std::cout << "Preparing output MS " << OutName << std::endl;
  // Open the MS and obtain the description.
  InMS = new MeasurementSet(InName);
  TableDesc tdesc = InMS->tableDesc();
  // Determine the output data shape.
  const ColumnDesc& desc = tdesc.columnDesc("DATA");
  IPosition data_ipos = DetermineDATAshape(*InMS);
  std::cout << "Old shape: " << data_ipos[0] << ":" <<  data_ipos[1] << std::endl;
  int old_nchan = data_ipos[1];
  int new_nchan = old_nchan;
  if (Squashing)
    { new_nchan     = Details.NChan/Details.Step;
      data_ipos[1]  = new_nchan;
    }
  else
    { Details.NChan = old_nchan;
      Details.Step  = 1;
      Details.Start = 0;
    }
  std::cout << "New shape: " << data_ipos[0] << ":" <<  data_ipos[1] << std::endl;

  // Create the output table without the data columns.
  Table temptable = InMS->project(SELECTblock);
  TableDesc tempdesc = temptable.tableDesc();
  // Remove possible hypercolumn definitions.
  tempdesc.adjustHypercolumns (SimpleOrderedMap<String,String>(String()));
  Record dminfo = temptable.dataManagerInfo();
  // Use TSM for UVW.
  TableCopy::setTiledStMan (dminfo, Vector<String>(1, "UVW"),
                            "TiledColumnStMan", "TiledUVW", IPosition(2, 3, 1024));
  // Replace all non-writable storage managers by SSM.
  dminfo = TableCopy::adjustStMan (dminfo);
  SetupNewTable newtab(OutName, tempdesc, Table::NewNoReplace);
  newtab.bindCreate (dminfo);
  Table outtable(newtab);
  {
    // Add DATA column using tsm.
    // Use a tilesize of 32 KB (4*8*128 values * 8 bytes).
    TiledColumnStMan tsm("TiledData", IPosition(3, data_ipos[0], 8, 128));

    TableResize(desc, data_ipos, &tsm, outtable);
  }
  {
    // Add FLAG column using tsm.
    TiledColumnStMan tsmf("TiledFlag", IPosition(3, data_ipos[0], 8, 128*8));
    TableResize(tdesc["FLAG"], data_ipos, &tsmf, outtable);
  }
  {
    // Add WEIGHT_SPECTRUM column using tsm.
    TiledColumnStMan tsmw("TiledWeightSpec", IPosition(3, data_ipos[0], 8, 128));
    TableResize(tdesc["WEIGHT_SPECTRUM"], data_ipos, &tsmw, outtable);
  }
  // If both present handle the CORRECTED_DATA and MODEL_DATA column.
  if (Details.Columns)
  {
    if (tdesc.isColumn("CORRECTED_DATA") && tdesc.isColumn("MODEL_DATA"))
    {
      cout << "MODEL_DATA detected for processing" << endl;
      ColumnDesc mdesc = tdesc.columnDesc("MODEL_DATA");
      TableRecord& keyset = mdesc.rwKeywordSet();
      // Redefine possible keywords used by the CASA VisSet classes (in imager).
      if (keyset.isDefined("CHANNEL_SELECTION")) {
        keyset.removeField("CHANNEL_SELECTION");
      }
      Matrix<Int> selection(2, Info.NumBands);
      selection.row(0) = 0;
      selection.row(1) = new_nchan;
      keyset.define("CHANNEL_SELECTION", selection);
      TiledColumnStMan tsmm("ModelData", IPosition(3, data_ipos[0], 8, 128));
      TableResize(mdesc, data_ipos, &tsmm, outtable);

      cout << "CORRECTED_DATA detected for processing" << endl;
      TiledColumnStMan tsmc("CorrectedData", IPosition(3, data_ipos[0], 8, 128));
      TableResize(tdesc["CORRECTED_DATA"], data_ipos, &tsmc, outtable);

      TiledColumnStMan tsmw("TiledWeight", IPosition(3, data_ipos[0], 8, 128));
      TableResize(tdesc["IMAGING_WEIGHT"], data_ipos, &tsmw, outtable);
    }
    else
    {
      if (tdesc.isColumn("CORRECTED_DATA") || tdesc.isColumn("MODEL_DATA"))
      {
        cout << "Only one of CORRECTED_DATA and MODEL_DATA columns is present; "
             << "it is ignored" << endl;
      }
      Details.Columns = false;
    }
  }
  cout << " copying info and subtables ..." << endl;
  // Copy the info and subtables.
  TableCopy::copyInfo(outtable, temptable);
  TableCopy::copySubTables(outtable, temptable);

  // All columns are present, so now it can be opened as an MS.
  outtable.flush();
  OutMS = new MeasurementSet(OutName, Table::Update);

  //Fix the SpectralWindow values
  IPosition spw_ipos(1,new_nchan);
  MSSpectralWindow inSPW = InMS->spectralWindow();
  //ugly workaround MSSpectral window does no allow deleting and then recreating columns
  Table outSPW = Table(OutName + "/SPECTRAL_WINDOW", Table::Update);
  ScalarColumn<Int> channum(outSPW, "NUM_CHAN");
  channum.fillColumn(new_nchan);

  TableDesc SPWtdesc = inSPW.tableDesc();
  TableResize(SPWtdesc["CHAN_FREQ"], spw_ipos, 0, outSPW);

  TableResize(SPWtdesc["CHAN_WIDTH"], spw_ipos, 0, outSPW);

  TableResize(SPWtdesc["EFFECTIVE_BW"], spw_ipos, 0, outSPW);

  TableResize(SPWtdesc["RESOLUTION"], spw_ipos, 0, outSPW);

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
  // Usually a needless sort on TIME does not harm so much.
  // However, for LofarStMan is can be quite costly. As we know its data
  // is always in time order, we do not sort for LofarStMan.
  // Determine that by looking at the storage manager type.
  TableIterator::Option opt = TableIterator::QuickSort;
  Record dminfo = InMS->dataManagerInfo();
  for (unsigned i=0; i<dminfo.nfields(); ++i) {
    Record subrec = dminfo.subRecord(i);
    if (subrec.asString("TYPE") == "LofarStMan") {
      opt = TableIterator::NoSort;
      break;
    }
  }
  Block<String> ms_iteration_variables(1);
  ms_iteration_variables[0] = "TIME";

  return TableIterator((*InMS), ms_iteration_variables,
                       TableIterator::Ascending, opt);
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
    TimeData.BufTime[index].push_front(time(i));
    TimeData.BufTimeCentroid[index].push_front(time_centroid(i));
    TimeData.BufInterval[index].push_front(interval(i));
    TimeData.BufExposure[index].push_front(exposure(i));
    TimeData.BufUvw[index].push_front(uvw(i));
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
  { modeldata.attach(                     DataTable, "MODEL_DATA");
    correcteddata.attach(                 DataTable, "CORRECTED_DATA");
  }
  //cout << "Processing: " << MVTime(temp(0)/(24*3600)).string(MVTime::YMD) << endl; //for testing purposes

  for (int i = 0; i < rowcount; i++)
  {
    int bi    = Info.BaselineIndex[baseline_t(antenna1(i), antenna2(i))];
    int band  = bandnr(i);
    int index = (band % Info.NumBands) * Info.NumPairs + bi;

    data.put(nrows + i, Buffer.Data[index].xyPlane(pos));
    flags.put(nrows + i, Buffer.Flags[index].xyPlane(pos));
    time.set(nrows + i, TimeData.Time[index][0]);
    time_centroid.set(nrows + i, TimeData.TimeCentroid[index][0]);
    exposure.set(nrows + i, TimeData.Exposure[index][0]);
    interval.set(nrows + i, TimeData.Interval[index][0]);
    uvw.put(nrows + i, TimeData.Uvw[index][0]);
    weights.put(nrows + i, Buffer.Weights[index].xyPlane(pos));
    if (columns)
    {
      modeldata.put(nrows + i, Buffer.ModelData[index].xyPlane(pos));
      correcteddata.put(nrows + i, Buffer.CorrectedData[index].xyPlane(pos));
    }
  }
  TimeData.Clear();
}

//===============>>> MsFile  <<<===============
