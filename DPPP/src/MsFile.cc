//# Copyright (C) 2006-8
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$
//#
//# @author Adriaan Renting

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
#include <DPPP/Package__Version.h>
#include <Common/LofarLogger.h>

using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  MsFile::MsFile  <<<===============

MsFile::MsFile(const std::string& msin, const std::string& msout):
  SELECTblock(20),
  InName  (msin),
  OutName (msout),
  InMS    (NULL),
  OutMS   (NULL),
  itsHasWeightSpectrum(false),
  itsIsOrdered        (false)
{
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
  // Open the MS and obtain the description.
  // If no output MS, open it for update.
  InMS = new MeasurementSet(InName,
                            OutName.empty() ? Table::Update : Table::Old);
  // Test if WEIGHT_SPECTRUM is present.
  TableDesc tdesc = InMS->tableDesc();
  if (tdesc.isColumn("WEIGHT_SPECTRUM")) {
    // The column is there, but it might not contain values. Test row 0.
    itsHasWeightSpectrum = ROArrayColumn<Float>(*InMS, "WEIGHT_SPECTRUM").isDefined(0);
  }
  // Get the main table in TIME order.
  // Determine if stored using LofarStMan; If so, we know it is in time order.
  {
    Record dminfo = InMS->dataManagerInfo();
    for (unsigned i=0; i<dminfo.nfields(); ++i) {
      Record subrec = dminfo.subRecord(i);
      if (subrec.asString("TYPE") == "LofarStMan") {
        itsOrderedTable = *InMS;
        itsIsOrdered = true;
        break;
      }
    }
    // If not in time order, sort the main table.
    if (!itsIsOrdered) {
      itsOrderedTable = InMS->sort ("TIME");
    }
  }
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
  if (! itsIsOrdered) {
    checkGaps (Info, OutName.empty());
  }
  // Obtain the MS description.
  TableDesc tdesc = InMS->tableDesc();
  if (OutName.empty()) {
    std::cout << "No output MS; "
              << "flags will be updated in input MS" << std::endl;
    // The only DATA column to handle is the one used for flagging.
    // Check it exists.
    Details.DataColumns.resize (1);
    Details.DataColumns[0] = Details.FlagColumn;
    ASSERT (tdesc.isColumn(Details.FlagColumn));
    OutMS = new MeasurementSet(*InMS);
    // Write the parset info into the history.
    writeHistory (*OutMS, Details);
    // Create a TableIterator to have the correct input part.
    itsIterator = TimeIterator();
    return;
  }
  // A new MS will be created.
  std::cout << "Preparing output MS " << OutName << std::endl;
  // At least the DATA column needs to be handled.
  // All columns are needed if not flagging on DATA.
  Details.DataColumns.resize (1);
  Details.DataColumns[0] = "DATA";
  if (Details.FlagColumn != "DATA") {
    Details.AllColumns = True;
  }
  // Determine the output data shape.
  const ColumnDesc& desc = tdesc.columnDesc("DATA");
  IPosition data_ipos = DetermineDATAshape(*InMS);
  std::cout << "Old shape: " << data_ipos[0] << ":" <<  data_ipos[1] << std::endl;
  int old_nchan = data_ipos[1];
  int new_nchan = old_nchan;
  int nchan = Details.NChan;
  int nchanmax = old_nchan - Details.Start;
  if (nchan <= 0  ||  nchan > nchanmax) {
    nchan = nchanmax;
  }
  if (nchan <= 0) {
    THROW(PipelineException, "Error, no channels left (incorrect Start or NChan)");
  }
  if (Squashing)
    { Details.NChan = nchan;
      new_nchan     = (nchan + Details.Step - 1) / Details.Step;
      if (new_nchan == 1) {
        Details.Step = nchan;
      }
      // For the time being throw an exception if not divisible.
      // If not divisible, it works fine in the CHAN_FREQ averaging below,
      // but averaging in DataSquasher might fail.
      // I (Ger) do not have time to investigate that.
      // Another consideration: people may not like channels of unequal size.
      if (nchan != new_nchan*int(Details.Step)) {
        THROW(PipelineException, "Error: Step does not evenly divide nchan");
      }
      data_ipos[1]  = new_nchan;
    }
  else
    { Details.NChan = old_nchan;
      Details.Step  = 1;
      Details.Start = 0;
    }
  std::cout << "New shape: " << data_ipos[0] << ":" <<  data_ipos[1] << std::endl;

  // Check FreqWindow size.
  if (int(Details.FreqWindow) > old_nchan) {
    Details.FreqWindow = old_nchan;
  }

  // Create the output table without the data columns.
  Table temptable = InMS->project(SELECTblock);
  TableDesc tempdesc = temptable.tableDesc();
  // Remove possible hypercolumn definitions.
  tempdesc.adjustHypercolumns (SimpleOrderedMap<String,String>(String()));
  Record dminfo = temptable.dataManagerInfo();
  // Determine the DATA tile shape. Use all pols and the given #channels.
  // The nr of rows in a tile is determined by the given tile size (in kbytes).
  IPosition tileShape(3, data_ipos[0], Details.TileNChan, 1);
  tileShape[2] = Details.TileSize * 1024 / (8 * tileShape[0] * tileShape[1]);
  if (tileShape[2] < 1) {
    tileShape[2] = 1;
  }
  // Use TSM for UVW.
  // Use as many rows as used for the DATA columns, but minimal 1024.
  int tsmnrow = tileShape[2];
  if (tsmnrow < 1024) {
    tsmnrow = 1024;
  }
  TableCopy::setTiledStMan (dminfo, Vector<String>(1, "UVW"),
                            "TiledColumnStMan", "TiledUVW",
                            IPosition(2, 3, tsmnrow));
  // Replace all non-writable storage managers by SSM.
  dminfo = TableCopy::adjustStMan (dminfo);
  SetupNewTable newtab(OutName, tempdesc, Table::NewNoReplace);
  newtab.bindCreate (dminfo);
  Table outtable(newtab);
  {
    // Add DATA column using tsm.
    TiledColumnStMan tsm("TiledData", tileShape);
    TableResize(desc, data_ipos, &tsm, outtable);
  }
  {
    // Add FLAG column using tsm.
    // Use larger tile shape because flags are stored as bits.
    IPosition tileShapeF(tileShape);
    tileShapeF[2] *= 8;
    TiledColumnStMan tsmf("TiledFlag", tileShapeF);
    TableResize(tdesc["FLAG"], data_ipos, &tsmf, outtable);
  }
  {
    // Add WEIGHT_SPECTRUM column using tsm.
    TiledColumnStMan tsmw("TiledWeightSpectrum", tileShape);
    ArrayColumnDesc<Float> wsdesc("WEIGHT_SPECTRUM", "weight per pol/chan",
                                  data_ipos, ColumnDesc::FixedShape);
    TableResize(wsdesc, data_ipos, &tsmw, outtable);
  }
  // If both present handle the CORRECTED_DATA and MODEL_DATA column.
  if (Details.AllColumns)
  {
    if (tdesc.isColumn("CORRECTED_DATA") && tdesc.isColumn("MODEL_DATA"))
    {
      // Add them to the columns to handle.
      Details.DataColumns.resize (3);
      Details.DataColumns[1] = "CORRECTED_DATA";
      Details.DataColumns[2] = "MODEL_DATA";
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
      TiledColumnStMan tsmm("ModelData", tileShape);
      TableResize(mdesc, data_ipos, &tsmm, outtable);

      cout << "CORRECTED_DATA detected for processing" << endl;
      TiledColumnStMan tsmc("CorrectedData", tileShape);
      TableResize(tdesc["CORRECTED_DATA"], data_ipos, &tsmc, outtable);

      IPosition iwShape(1, data_ipos[1]);
      IPosition iwShapeTile(2, tileShape[1], tileShape[2]);
      TiledColumnStMan tsmw("TiledImagingWeight", iwShapeTile);
      ColumnDesc iwdesc(ArrayColumnDesc<float>("IMAGING_WEIGHT"));
      TableResize(iwdesc, iwShape, &tsmw, outtable);
    }
    else
    {
      if (tdesc.isColumn("CORRECTED_DATA") || tdesc.isColumn("MODEL_DATA"))
      {
        cout << "Only one of CORRECTED_DATA and MODEL_DATA columns is present; "
             << "it is ignored" << endl;
      }
    }
  }
  ASSERTSTR (Details.FlagColumn == "DATA"  ||  Details.DataColumns.size() > 1,
             "DataColumn " << Details.FlagColumn
             << " to be used in flagging does not exist");
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
  ScalarColumn<Double> outTOTALBW(outSPW, "TOTAL_BANDWIDTH");

  Vector<Double> old_freq(old_nchan, 0.0);
  Vector<Double> new_freq(new_nchan, 0.0);
  Vector<Double> old_width(old_nchan, 0.0);
  Vector<Double> new_width(new_nchan, 0.0);
  Vector<Double> old_bw(old_nchan, 0.0);
  Vector<Double> new_bw(new_nchan, 0.0);
  Vector<Double> old_res(old_nchan, 0.0);
  Vector<Double> new_res(new_nchan, 0.0);

  for (unsigned int i = 0; i < inSPW.nrow(); i++)
  {
    inFREQ.get(i, old_freq);
    inWIDTH.get(i, old_width);
    inBW.get(i, old_bw);
    inRESOLUTION.get(i, old_res);
    double totalbw = 0;
    uint first = Details.Start;
    uint last  = first + Details.Step;
    // This loops assumes regularly spaced, adjacent frequency channels.
    for (int j = 0; j < new_nchan; j++)
    { 
      if (last > old_freq.size()) {
        last = old_freq.size();
      }
      int nrch = last-first;
      new_freq[j]  = 0.5 * (old_freq[first] + old_freq[last-1]);
      new_width[j] = nrch * old_width[0];
      new_bw[j]    = nrch * old_bw[0];
      new_res[j]   = nrch * old_res[0];
      totalbw += new_bw[j];
      first = last;
      last += Details.Step;
    }
    outFREQ.put(i, new_freq);
    outWIDTH.put(i, new_width);
    outBW.put(i, new_bw);
    outRESOLUTION.put(i, new_res);
    outTOTALBW.put(i, totalbw);
  }
  // Write the parset info into the history.
  writeHistory (*OutMS, Details);
  OutMS->flush(true);
  cout << "Finished preparing output MS" << endl;
}

void MsFile::writeHistory (MeasurementSet& ms, const RunDetails& details)
{
  MSHistory histtab(ms.history());
  cout << "iswritable: " << ms.isWritable()<<' '<<histtab.isWritable()<<' '<<ms.history().isWritable()<<endl;
  histtab.reopenRW();
  cout << "iswritable: " << ms.isWritable()<<' '<<histtab.isWritable()<<' '<<ms.history().isWritable()<<endl;
  uInt rownr = histtab.nrow();
  histtab.addRow();
  MSHistoryColumns histcols(histtab);
  histcols.observationId().put (rownr, 0);
  histcols.application().put   (rownr, "DPPP");
  histcols.message().put       (rownr, details.AllParms);
  histcols.priority().put      (rownr, "NORMAL");
  histcols.time().put          (rownr, Time().modifiedJulianDay()*24.*3600.);
  histcols.origin().put        (rownr, Version::getInfo<DPPPVersion>("IDPPP",
                                                                     "full"));
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
  return TableIterator(itsOrderedTable, ms_iteration_variables,
                       TableIterator::Ascending, TableIterator::NoSort);
}

//===============>>> MsFile::UpdateTimeslotData  <<<===============
void MsFile::UpdateTimeslotData(casa::TableIterator& Data_iter,
                                MsInfo& Info,
                                DataBuffer& Buffer,
                                TimeBuffer& TimeData,
                                bool missingTime,
                                double timeValue)
{
  Table         TimeslotTable = Data_iter.table();
  int           rowcount      = TimeslotTable.nrow();
  const vector<string>& dataColumns = Buffer.dataColumns();
  ROTableVector<Int>            antenna1      (TimeslotTable, "ANTENNA1");
  ROTableVector<Int>            antenna2      (TimeslotTable, "ANTENNA2");
  ROTableVector<Int>            bandnr        (TimeslotTable, "DATA_DESC_ID");
  ROArrayColumn<Double>         uvw           (TimeslotTable, "UVW");
  ROTableVector<Double>         time_centroid (TimeslotTable, "TIME_CENTROID");
  ROTableVector<Double>         interval      (TimeslotTable, "INTERVAL");
  ROTableVector<Double>         exposure      (TimeslotTable, "EXPOSURE");
  ROArrayColumn<Bool>           flags         (TimeslotTable, "FLAG");
  ROArrayColumn<Float>          weights;
  if (itsHasWeightSpectrum) {
    weights.attach(TimeslotTable, "WEIGHT_SPECTRUM");
  } else {
    weights.attach(TimeslotTable, "WEIGHT");
  }
  vector<ROArrayColumn<Complex> > data(dataColumns.size());
  vector<Cube<Complex> >        tempData(dataColumns.size());
  for (uint i=0; i<data.size(); ++i) {
    data[i].attach(TimeslotTable, dataColumns[i]);
    tempData[i].resize (Info.NumPolarizations, Info.NumChannels, rowcount);
  }
  Cube<Bool>                    tempFlags(Info.NumPolarizations, Info.NumChannels, rowcount);
  Matrix<Float>                 tempWeights(Info.NumPolarizations, rowcount);
  Cube<Float>                   tempWeightSpectrum(Info.NumPolarizations, Info.NumChannels, rowcount);

  if (missingTime) {
    // Clear the array; note that tempData is already initialized to Complex().
    tempFlags = True;
    tempWeights = 0;
    tempWeightSpectrum = 0;
  } else {
    // Get all DATA columns
    for (uint i=0; i<data.size(); ++i) {
      data[i].getColumn(tempData[i]);
    }
    flags.getColumn(tempFlags);
    if (itsHasWeightSpectrum) {
      weights.getColumn(tempWeightSpectrum);
    } else {
      weights.getColumn(tempWeights);
    }
  }

  for (int i = 0; i < rowcount; i++)
  {
    int bi    = Info.getBaselineIndex(antenna1(i), antenna2(i));
    ASSERT (bi>=0);
    int band  = bandnr(i);
    int index = (band % Info.NumBands) * Info.NumPairs + bi;
    for (uint j=0; j<data.size(); ++j) {
      Buffer.Data[j][index].xyPlane(Buffer.Position) = tempData[j].xyPlane(i);
    }
    Buffer.Flags[index].xyPlane(Buffer.Position) = tempFlags.xyPlane(i);
    // Create reference to slice in weigth array.
    // Fill it with the input weight-spectrum if given.
    // Otherwise use the weight per pol and copy for all channels.
    Array<Float> curWeight (Buffer.Weights[index].xyPlane(Buffer.Position));
    if (itsHasWeightSpectrum) {
      curWeight = tempWeightSpectrum.xyPlane(i);
    } else {
      // Only a weight per polarization, so copy for all channels
      AlwaysAssert (curWeight.contiguousStorage(), AipsError);
      Float* dstp = curWeight.data();
      const Float* srcp = tempWeights.data() + i*Info.NumPolarizations;
      for (int j=0; j<Info.NumChannels; ++j) {
        for (int k=0; k<Info.NumPolarizations; ++k) {
          *dstp++ = srcp[k];
        }
      }
    }

    TimeData.BufTime[index].push_front(timeValue);
    TimeData.BufTimeCentroid[index].push_front(time_centroid(i));
    TimeData.BufInterval[index].push_front(interval(i));
    TimeData.BufExposure[index].push_front(exposure(i));
    TimeData.BufUvw[index].push_front(uvw(i));
  }
}

//===============>>> MsFile::WriteData  <<<===============

void MsFile::WriteData(casa::TableIterator& Data_iter,
                       MsInfo& Info,
                       DataBuffer& Buffer,
                       TimeBuffer& TimeData)
{
  uint  rowcount  = Data_iter.table().nrow();
  int   pos       = (Buffer.Position+1) % Buffer.WindowSize;
  Table DataTable = *OutMS;
  uint  strow     = DataTable.nrow();
  const vector<string>& dataColumns = Buffer.dataColumns();
  if (OutName.empty()) {
    ASSERT (! itsIterator.pastEnd());
    DataTable = itsIterator.table();
    itsIterator.next();
    ASSERT (DataTable.nrow() == rowcount);
    strow = 0;
  } else {
    Table temptable = Data_iter.table().project(SELECTblock);
    DataTable.addRow(rowcount);
    Table dummy = DataTable.project(SELECTblock);
    TableCopy::copyRows(dummy, temptable, strow, 0, rowcount, False);
  }

  ROTableVector<Int>        antenna1     (DataTable, "ANTENNA1");
  ROTableVector<Int>        antenna2     (DataTable, "ANTENNA2");
  ROTableVector<Int>        bandnr       (DataTable, "DATA_DESC_ID");
  TableVector<Double>       time         (DataTable, "TIME");
  TableVector<Double>       time_centroid(DataTable, "TIME_CENTROID");
  TableVector<Double>       exposure     (DataTable, "EXPOSURE");
  TableVector<Double>       interval     (DataTable, "INTERVAL");
  ArrayColumn  <Bool>       flags        (DataTable, "FLAG");
  vector<ArrayColumn<Complex> > data(dataColumns.size());
  ArrayColumn  <Double>     uvw;
  ArrayColumn  <Float>      weights;
  if (! OutName.empty()) {
    uvw.attach     (DataTable, "UVW");
    weights.attach (DataTable, "WEIGHT_SPECTRUM");
    for (uint i=0; i<data.size(); ++i) {
      data[i].attach (DataTable, dataColumns[i]);
    }
  }
  //cout << "Processing: " << MVTime(temp(0)/(24*3600)).string(MVTime::YMD) << endl; //for testing purposes

  for (uint i = 0; i < rowcount; i++)
  {
    int bi    = Info.getBaselineIndex(antenna1(i), antenna2(i));
    ASSERT (bi>=0);
    int band  = bandnr(i);
    int index = (band % Info.NumBands) * Info.NumPairs + bi;

    flags.put(strow + i, Buffer.Flags[index].xyPlane(pos));
    if (! OutName.empty()) {
      for (uint j=0; j<data.size(); ++j) {
        data[j].put(strow + i, Buffer.Data[j][index].xyPlane(pos));
      }
      time.set(strow + i, TimeData.Time[index][0]);
      time_centroid.set(strow + i, TimeData.TimeCentroid[index][0]);
      exposure.set(strow + i, TimeData.Exposure[index][0]);
      interval.set(strow + i, TimeData.Interval[index][0]);
      uvw.put(strow + i, TimeData.Uvw[index][0]);
      weights.put(strow + i, Buffer.Weights[index].xyPlane(pos));
    }
  }
  TimeData.Clear();
}


void MsFile::checkGaps(const MsInfo& info, bool updateMS) const
{
  Vector<Double> times = ROScalarColumn<Double>(itsOrderedTable,"TIME").getColumn();
  // Make unique; use insertion sort, because it is already in time order.
  int nrtim = genSort (times, Sort::InsSort + Sort::NoDuplicates);
  // Now check if data set is regular.
  int nrbasel = info.NumPairs;
  int nrsample = InMS->nrow();
  int nrband = nrsample / (nrbasel * nrtim);
  ASSERTSTR (nrband * nrtim * nrbasel == nrsample,
             "The MS cannot be handled by DPPP; it should contain:\n"
             " - cross and auto-correlations for all antennae in ANTENNA table\n"
             " - no missing time slots\n"
             " - the same interval length for each time slot\n"
             " #bands=" << nrband << " #timeslots=" << nrtim
             << " #baselines=" << nrbasel << " #samples=" << nrsample);
  // Now test if times are regular.
  ROScalarColumn<Double> intvCol(itsOrderedTable, "INTERVAL");
  Double intv = intvCol(0);
  bool correct = false;
  for (int i=1; i<nrtim-1; ++i) {
    Double diff = times[i] - times[i-1];
    if (!near (diff, intv, 1e-3)) {
      if (!correct) {
        cout << "Normal time interval is " << intv << " seconds" << endl;
      }
      cout << "Time slot " << i << " is " << diff
           << " seconds after previous";
      if (!updateMS) {
        cout << "; will add " << int(diff/intv + 0.1)-1
             << " flagged time slots";
      }
      cout << endl;
    }
  }
}

//===============>>> MsFile  <<<===============
