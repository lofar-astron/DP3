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
#include <iostream>

#include <DPPP/Pipeline.h>
#include <DPPP/MsInfo.h>
#include <DPPP/MsFile.h>
#include <DPPP/RunDetails.h>
#include <DPPP/BandpassCorrector.h>
#include <DPPP/Flagger.h>
#include <DPPP/DataSquasher.h>
#include <DPPP/DataBuffer.h>
#include <DPPP/TimeBuffer.h>
#include <DPPP/FlaggerStatistics.h>

#include <DPPP/ProgressMeter.h>

using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  Pipeline::Pipeline  <<<===============

Pipeline::Pipeline(MsInfo* info, MsFile* msfile, RunDetails* details,
                   BandpassCorrector* bandpass, Flagger* flagger, DataSquasher* squasher):
  myInfo(info),
  myFile(msfile),
  myDetails(details),
  myBandpass(bandpass),
  myFlagger(flagger),
  mySquasher(squasher),
  BandpassData(NULL),
  FlaggerData(NULL),
  SquasherData(NULL),
  myStatistics(NULL),
  TimeData(NULL)
{
  myStatistics = new FlaggerStatistics(*myInfo);
}

//===============>>>  Pipeline::~Pipeline  <<<===============

Pipeline::~Pipeline()
{
  delete TimeData;
  delete myStatistics;
  delete BandpassData;
  if (FlaggerData != BandpassData)
  { delete FlaggerData;
  }
  if (SquasherData != FlaggerData)
  { delete SquasherData;
  }
}

//===============>>>  Pipeline::~Pipeline  <<<===============

void Pipeline::MirrorBuffer(DataBuffer& buffer, MsInfo& info, int step)
{
  int from_pos;
  int to_pos;
  if (step) //end of the read sequence
  { from_pos = (buffer.WindowSize + buffer.Position - 2*step) % buffer.WindowSize;
    to_pos   = buffer.Position;
  }
  else //start of the read sequence
  { from_pos = buffer.Position;
    to_pos   = (buffer.WindowSize - buffer.Position) % buffer.WindowSize;
  }
  for (int i = 0; i < info.NumBands * info.NumPairs; i++)
  { buffer.Data[i].xyPlane(to_pos) = buffer.Data[i].xyPlane(from_pos);
  }
}

//===============>>> ComplexMedianFlagger::UpdateTimeslotData  <<<===============
void Pipeline::Run(MsInfo* SquashedInfo, bool Columns)
{
  BandpassData = new DataBuffer(myInfo, myDetails->TimeWindow, Columns);
  // Not needed unless Flagger starts altering data, or Bandpass starts altering flags
  //  if (myFlagger && myBandpass)
  //  { FlaggerData = new DataBuffer(info, myDetails->TimeWindow);
  //  }
  //  else
  //  { FlaggerData = BandpassData;
  //  }
  FlaggerData = BandpassData;
  TimeData = new TimeBuffer(FlaggerData->NumSlots);
  TimeData->Clear();
  if (mySquasher)
  { SquasherData = new DataBuffer(SquashedInfo, 2, Columns); //two buffers, one for unflagged one for all data
  }
  else
  { SquasherData = FlaggerData;
  }

  cout << "Processing " << myFile->nrow() << " rows from input MS ..." << endl;
  ProgressMeter progress(0.0, myFile->nrow(), "DPPP", "Rows processed", "", "", True, 1);
  Double nrowProc = 0;
  progress.update (nrowProc, True);

  TableIterator time_iter    = (*myFile).TimeIterator();
  int           TimeCounter  = 0;
  int           WriteCounter = 0;
  double lastTime = 0;
  double interval = -1;
  while (!time_iter.pastEnd())
  {
    double time = ROScalarColumn<Double>(time_iter.table(), "TIME")(0);
    // If first time, fill in interval size.
    // Otherwise check if there are missing time slots.
    bool missingTime = false;
    if (interval < 0) {
      interval = ROScalarColumn<Double>(time_iter.table(), "INTERVAL")(0);
    } else {
      if (lastTime + interval + 0.1 < time) {
        time = lastTime + interval;
        missingTime = true;
      }
    }
    lastTime = time;
    BandpassData->Position = ++(BandpassData->Position) % BandpassData->WindowSize;
    myFile->UpdateTimeslotData(time_iter, *myInfo, *BandpassData, *TimeData,
                               missingTime, time);
    if (myBandpass)
    { myBandpass->ProcessTimeslot(*BandpassData, *myInfo, *myDetails);
    }
    if (TimeCounter > 0 && TimeCounter <= (BandpassData->WindowSize - 1)/2)
    { MirrorBuffer(*BandpassData, *myInfo, 0);
    }
    if (TimeCounter >= (BandpassData->WindowSize - 1)/2)
    { if (myFlagger)
      { myFlagger->ProcessTimeslot(*FlaggerData, *myInfo, *myDetails, *myStatistics);
      }
    }
    if (TimeCounter >= (FlaggerData->WindowSize - 1))
    { TimeData->ShiftBuffer();
      if (mySquasher)
      { mySquasher->ProcessTimeslot(*FlaggerData, *SquasherData, *myInfo, *myDetails, *TimeData);
      }
      if ((TimeCounter+1) % myDetails->TimeStep == 0)
      { TimeData->Squash();
        myFile->WriteData(time_iter, *SquashedInfo, *SquasherData, *TimeData);
        TimeData->Clear();
        WriteCounter++;
      }
    }
    nrowProc += time_iter.table().nrow();
    progress.update(nrowProc, True);
    TimeCounter++;
    if (!missingTime) {
      time_iter++;
    }
  }
  time_iter.reset(); //we still need a table as a template for writing the data
  for (int i = 1; i <= (FlaggerData->WindowSize - 1); i++) //write the last couple of values
  {
    FlaggerData->Position = ++(FlaggerData->Position) % FlaggerData->WindowSize;
    if (myFlagger && (i <= (FlaggerData->WindowSize - 1)/2))
    { MirrorBuffer(*FlaggerData, *myInfo, i);
      myFlagger->ProcessTimeslot(*FlaggerData, *myInfo, *myDetails, *myStatistics);
    }
    TimeData->ShiftBuffer();
    if (mySquasher)
    { mySquasher->ProcessTimeslot(*FlaggerData, *SquasherData, *myInfo, *myDetails, *TimeData);
    }
    if ((TimeCounter+1) % myDetails->TimeStep == 0)
    { TimeData->Squash();
      myFile->WriteData(time_iter, *SquashedInfo, *SquasherData, *TimeData);
      TimeData->Clear();
      WriteCounter++;
    }
    TimeCounter++;
  }
  cout << "Written timeslots: " << WriteCounter << endl;
  cout << "Processed timeslots: " << TimeCounter - FlaggerData->WindowSize + 1 << endl;
  myStatistics->PrintStatistics(std::cout);
}

//===============>>> Pipeline  <<<===============

