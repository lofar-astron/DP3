/***************************************************************************
 *   Copyright (C) 2006-8 by Adriaan Renting and ASTRON/LOFAR              *
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
#ifndef LOFARPIPELINEPROCESSCONTROL_H
#define LOFARPIPELINEPROCESSCONTROL_H

#include <PLC/ProcessControl.h>
#include <Common/ParameterSet.h>

/**
@author Adriaan Renting
*/
/// @file
/// @brief Class to hold code for inititalising IDPP and starting the processing in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    ///Foreward declarations
    class BandpassCorrector;
    class Flagger;
    class DataSquasher;
    class MsFile;
    class MsInfo;
    class RunDetails;
    class Pipeline;

    class PipelineProcessControl : public LOFAR::ACC::PLC::ProcessControl
    {
    private:
      std::string  itsInMS; ///< Name and location of input MS
      std::string  itsOutMS; ///< Name and location of output MS
      std::string  itsClusterDesc; ///< Name and location of ClusterDesc file
      std::string  itsVdsDir; ///< Optional directory where to store VDS files
      unsigned int itsBandpass; ///< from ParameterSet
      unsigned int itsFlagger; ///< From ParameterSet
      unsigned int itsSquasher; ///< From ParameterSet

      Pipeline*          myPipeline; ///< the actual pipeline, all others need to be initialised first
      MsFile*            myFile;
      MsInfo*            myInfo;
      BandpassCorrector* myBandpass;
      Flagger*           myFlagger;
      DataSquasher*      mySquasher;
      RunDetails*        myDetails;
    public:
      PipelineProcessControl();

      ~PipelineProcessControl(void);
      /// \name Command to control the processes.
      /// There are a dozen commands that can be sent to a application process
      /// to control its flow. The return values for these command are:<br>
      /// - True   - Command executed succesfully.
      /// - False  - Command could not be executed.
      ///
      /// @{

      /// During the \c define state the process check the contents of the
      /// ParameterSet it received during start-up. When everthing seems ok the
      /// process constructs the communication channels for exchanging data
      /// with the other processes. The connection are NOT made in the stage.
      tribool define   (void);

      /// When a process receives an \c init command it allocates the buffers it
      /// needs an makes the connections with the other processes. When the
      /// process succeeds in this it is ready for dataprocessing (or whatever
      /// task the process has).
      tribool init     (void);

      /// During the \c run phase the process does the work it is designed for.
      /// The run phase stays active until another command is send.
      tribool run      (void);

      tribool pause(const std::string&);
      tribool release(void);
      tribool quit(void);
      tribool recover(const std::string&);
      tribool reinit(const  std::string&);
      tribool snapshot(const std::string&);
      std::string askInfo(const std::string&);

    }; //class PipelineProcessControl
  } //namespace CS1
}; //namespace LOFAR

#endif
