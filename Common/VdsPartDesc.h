//# VdsPartDesc.h: Description of a visibility data set or part thereof
//#
//# Copyright (C) 2005
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
//# $Id: VdsPartDesc.h 16886 2010-12-08 10:43:17Z diepen $

#ifndef LOFAR_LMWCOMMON_VDSPARTDESC_H
#define LOFAR_LMWCOMMON_VDSPARTDESC_H

// @file
// @brief Description of a visibility data set or part thereof.
// @author Ger van Diepen (diepen AT astron nl)

//# Includes
#include "ParameterSet.h"

#include "../Blob/BlobOStream.h"
#include "../Blob/BlobIStream.h"

#include <string>
#include <vector>
#include <iosfwd>

namespace DP3 { namespace CEP {

  // @ingroup LMWCommon
  // @brief Description of a visibility data set or part thereof.

  // This class holds the description of a visibility data set (VDS) part.
  // It defines the name of the part and on which file system it is located.
  // Using the ClusterDesc object it can be derived on which node this
  // VDS part can be processed best. This is done by the WorkersDesc
  // class.
  //
  // The description of the VDS also contains info about the time,
  // frequency, and baseline domain of the visibility data.
  //
  // The information is made persistent in a LOFAR .parset file.

  class VdsPartDesc
  {
  public:
    // Construct an empty object.
    VdsPartDesc()
      : itsStartTime(0), itsEndTime(1), itsStepTime(1)
      {}

    // Construct from the given parameterset.
    explicit VdsPartDesc (const ParameterSet&);

    // Set VDS name and file system.
    void setName (const std::string& name, const std::string& fileSys);

    // Set the original file name.
    void setFileName (const std::string& name)
      { itsFileName = name; }

    // Set the name of the ClusterDesc file used.
    void setClusterDescName (const std::string& cdName)
      { itsCDescName = cdName; }

    // Change the base part of the name.
    void changeBaseName (const std::string& newBaseName);

    // Set the observation start and end time.
    // Optionally the start and end per time interval can be set.
    void setTimes (double startTime, double endTime, double stepTime,
                   const std::vector<double>& starts = std::vector<double>(),
                   const std::vector<double>& ends   = std::vector<double>());

    // Add a band.
    // <group>
    void addBand (int nchan, double startFreq, double endFreq);
    void addBand (int nchan, const std::vector<double>& startFreq,
		  const std::vector<double>& endFreq);
    // </group>

    // Add an extra parameter. It is added to the subset 'Extra.'.
    // If the parameter already exists, it is replaced.
    void addParm (const std::string& key, const std::string& value)
      { return itsParms.add (key, value); }

    // Get access to the extra parameters.
    const ParameterSet& getParms() const
      { return itsParms; }

    // Clear the extra parameters.
    void clearParms()
      { itsParms.clear(); }

    // Write the VdsPartDesc object in parset format.
    void write (std::ostream& os, const std::string& prefix) const;

    // Get the values.
    // @{
    const std::string& getName() const
      { return itsName; }
    const std::string& getFileName() const
      { return itsFileName; }
    const std::string& getFileSys() const
      { return itsFileSys; }
    const std::string& getClusterDescName() const
      { return itsCDescName; }
    double getStartTime() const
      { return itsStartTime; }
    double getEndTime() const
      { return itsEndTime; }
    double getStepTime() const
      { return itsStepTime; }
    const std::vector<double>& getStartTimes() const
      { return itsStartTimes; }
    const std::vector<double>& getEndTimes() const
      { return itsEndTimes; }
    int getNBand() const
      { return itsNChan.size(); }
    const std::vector<int>& getNChan() const
      { return itsNChan; }
    const std::vector<double>& getStartFreqs() const
      { return itsStartFreqs; }
    const std::vector<double>& getEndFreqs() const
      { return itsEndFreqs; }
    // @}

  // Put/get the object to/from a blob.
  // <group>
  BlobOStream& toBlob (BlobOStream&) const;
  BlobIStream& fromBlob (BlobIStream&);
  // </group>

  private:
    std::string itsName;       //# full name of the VDS desc
    std::string itsFileName;   //# full name of the VDS (data set name)
    std::string itsFileSys;    //# name of file system the VDS resides on
    std::string itsCDescName;  //# name of ClusterDesc file used
    double      itsStartTime;
    double      itsEndTime;
    double      itsStepTime;
    std::vector<double> itsStartTimes;
    std::vector<double> itsEndTimes;
    std::vector<int32_t>  itsNChan;        //# nr of channels per band
    std::vector<double> itsStartFreqs;   //# start freq of each channel
    std::vector<double> itsEndFreqs;     //# end freq of each channel
    ParameterSet itsParms;     //# extra parameters
  };
    
  // Put/get the object to/from a blob.
  // <group>
    inline BlobOStream& operator<< (BlobOStream& bs, const VdsPartDesc& vpd)
    { return vpd.toBlob (bs); }
    inline BlobIStream& operator>> (BlobIStream& bs, VdsPartDesc& vpd)
    { return vpd.fromBlob (bs); }
  // </group>

}} //# end namespaces

#endif
