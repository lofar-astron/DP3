// VdsMaker.cc: Class to create the description of an MS
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "VdsMaker.h"
#include "VdsDesc.h"
#include "ClusterDesc.h"
#include "StreamUtil.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSMainColumns.h>
#include <casacore/ms/MeasurementSets/MSAntenna.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/ms/MeasurementSets/MSField.h>
#include <casacore/ms/MeasurementSets/MSFieldColumns.h>
#include <casacore/ms/MeasurementSets/MSPolarization.h>
#include <casacore/ms/MeasurementSets/MSPolColumns.h>
#include <casacore/ms/MeasurementSets/MSDataDescription.h>
#include <casacore/ms/MeasurementSets/MSDataDescColumns.h>
#include <casacore/ms/MeasurementSets/MSSpectralWindow.h>
#include <casacore/ms/MeasurementSets/MSSpWindowColumns.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Utilities/LinearSearch.h>
#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/OS/HostInfo.h>
#include <casacore/casa/Exceptions/Error.h>

#include <boost/make_unique.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace dp3;
using namespace casacore;
using namespace std;

namespace dp3 {
namespace common {

void VdsMaker::getFreqInfo(MS& ms, vector<int>& nrchan,
                           vector<Vector<double>>& startFreq,
                           vector<Vector<double>>& endFreq) {
  MSDataDescription mssub1(ms.dataDescription());
  ROMSDataDescColumns mssub1c(mssub1);
  MSSpectralWindow mssub(ms.spectralWindow());
  ROMSSpWindowColumns mssubc(mssub);
  int nrspw = mssub1.nrow();
  nrchan.reserve(nrspw);
  startFreq.reserve(nrspw);
  endFreq.reserve(nrspw);
  for (int spw = 0; spw < nrspw; ++spw) {
    Vector<double> chanFreq = mssubc.chanFreq()(spw);
    Vector<double> chanWidth = mssubc.chanWidth()(spw);
    nrchan.push_back(chanFreq.size());
    startFreq.push_back(chanFreq - chanWidth / 2.);
    endFreq.push_back(chanFreq + chanWidth / 2.);
  }
}

void VdsMaker::getFields(MS& ms, vector<double>& ra, vector<double>& dec,
                         vector<string>& refType) {
  MSField mssub(ms.field());
  ROMSFieldColumns mssubc(mssub);
  int nrf = mssub.nrow();
  ra.resize(nrf);
  dec.resize(nrf);
  refType.resize(nrf);
  for (int i = 0; i < nrf; ++i) {
    Array<MDirection> mds = mssubc.referenceDirMeasCol()(i);
    ra[i] = mds.data()[i].getValue().get()[0];
    dec[i] = mds.data()[i].getValue().get()[1];
    refType[i] = mds.data()[i].getRefString();
  }
}

void VdsMaker::getAntNames(MS& ms, vector<string>& antNames) {
  MSAntenna mssub(ms.antenna());
  ROMSAntennaColumns mssubc(mssub);
  Vector<String> names = mssubc.name().getColumn();
  antNames.resize(names.size());
  for (unsigned int i = 0; i < names.size(); ++i) {
    antNames[i] = names[i];
  }
}

void VdsMaker::getCorrInfo(MS& ms, vector<string>& corrTypes) {
  MSPolarization mssub(ms.polarization());
  if (mssub.nrow() > 0) {
    ROMSPolarizationColumns mssubc(mssub);
    Vector<Int> ctvec = mssubc.corrType()(0);
    int nrp = ctvec.nelements();
    corrTypes.resize(nrp);
    for (int i = 0; i < nrp; ++i) {
      corrTypes[i] = Stokes::name(Stokes::type(ctvec(i)));
    }
  }
}

void VdsMaker::getDataFileInfo(MS& ms, string& name, bool& regular,
                               vector<int>& tileShape, vector<int>& cubeShape) {
  regular = false;
  Record rec = ms.dataManagerInfo();
  // Find the subrecord containing the DATA column.
  for (unsigned int i = 0; i < rec.nfields(); ++i) {
    const Record& subrec = rec.subRecord(i);
    Vector<String> colNames(subrec.asArrayString("COLUMNS"));
    int inx = linearSearch1(colNames, String("DATA"));
    if (inx >= 0) {
      const Record& specrec = subrec.subRecord("SPEC");
      // Stored as a hypercube?
      if (specrec.isDefined("HYPERCUBES")) {
        // Non-regular if multiple columns in storage manager.
        regular = (colNames.size() == 1);
        const Record& hcrec = specrec.subRecord("HYPERCUBES");
        if (hcrec.nfields() != 1) {
          regular = false;
        } else {
          const Record& tsmrec = hcrec.subRecord(0);
          tsmrec.asArrayInt("TileShape").tovector(tileShape);
          tsmrec.asArrayInt("CubeShape").tovector(cubeShape);
          // Find the name of the file containing the data.
          // It can be TSM1 or TSM0.
          ostringstream str;
          str << specrec.asInt("SEQNR");
          name = string(ms.tableName()) + "/table.f" + str.str() + "_TSM1";
          if (!File(name).exists()) {
            name = string(ms.tableName()) + "/table.f" + str.str() + "_TSM0";
            if (!File(name).exists()) {
              regular = false;
            }
          }
        }
      }
      break;
    }
  }
}

string VdsMaker::findFileSys(const string& fileName,
                             const dp3::common::ClusterDesc& cdesc,
                             const string& hostName) {
  // Find the file system by looking for a matching mountpoint.
  const vector<dp3::common::NodeDesc>& nodes = cdesc.getNodes();
  // First find the NodeDesc for this node.
  string nodeName(hostName);
  if (nodeName.empty()) {
    nodeName = "localhost";
  }
  unsigned int i = 0;
  // If no hostname is given, try localhost and the real hostname.
  for (int j = 0; j < 2; ++j) {
    i = 0;
    for (; i < nodes.size(); ++i) {
      if (nodes[i].getName() == nodeName) {
        break;
      }
    }
    if (i < nodes.size() || !hostName.empty()) {
      break;
    }
    nodeName = HostInfo::hostName();
  }
  if (i >= nodes.size())
    throw std::runtime_error("Hostname '" + nodeName +
                             "' not found in "
                             "ClusterDesc file ");
  return nodes[i].findFileSys(fileName);
}

void VdsMaker::create(const string& msName, const string& outName,
                      const string& clusterDescName, const string& hostName,
                      bool fillTimes) {
  // Open the table.
  MS ms(msName);
  // Create and fill the Vds object.
  dp3::common::VdsPartDesc msd;
  ostringstream oss;
  // Fill in MS path and name.
  Path mspr(msName);
  string absName = mspr.absoluteName();
  // If the ClusterDesc file is given, try to find filesys and put its
  // absolute path into the VdsPartDesc.
  // Otherwise it is unknown.
  if (clusterDescName.empty()) {
    msd.setName(absName, "unknown");
  } else {
    Path cdpath(clusterDescName);
    msd.setClusterDescName(cdpath.absoluteName());
    dp3::common::ClusterDesc cdesc(clusterDescName);
    msd.setName(absName, findFileSys(absName, cdesc, hostName));
  }
  msd.setFileName(absName);
  // Get freq info.
  // Fill in correlation info.
  vector<string> corrNames;
  getCorrInfo(ms, corrNames);
  ostringstream oss1;
  oss1 << corrNames;
  msd.addParm("CorrNames", oss1.str());
  // Fill in freq info.
  vector<int> nchan;
  vector<Vector<double>> startFreq, endFreq;
  getFreqInfo(ms, nchan, startFreq, endFreq);
  for (unsigned int i = 0; i < nchan.size(); ++i) {
    vector<double> sfreq, efreq;
    startFreq[i].tovector(sfreq);
    endFreq[i].tovector(efreq);
    msd.addBand(nchan[i], sfreq, efreq);
  }
  // Write the field directions (in J2000).
  vector<double> ra, dec;
  vector<string> refType;
  getFields(ms, ra, dec, refType);
  int nrfield = ra.size();
  ostringstream oss2a, oss2b, oss2c;
  oss2a << '[';
  oss2b << '[';
  oss2c << '[';
  for (int i = 0; i < nrfield; ++i) {
    if (i > 0) {
      oss2a << ',';
      oss2b << ',';
      oss2c << ',';
    }
    oss2a << MVAngle::Format(MVAngle::TIME, 12)
          << MVAngle(Quantity(ra[i], "rad"));
    oss2b << MVAngle::Format(MVAngle::ANGLE, 12)
          << MVAngle(Quantity(dec[i], "rad"));
    oss2c << refType[i];
  }
  oss2a << ']';
  oss2b << ']';
  oss2c << ']';
  msd.addParm("FieldDirectionRa", oss2a.str());
  msd.addParm("FieldDirectionDec", oss2b.str());
  msd.addParm("FieldDirectionType", oss2c.str());
  // Fill in station names.
  vector<string> antNames;
  getAntNames(ms, antNames);
  ostringstream oss2;
  oss2 << antNames;
  msd.addParm("StationNames", oss2.str());
  // Get the data file name.
  string dfName;
  bool dfRegular;
  vector<int> tileShape, cubeShape;
  getDataFileInfo(ms, dfName, dfRegular, tileShape, cubeShape);
  msd.addParm("DataFileName", dfName);
  ostringstream oss3;
  oss3 << dfRegular;
  msd.addParm("DataFileIsRegular", oss3.str());
  ostringstream oss4;
  oss4 << tileShape;
  msd.addParm("DataTileShape", oss4.str());
  std::ostringstream oss5;
  oss5 << cubeShape;
  msd.addParm("DataCubeShape", oss5.str());

  // Fill in times.
  ROMSMainColumns mscol(ms);
  uInt nrow = ms.nrow();
  if (nrow <= 0) {
    std::cout << "MeasurementSet " + absName + " is empty";
  } else {
    // Get start and end time. Get the step time from the middle one.
    double stepTime = mscol.exposure()(nrow / 2);
    double startTime = mscol.time()(0) - mscol.exposure()(0) / 2;
    double endTime = mscol.time()(nrow - 1) + mscol.exposure()(nrow - 1) / 2;
    if (fillTimes) {
      // Get all unique times.
      Table msuniq = ms.sort("TIME", Sort::Ascending,
                             Sort::QuickSort + Sort::NoDuplicates);
      Vector<double> tims = ROScalarColumn<double>(msuniq, "TIME").getColumn();
      Vector<double> intv =
          ROScalarColumn<double>(msuniq, "INTERVAL").getColumn();
      vector<double> stimes(tims.size());
      vector<double> etimes(tims.size());
      for (unsigned int i = 0; i < tims.size(); ++i) {
        stimes[i] = tims[i] - intv[i] * 0.5;
        etimes[i] = tims[i] + intv[i] * 0.5;
      }
      msd.setTimes(startTime, endTime, stepTime, stimes, etimes);
    } else {
      msd.setTimes(startTime, endTime, stepTime);
    }
  }
  // Write into the vds file.
  ofstream ostr(outName.c_str());
  if (!ostr)
    throw std::runtime_error("File " + outName + " could not be created");
  msd.write(ostr, "");
}

void VdsMaker::combine(const string& gdsName, const vector<string>& vdsNames) {
  // Form the global desc.
  dp3::common::VdsPartDesc globalvpd;
  globalvpd.setName(gdsName, string());
  vector<double> sfreq(1);
  vector<double> efreq(1);

  // Read all parts.
  // Add a band, but with only its start and end freq (not all freqs).
  // Determine the minimum and maximum time.
  double startTime = 1e30;
  double endTime = 0;
  std::vector<std::unique_ptr<dp3::common::VdsPartDesc>> vpds;
  vpds.reserve(vdsNames.size());
  for (unsigned int j = 0; j < vdsNames.size(); ++j) {
    auto vpd =
        boost::make_unique<dp3::common::VdsPartDesc>(ParameterSet(vdsNames[j]));
    // Skip a VDS with an empty time (it has no data).
    casacore::Path path(vdsNames[j]);
    // File name gets the original MS name.
    // Name gets the name of the VDS file.
    vpd->setFileName(vpd->getName());
    vpd->setName(path.absoluteName(), vpd->getFileSys());
    const vector<int>& chans = vpd->getNChan();
    const vector<double>& sf = vpd->getStartFreqs();
    const vector<double>& ef = vpd->getEndFreqs();
    int inxf = 0;
    for (unsigned int i = 0; i < chans.size(); ++i) {
      int nchan = chans[i];
      sfreq[0] = sf[inxf];
      // A band can be given with individual freqs or a single freq range.
      if (chans.size() == sf.size()) {
        ++inxf;
      } else {
        inxf += nchan;
      }
      efreq[0] = ef[inxf - 1];
      globalvpd.addBand(nchan, sfreq, efreq);
    }
    // Get minimum/maximum time.
    if (vpd->getStartTime() == 0) {
      std::cout << "Dataset " << vdsNames[j] << " is completely empty";
    } else {
      startTime = std::min(startTime, vpd->getStartTime());
      endTime = std::max(endTime, vpd->getEndTime());
    }
    vpds.push_back(std::move(vpd));
  }
  // Exit if no valid VDS files.
  if (vpds.empty()) throw std::runtime_error("No VDS files are given");
  if (startTime == 0)
    throw std::runtime_error("All datasets seems to be empty");

  // Set the times in the global desc (using the first part).
  // Set the clusterdesc name.
  // If defined, set the Extra parameters giving the field directions.
  // Form the global desc.
  globalvpd.setTimes(startTime, endTime, vpds[0]->getStepTime(),
                     vpds[0]->getStartTimes(), vpds[0]->getEndTimes());
  globalvpd.setClusterDescName(vpds[0]->getClusterDescName());
  if (vpds[0]->getParms().isDefined("FieldDirectionRa")) {
    globalvpd.addParm("FieldDirectionRa",
                      vpds[0]->getParms().getString("FieldDirectionRa"));
    globalvpd.addParm("FieldDirectionDec",
                      vpds[0]->getParms().getString("FieldDirectionDec"));
    globalvpd.addParm("FieldDirectionType", vpds[0]->getParms().getString(
                                                "FieldDirectionType", "J2000"));
  }
  dp3::common::VdsDesc gdesc(globalvpd);

  // Now add all parts to the global desc and write it.
  // Print a warning if times differ.
  // Also cleanup the objects.
  for (unsigned int i = 0; i < vpds.size(); ++i) {
    vpds[i]->clearParms();
    gdesc.addPart(*vpds[i]);
    if (vpds[i]->getStartTime() != globalvpd.getStartTime() ||
        vpds[i]->getEndTime() != globalvpd.getEndTime() ||
        vpds[i]->getStepTime() != globalvpd.getStepTime()) {
      cerr << "The times of part " << i << " (" << vpds[i]->getName()
           << ") are different" << endl;
    }
    vpds[i].reset();
  }
  ofstream ostr(gdsName.c_str());
  if (!ostr)
    throw std::runtime_error("File " + gdsName + " could not be created");
  gdesc.write(ostr);
}

}  // namespace common
}  // namespace dp3
