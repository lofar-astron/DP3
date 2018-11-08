#include "H5Parm.h"

#include "Exceptions.h"

#include "../Common/StringUtil.h"

#include <cstring>
#include <complex>
#include <boost/filesystem/operations.hpp>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>

#include <hdf5.h>

using namespace std;

namespace DP3 {
  H5Parm::H5Parm(const std::string& filename, bool forceNew,
                 bool forceNewSolSet, const std::string& solSetName) :
      H5::H5File(filename, forceNew?H5F_ACC_TRUNC:H5F_ACC_RDONLY)
  {
    if (forceNewSolSet || getNumObjs()==0) { // Create a new solSet
      if (solSetName=="") {
        // Get the name of first non-existing solset
        stringstream newSolSetName;
        H5::Group tryGroup;
        for (uint solSetNum=0; solSetNum<100; ++solSetNum) {
          try {
            H5::Exception::dontPrint();
            newSolSetName<<"sol"<<setfill('0')<<setw(3)<<solSetNum;
            tryGroup = openGroup(newSolSetName.str());
            newSolSetName.str("");
          }
          catch (H5::FileIException& not_found_error ) {
            // solSetName does not exist yet
            break;
          }
          tryGroup.close();
        }
        _solSet = createGroup("/"+newSolSetName.str(), H5P_DEFAULT);
      } else {
        // Create solset with the given name
        _solSet = createGroup("/"+solSetName, H5P_DEFAULT);
      }
      addVersionStamp(_solSet);
    } else {
      string solSetNameToOpen=solSetName;
      if (solSetNameToOpen=="") {
        if (this->getNumObjs()==1) {
          solSetNameToOpen=this->getObjnameByIdx(0);
        } else {
          throw Exception("H5Parm " + filename + " contains more than one SolSet, " +
              "please specify which one to use.");
        }
      }

      _solSet = openGroup(solSetNameToOpen);

      vector<string> solTabNames;
      for (uint i=0; i<_solSet.getNumObjs();++i) {
        if (_solSet.getObjTypeByIdx(i)==H5G_GROUP) {
          solTabNames.push_back(_solSet.getObjnameByIdx(i));
        }
      }

      for (vector<string>::iterator solTabName=solTabNames.begin();
           solTabName!=solTabNames.end(); ++solTabName) {
        H5::Group group = _solSet.openGroup(*solTabName);
        _solTabs.insert(
            std::map<std::string, SolTab>::value_type (*solTabName, SolTab(group)));
      }
    }
  }


  H5Parm::H5Parm() {
  }

  H5Parm::~H5Parm() {
    // Throw an error if the antenna or source table is not present
    //_solSet.openDataSet("antenna");
    //_solSet.openDataSet("source");
    _solSet.close();
  }

  string H5Parm::getSolSetName() const {
    char buffer[100];
    hsize_t namelen = H5Iget_name(_solSet.getId(),buffer,100);
    buffer[namelen+1]=0;
    // Strip leading '/'
    return buffer+1;
  }

  void H5Parm::addVersionStamp(H5::Group &node) {
    // Write an attribute with the h5parm version
    H5::Attribute attr = node.createAttribute("h5parm_version",
                                              H5::StrType(H5::PredType::C_S1, 3),
                                              H5::DataSpace());
    attr.write(H5::StrType(H5::PredType::C_S1, 3), "1.0");
  }

  void H5Parm::addSources (const std::vector<std::string>& names,
                           const std::vector<std::pair<double, double> >& dirs) {
    hsize_t dims[1];

    // Create data type
    dims[0]=2;  // For ra, dec in directions
    H5::CompType sourceType(sizeof(source_t));
    sourceType.insertMember("name", HOFFSET(antenna_t, name), H5::StrType(H5::PredType::C_S1, 128));
    sourceType.insertMember("dir", HOFFSET(source_t, dir), H5::ArrayType(H5::PredType::NATIVE_FLOAT, 1, dims));

    // Create dataset
    dims[0] = names.size();
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = _solSet.createDataSet("source", sourceType, dataspace);

    // Prepare data
    vector<source_t> sources(names.size());
    for (uint src=0; src<sources.size(); ++src) {
      std::strncpy(sources[src].name, names[src].c_str(), 128);
      sources[src].dir[0] = dirs[src].first;
      sources[src].dir[1] = dirs[src].second;
    }

    // Write data
    dataset.write(&(sources[0]), sourceType);
  }

  void H5Parm::addAntennas (const std::vector<std::string>& names,
                            const std::vector<std::vector<double> >& positions) {
    hsize_t dims[1];

    // Create data type
    dims[0]=3;  // For x,y,z in positions
    H5::CompType antennaType(sizeof(antenna_t));
    antennaType.insertMember("name", HOFFSET(antenna_t, name), H5::StrType(H5::PredType::C_S1, 16));
    antennaType.insertMember("position", HOFFSET(antenna_t, position), H5::ArrayType(H5::PredType::NATIVE_FLOAT, 1, dims));

    // Create dataset
    dims[0] = names.size();
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = _solSet.createDataSet("antenna", antennaType, dataspace);

    // Prepare data
    vector<antenna_t> ants(names.size());
    for (uint ant=0; ant<ants.size(); ++ant) {
      std::strncpy(ants[ant].name, names[ant].c_str(), 16);
      const std::vector<double>& pos = positions[ant];
      ants[ant].position[0] = pos[0];
      ants[ant].position[1] = pos[1];
      ants[ant].position[2] = pos[2];
    }

    dataset.write(&(ants[0]), antennaType);
  }

  H5Parm::SolTab& H5Parm::getSolTab(const std::string& name) {
    std::map<std::string, SolTab>::iterator item =
      _solTabs.find(name);
    if (item == _solTabs.end()) {
      throw Exception("SolTab " + name + " does not exist in solset " +
                       getSolSetName());
    }
    return item->second;
  }

  bool H5Parm::hasSolTab(const string& solTabName) const {
    return _solTabs.find(solTabName) != _solTabs.end();
  }

  H5Parm::SolTab& H5Parm::createSolTab(const std::string& name,
                                       const std::string& type,
                                       const std::vector<H5Parm::AxisInfo> axes) {
    H5::Group newgroup = _solSet.createGroup(name);
    std::map<std::string, SolTab>::iterator newItem =
      _solTabs.insert(std::make_pair(name, SolTab(newgroup, type, axes))).first;
    return newItem->second;
  }
}
