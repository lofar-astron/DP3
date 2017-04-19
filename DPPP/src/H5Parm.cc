#include <lofar_config.h>
#include <DPPP/H5Parm.h>
#include <Common/Exception.h>
#include <cstring>
#include <iostream>

namespace LOFAR {
  H5Parm::H5Parm(const std::string& filename) {
    // ACC_EXCL to give error when file already exists.
    try {
//      H5::Exception::dontPrint();
      _hdf5file = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException err) {
      throw Exception ("H5Parm: not able to open file "+filename+" for writing, does it exist already?");
    }

    // Make a new group
    _solSet = _hdf5file.createGroup("/sol000", H5P_DEFAULT);

    addVersionStamp(_solSet);
  }

  H5Parm::~H5Parm() {
    _solSet.close();
    _hdf5file.close();
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
    source_t tmp;
    H5::CompType sourceType(sizeof tmp);
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
    antenna_t tmp;
    H5::CompType antennaType(sizeof tmp);
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
      std::vector<double> pos = positions[ant];
      ants[ant].position[0]=pos[0];
      ants[ant].position[1]=pos[1];
      ants[ant].position[2]=pos[2];
    }

    dataset.write(&(ants[0]), antennaType);
  }

  void H5Parm::addSolution (const std::string& solName,
                            const std::string& solType,
                            const std::string& axesstr,
                            const std::vector<hsize_t>& dims,
                            const std::vector<double>& vals,
                            const std::vector<double>& weights) {
    H5::Group solTab = _solSet.createGroup(solName);
    H5::Attribute attr = solTab.createAttribute("TITLE",
                                H5::StrType(H5::PredType::C_S1, solType.size()),
                                H5::DataSpace());
    attr.write(H5::StrType(H5::PredType::C_S1, solType.size()), solType);
    addVersionStamp(solTab);

    // ASSERT dims product == vals.size()
    // ASSERT dims.size() - 1 == aantal komma's in axesstr

    H5::DataSpace dataspace(dims.size(), &(dims[0]), NULL);
    H5::DataSet dataset = solTab.createDataSet("val", 
                                           H5::PredType::IEEE_F64LE, dataspace);

    dataset.write(&(vals[0]), H5::PredType::IEEE_F64LE);

    // Write an attribute with the axes
    attr = dataset.createAttribute("AXES",
                             H5::StrType(H5::PredType::C_S1, axesstr.size()),
                             H5::DataSpace());
    attr.write(H5::StrType(H5::PredType::C_S1, axesstr.size()), axesstr);

    // Add weights
    H5::DataSet weightset = solTab.createDataSet("weight", H5::PredType::IEEE_F64LE, dataspace);

    // If weights are empty, write ones everywhere
    if (weights.empty()) {
      vector<double> fullweights(vals.size(), 1);
      weightset.write(&(fullweights[0]), H5::PredType::IEEE_F64LE);
    } else {
      weightset.write(&(weights[0]), H5::PredType::IEEE_F64LE);
    }

    attr = weightset.createAttribute("AXES",
                             H5::StrType(H5::PredType::C_S1, axesstr.size()),
                             H5::DataSpace());
    attr.write(H5::StrType(H5::PredType::C_S1, axesstr.size()), axesstr);
  }

  void H5Parm::addSolution (const std::string& solName,
                            const std::string& solType,
                            const std::string& axesstr,
                            const std::vector<hsize_t>& dims,
                            const std::vector<std::complex<double> >& vals,
                            const std::vector<double>& weights,
                            bool toAmplitudes) {
    // Convert values to real numbers by taking amplitude or argument
    vector<double> realvals(vals.size());

    if (toAmplitudes) {
      transform(vals.begin(), vals.end(), realvals.begin(), takeAbs);
    } else { // Phase only
      transform(vals.begin(), vals.end(), realvals.begin(), takeArg);
    }

    addSolution(solName, solType, axesstr, dims, realvals, weights);
  }

  void H5Parm::setSolAntennas(const std::string& solName,
                              const std::vector<std::string>& solAntennas) {
    H5::Group solTab(_solSet.openGroup(solName));

    // TODO: assert that antenna is present in antenna table in solution set
    hsize_t dims[1];
    dims[0]=solAntennas.size();

    // Create dataset
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = solTab.createDataSet("ant", H5::StrType(H5::PredType::C_S1, 16), dataspace);

    // Prepare data
    char antArray[solAntennas.size()][16];
    for (uint i=0; i<solAntennas.size(); ++i) {
      std::strncpy(antArray[i], solAntennas[i].c_str(), 16);
    }

    dataset.write(antArray, H5::StrType(H5::PredType::C_S1, 16));
  }

  void H5Parm::setSolSources(const std::string& solName,
                             const std::vector<std::string>& solSources) {
    H5::Group solTab(_solSet.openGroup(solName));
    // TODO: assert that source is present in source table of solution set
    hsize_t dims[1];
    dims[0]=solSources.size();

    // Create dataset
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = solTab.createDataSet("dir", H5::StrType(H5::PredType::C_S1, 128), dataspace);

    // Prepare data
    char srcArray[solSources.size()][128];
    for (uint i=0; i<solSources.size(); ++i) {
      std::strncpy(srcArray[i], solSources[i].c_str(), 128);
    }

    dataset.write(srcArray, H5::StrType(H5::PredType::C_S1, 128));
  }

  void H5Parm::setFreqs(const std::string& solName,
                        const std::vector<double>& freqs) {
    H5::Group solTab(_solSet.openGroup(solName));
    hsize_t dims[1];
    dims[0]=freqs.size();

    // Create dataset
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = solTab.createDataSet("freq", H5::PredType::IEEE_F64LE, dataspace);

    dataset.write(&(freqs[0]), H5::PredType::IEEE_F64LE);
  }

  void H5Parm::setTimes(const std::string& solName,
                        const std::vector<double>& times) {
    H5::Group solTab(_solSet.openGroup(solName));
    hsize_t dims[1];
    dims[0]=times.size();

    // Create dataset
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = solTab.createDataSet("time", H5::PredType::IEEE_F64LE, dataspace);

    dataset.write(&(times[0]), H5::PredType::IEEE_F64LE);
  }
}
