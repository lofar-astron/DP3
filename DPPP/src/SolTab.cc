#include <lofar_config.h>
#include <DPPP/H5Parm.h>
#include <Common/Exception.h>
#include <Common/StringUtil.h>
#include <Common/LofarLogger.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <iomanip>

#include <hdf5.h>

using namespace std;

namespace LOFAR {
  H5Parm::SolTab::SolTab(H5::Group group,
                         const std::string& type,
                         const std::vector<AxisInfo> axes):
                           H5::Group(group),
                           _type(type),
                           _axes(axes)
  {
    H5::Attribute attr = createAttribute("TITLE",
                             H5::StrType(H5::PredType::C_S1, _type.size()),
                             H5::DataSpace());
    attr.write(H5::StrType(H5::PredType::C_S1, _type.size()), _type);
    addVersionStamp(*this);
  }

  H5Parm::SolTab::SolTab(H5::Group& group):
    H5::Group(group)
  {
    // Read the type from the "TITLE" attribute
    H5::Attribute typeattr = openAttribute("TITLE");
    hsize_t typenamelen = typeattr.getDataType().getSize();
    char typecstr[typenamelen+1];
    typecstr[typenamelen]='\0';
    typeattr.read(typeattr.getDataType(), &typecstr);
    _type = typecstr;

    readAxes();
  }

  H5Parm::SolTab::~SolTab() {
  }

  H5Parm::AxisInfo H5Parm::SolTab::getAxis(uint i) {
    return _axes[i];
  }

  H5Parm::AxisInfo H5Parm::SolTab::getAxis(const std::string& axisName) {
    for (uint i=0; i<_axes.size(); ++i) {
      if (_axes[i].name == axisName) {
        return _axes[i];
      }
    }
    THROW(Exception, "Axis "<<axisName<<" does not exist in "<<getName());
  }

  bool H5Parm::SolTab::hasAxis(const string& axisName) {
    for (size_t i=0; i<_axes.size(); ++i) {
      if (_axes[i].name==axisName)
        return true;
    }
    return false;
  }

  size_t H5Parm::SolTab::getAxisIndex(const std::string& axisName) {
    for (uint i=0; i<_axes.size(); ++i) {
      if (_axes[i].name == axisName) {
        return i;
      }
    }
    THROW(Exception, "Axis "<<axisName<<" does not exist in "<<getName());
  }

  void H5Parm::SolTab::setValues(const std::vector<double>& vals,
                                 const std::vector<double>& weights) {
    // Convert axes to comma separated string, fill dims
    string axesstr = _axes[0].name;
    vector<hsize_t> dims(_axes.size());
    for (uint i=0; i<_axes.size(); ++i) {
      dims[i] = _axes[i].size;
      if (i>0) {
        axesstr += ","+_axes[i].name;
      }
    }

    H5::DataSpace dataspace(dims.size(), &(dims[0]), NULL);
    H5::DataSet dataset = createDataSet("val", H5::PredType::IEEE_F64LE,
                                        dataspace);

    dataset.write(&(vals[0]), H5::PredType::IEEE_F64LE);

    H5::Attribute attr = dataset.createAttribute("AXES",
                             H5::StrType(H5::PredType::C_S1, axesstr.size()),
                             H5::DataSpace());
    attr.write(H5::StrType(H5::PredType::C_S1, axesstr.size()), axesstr);

    // Add weights
    H5::DataSet weightset = createDataSet("weight", H5::PredType::IEEE_F64LE,
                                          dataspace);

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

  void H5Parm::SolTab::setComplexValues(const std::vector<std::complex<double> >& vals,
                                       const std::vector<double>& weights,
                                       bool toAmplitudes) {
    // Convert values to real numbers by taking amplitude or argument
    vector<double> realvals(vals.size());

    if (toAmplitudes) {
      transform(vals.begin(), vals.end(), realvals.begin(), takeAbs);
    } else { // Phase only
      transform(vals.begin(), vals.end(), realvals.begin(), takeArg);
    }

    setValues(realvals, weights);
  }

  void H5Parm::SolTab::readAxes() {
    H5::DataSet val;
    try {
      val = openDataSet("val");
    } catch (H5::GroupIException& e) {
      THROW(Exception, "SolTab "<<getName()<<" has no values");
    }

    H5::Attribute axesattr;
    try {
      axesattr = val.openAttribute("AXES");
    } catch (H5::AttributeIException& e) {
      THROW(Exception, "Values of SolTab "<<getName()<<" has no AXIS attribute");
    }

    hsize_t axesstrlen = axesattr.getDataType().getSize();
    char axescstr[axesstrlen+1];
    axescstr[axesstrlen]='\0';
    axesattr.read(axesattr.getDataType(), &axescstr);
    vector<string> axesnames = StringUtil::tokenize(axescstr,",");

    uint ndims = axesnames.size();

    // Get number of dimensions and size of all dimensions
    H5::DataSpace ds = val.getSpace();
    ASSERT (ds.getSimpleExtentNdims() == int(ndims));
    hsize_t dims_out[ndims];
    ds.getSimpleExtentDims(dims_out);

    for (uint i=0; i<axesnames.size(); ++i) {
      AxisInfo a(axesnames[i], dims_out[i]);
      _axes.push_back(a);
    }
  }

  string H5Parm::SolTab::getName() const {
    size_t len = H5Iget_name(getId(),NULL,0);
    char buffer[len];
    H5Iget_name(getId(),buffer,len+1);
    // Strip leading /
    return buffer+1;
  }

  vector<double> H5Parm::SolTab::getValues(
              const std::string& antName,
              uint starttimeslot, uint ntime, uint timestep,
              uint startfreq, uint nfreq, uint freqstep,
              uint pol, uint dir) {
    vector<double> res(ntime*nfreq);
    H5::DataSet val = openDataSet("val");

    // Set offsets and strides
    hsize_t memdims[_axes.size()];
    hsize_t offset[_axes.size()];
    hsize_t count[_axes.size()];
    hsize_t stride[_axes.size()];

    bool timeAxisFound = false;
    for (uint i=0; i<_axes.size(); ++i) {
      stride[i] = 1;
      count[i] = 1;
      memdims[i] = 1;
      if (_axes[i].name=="time") {
        offset[i] = starttimeslot;
        stride[i] = timestep;
        count[i] = ntime;
        memdims[i] = ntime;
        timeAxisFound = true;
      } else if (_axes[i].name=="freq") {
        offset[i] = startfreq;
        stride[i] = freqstep;
        count[i] = nfreq;
        memdims[i] = nfreq;
      } else if (_axes[i].name=="ant") {
        offset[i] = getAntIndex(antName);
      } else if (_axes[i].name=="dir") {
        offset[i] = dir;
      } else if (_axes[i].name=="pol") {
        offset[i] = pol;
      } else {
        ASSERT(_axes[i].size == 1);
        offset[i] = 0;
      }
    }

    ASSERT(timeAxisFound);

    H5::DataSpace dataspace = val.getSpace();

    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride);

    // Setup memory dataspace
    H5::DataSpace memspace(_axes.size(), memdims);
    try {
      val.read(&(res[0]), H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
    } catch (H5::DataSetIException& e) {
      e.printError();
      THROW(Exception, "Could not read data");
    }
    return res;
  }

  void H5Parm::SolTab::setAntennas(const std::vector<std::string>& solAntennas) {
    // TODO: assert that antenna is present in antenna table in solset
    hsize_t dims[1];
    dims[0]=solAntennas.size();

    // Create dataset
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = createDataSet("ant", H5::StrType(H5::PredType::C_S1, 16), dataspace);

    // Prepare data
    char antArray[solAntennas.size()][16];
    for (uint i=0; i<solAntennas.size(); ++i) {
      std::strncpy(antArray[i], solAntennas[i].c_str(), 16);
    }

    dataset.write(antArray, H5::StrType(H5::PredType::C_S1, 16));
  }

  void H5Parm::SolTab::setSources(const std::vector<std::string>& solSources) {
    // TODO: assert that source is present in source table of solset
    hsize_t dims[1];
    dims[0]=solSources.size();

    // Create dataset
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = createDataSet("dir",
                                        H5::StrType(H5::PredType::C_S1, 128),
                                        dataspace);

    // Prepare data
    char srcArray[solSources.size()][128];
    for (uint i=0; i<solSources.size(); ++i) {
      std::strncpy(srcArray[i], solSources[i].c_str(), 128);
    }

    dataset.write(srcArray, H5::StrType(H5::PredType::C_S1, 128));
  }

  void H5Parm::SolTab::setFreqs(const std::vector<double>& freqs) {
    setAxisMeta("freq", freqs);
  }

  void H5Parm::SolTab::setTimes(const std::vector<double>& times) {
    setAxisMeta("time", times);
  }

  void H5Parm::SolTab::setAxisMeta(const std::string& metaName,
                                   const std::vector<double>& metaVals) {
    hsize_t dims[1];
    dims[0]=metaVals.size();

    // Create dataset
    H5::DataSpace dataspace(1, dims, NULL);
    H5::DataSet dataset = createDataSet(metaName,
                                        H5::PredType::IEEE_F64LE, dataspace);

    dataset.write(&(metaVals[0]), H5::PredType::IEEE_F64LE);
  }

  hsize_t H5Parm::SolTab::getAntIndex(const string& antName) {
    return getNamedIndex(_antMap, "ant", antName);
  }

  hsize_t H5Parm::SolTab::getNamedIndex(map<string, hsize_t>& cache,
                                        const string& tableName,
                                        const string& elementName) {
    // Initialize _antMap on first use
    if (cache.empty()) {
      H5::DataSet dataset;
      H5::DataSpace dataspace;
      try {
        dataset = openDataSet(tableName);
        dataspace = dataset.getSpace();
      } catch (H5::GroupIException& e) {
        THROW(Exception, "SolTab has no table "<<tableName);
      }
      ASSERT(dataspace.getSimpleExtentNdims()==1);
      hsize_t dims[1];
      dataspace.getSimpleExtentDims(dims);

      // TODO: assert that DataType is String
      hsize_t strLen = dataset.getDataType().getSize();

      char elNames[strLen*dims[0]];
      dataset.read(elNames, H5::StrType(H5::PredType::C_S1, strLen));

      for (hsize_t elNum=0; elNum<dims[0];++elNum) {
        char elNamecstr[strLen+1];
        std::copy(elNames+elNum*strLen,elNames+(elNum+1)*strLen, elNamecstr);
        elNamecstr[strLen]='\0';
        cache[elNamecstr]=elNum;
      }
    }
    map<string, hsize_t>::iterator it = cache.find(elementName);
    if (it == cache.end()) {
      THROW(Exception, "SolTab has no element "<<elementName <<
                       " in "<<tableName);
    }
    return cache.find(elementName)->second;
  }

  hsize_t H5Parm::SolTab::getFreqIndex(double freq) const {
    H5::DataSet freqset;
    H5::DataSpace dataspace;

    try {
      freqset = openDataSet("freq");
      dataspace = freqset.getSpace();
    } catch (H5::GroupIException& e) {
      THROW(Exception, "SolTab "<<getName()<<" has no freq table");
    }
    ASSERT(dataspace.getSimpleExtentNdims()==1);

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);

    if (dims[0]==1) {
      return 0;
    }
    vector<double> freqs(dims[0]);
    freqset.read(&(freqs[0]), H5::PredType::NATIVE_DOUBLE);

    for (size_t i = 0; i<freqs.size()-1; ++i) {
      if (freq>freqs[i] && freq<freqs[i+1]) {
        // Nearest neighbor: i or i+1
        if (freq-freqs[i] < freqs[i+1]-freq) {
          return i;
        }
        else {
          return i+1;
        }
      }
    }
    THROW(Exception,"Frequency "<<fixed<<freq<<" not found in "<<getName());
    return 0;
  }

  hsize_t H5Parm::SolTab::getTimeIndex(double time) const {
    H5::DataSet timeset;
    H5::DataSpace dataspace;
    try {
      timeset = openDataSet("time");
      dataspace = timeset.getSpace();
    } catch (H5::GroupIException& e) {
      THROW(Exception, "SolTab "<<getName()<<" has no time table");
    }
    ASSERT(dataspace.getSimpleExtentNdims()==1);

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);

    vector<double> times(dims[0]);
    timeset.read(&(times[0]), H5::PredType::NATIVE_DOUBLE);

    double timeInterval = getTimeInterval();

    for (size_t i = 0; i<times.size(); ++i) {
      if (abs(times[i]-time)<timeInterval*0.5)
        return i;
    }
    THROW(Exception,"Time "<<fixed<<time<<" not found in "<<getName());
    return 0;
  }

  hsize_t H5Parm::SolTab::getDirIndex(const string& dirName) {
    return getNamedIndex(_dirMap, "dir", dirName);
  }

  double H5Parm::SolTab::getInterval(const std::string& axisName) const {
    H5::DataSet dataset;
    H5::DataSpace dataspace;
    try {
      dataset = openDataSet(axisName);
      dataspace = dataset.getSpace();
    } catch (H5::GroupIException& e) {
      THROW(Exception, "SolTab "<<getName()<<" has no axis table for "<<axisName);
    }
    ASSERT(dataspace.getSimpleExtentNdims()==1);

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);
    ASSERTSTR(dims[0]>1, "For reading the interval, more than one value is required.");

    hsize_t count[1], offset[1];
    count[0]=2; offset[0]=0;
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    // Get only two values
    vector<double> values(2);
    dataset.read(&(values[0]), H5::PredType::NATIVE_DOUBLE, dataspace, dataspace);

    return values[1]-values[0];
  }
}
