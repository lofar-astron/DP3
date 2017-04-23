#ifndef DPPP_H5PARM_H
#define DPPP_H5PARM_H


#include <string>

#include <H5Cpp.h>

#include <utility>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/measures/Measures/MPosition.h>

namespace LOFAR {
  class H5Parm
  {
    typedef struct antenna_t {
      char name[16];
      float position[3];
    } antenna_t;

    typedef struct source_t {
      char name[128];
      float dir[2];
    } source_t;

    public:
      H5Parm(const std::string& name);

      ~H5Parm();

      void addSources (const std::vector<std::string>& names,
                       const std::vector<std::pair<double, double> >& dirs);

      void addAntennas(const std::vector<std::string>& names,
                       const std::vector<std::vector<double> >& positions);

      // Add a solution
      // If weights are emtpy, write ones everywhere
      void addSolution(const std::string& solName,
                       const std::string& solType,
                       const std::string& axesstr, // String with axes names, fastest varying last
                       const std::vector<hsize_t>& dims,
                       const std::vector<double>& vals,
                       const std::vector<double>& weights);

      // Add a complex solution, taking either amplitude or phase
      void addSolution(const std::string& solName,
                       const std::string& solType,
                       const std::string& axesstr, // String with axes names, fastest varying last
                       const std::vector<hsize_t>& dims,
                       const std::vector<std::complex<double> >& vals,
                       const std::vector<double>& weights, 
                       bool toAmplitudes);

      void setSolAntennas(const std::string& solName,
                          const std::vector<std::string>& solAntennas);

      void setSolSources(const std::string& solName,
                         const std::vector<std::string>& solSources);

      void setFreqs(const std::string& solName,
                    const std::vector<double>& freqs);

      void setTimes(const std::string& solName,
                    const std::vector<double>& times);

    private:
      void addVersionStamp(H5::Group &node);
      static double takeAbs(std::complex<double> c) {
        return std::abs(c);
      }
      static double takeArg(std::complex<double> c) {
        return std::arg(c);
      }

      H5::H5File _hdf5file;
      H5::Group _solSet;
  };
}

#endif
