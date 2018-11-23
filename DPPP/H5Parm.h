#ifndef DPPP_H5PARM_H
#define DPPP_H5PARM_H

#include <string>
#include <vector>
#include <complex>
#include <map>

#include <H5Cpp.h>

#include <utility>

namespace DP3 {
  class H5Parm : private H5::H5File
  {
    typedef struct antenna_t {
      char name[16];
      float position[3];
    } antenna_t;

    typedef struct source_t {
      char name[128];
      float dir[2];
    } source_t;
    
    typedef struct polarization_t {
      char name[2];
    } polarization_t;

    public:
    // A name and the length of an exis, e.g. ('freq', 800) for 800 frequencies
    struct AxisInfo {
      public: AxisInfo(const std::string& name, uint size) :
        name(name), size(size) {};

        std::string name;
        uint size;
      };

      // SolTab is a solution table as defined in the H5Parm standard. It
      // contains one solution, e.g. all TEC values, with different axes
      // for that solution (e.g. time, freq, pol).
      class SolTab : private H5::Group {
        public:
          // Default constructor
          SolTab() {};

          // Create a new soltab, add it to its parent
          SolTab(H5::Group group,
                 const std::string& type,
                 const std::vector<AxisInfo> axes // Axes, fastest varying last
                 );

          // Create a soltab from a H5::Group (for reading existing files)
          SolTab(H5::Group& group);

          // The destructor could check for valid subtables
          virtual ~SolTab();

          std::vector<AxisInfo>& getAxes() {return _axes;}

          AxisInfo getAxis(uint i) const;

          // Get an axis, throw an exception if it does not exist
          AxisInfo getAxis(const std::string& axisName) const;

          size_t nAxes() { return _axes.size(); }

          bool hasAxis(const std::string& axisName);

          // Get the index of an axis
          size_t getAxisIndex(const std::string& axisname);

          void setAntennas(const std::vector<std::string>& solAntennas);

          void setSources(const std::vector<std::string>& solSources);

          void setPolarizations(const std::vector<std::string>& polarizations);

          void setFreqs(const std::vector<double>& freqs);

          // Get the values of a real-valued axis (e.g. "time" or "freq")
          std::vector<double> getRealAxis(const std::string& axisname) const;

          // Get the values of a string-valued axis (e.g. "dir" or "pol")
          std::vector<std::string> getStringAxis(const std::string& axisname);

          // Get the index of freq, using nearest neighbor
          // This assumes that the frequencies are in increasing order.
          hsize_t getFreqIndex(double freq) const;

          // Get the index of a time. Matches with 0.5*timeInterval
          // This assumes that all times are regularly spaced
          hsize_t getTimeIndex(double time) const;

          hsize_t getDirIndex(const std::string& dirName);

          // Gets the interval (in s.) between a time slot (default first) and
          // the next. Throws error if there is only one time slot.
          double getTimeInterval(size_t start=0) const {
            return getInterval("time", start);
          }

          // Gets the interval (in s.) between a channel (default first) and
          // the next. Throws error if there is only one frequency.
          double getFreqInterval(size_t start=0) const {
            return getInterval("freq", start);
          }

          void setTimes(const std::vector<double>& times);

          // Set metadata about an axis (like freq or time))
          void setAxisMeta(const std::string& metaName,
                           const std::vector<double>& metaVals);

          // Set metadata about an axis (like polarization, direction)
          void setAxisMeta(const std::string& metaName,
                           size_t strLen,
                           const std::vector<std::string>& metaVals);

          // Adds a real solution.
          // If weights are emtpy, write ones everywhere
          void setValues(const std::vector<double>& vals,
                         const std::vector<double>& weights,
                         const std::string& history="");

          // Add a complex solution, taking either amplitude or phase
          void setComplexValues(const std::vector<std::complex<double> >& vals,
                                const std::vector<double>& weights,
                                bool toAmplitudes, const std::string& history="");



          // Get the name of this SolTab
          std::string getName() const;

          std::string getType() const {return _type;}

          // Get the values of this SolTab for a given antenna.
          std::vector<double> getValues(
                                        const std::string& antName,
                                        uint starttimeslot, uint ntime, uint timestep,
                                        uint startfreq, uint nfreq, uint freqstep,
                                        uint pol, uint dir) {
            return getValuesOrWeights("val", antName,
                                      starttimeslot, ntime, timestep,
                                      startfreq, nfreq, freqstep,
                                      pol, dir);
          }

          // Get the weights of this SolTab for a given antenna.
          std::vector<double> getWeights(
                                        const std::string& antName,
                                        uint starttimeslot, uint ntime, uint timestep,
                                        uint startfreq, uint nfreq, uint freqstep,
                                        uint pol, uint dir) {
            return getValuesOrWeights("weight", antName,
                                      starttimeslot, ntime, timestep,
                                      startfreq, nfreq, freqstep,
                                      pol, dir);
          }

          std::vector<double> getValuesOrWeights(
                                        const std::string& valOrWeight,
                                        const std::string& antName,
                                        const std::vector<double>& times,
                                        const std::vector<double>& freqs,
                                        uint pol, uint dir, bool nearest);
        private:
          // Get the values or weights of this SolTab for a given antenna.
          std::vector<double> getValuesOrWeights(
                                        const std::string& valOrWeight,
                                        const std::string& antName,
                                        uint starttimeslot, uint ntime, uint timestep,
                                        uint startfreq, uint nfreq, uint freqstep,
                                        uint pol, uint dir);

          void readAxes();

          void fillCache(std::map<std::string, hsize_t>& cache,
                         const std::string& tableName);

          // Get the interval of the axis axisName
          double getInterval(const std::string& axisName, size_t start=0) const;
          hsize_t getAntIndex(const std::string& antName);
          hsize_t getNamedIndex(std::map<std::string, hsize_t>& cache,
                                const std::string& tableName,
                                const std::string& elementName);

          std::string _type;
          std::vector<AxisInfo> _axes;
          std::map<std::string, hsize_t> _antMap;
          std::map<std::string, hsize_t> _dirMap;
      };

      // Open existing H5Parm or create a new one
      // Default name is given by solSetName, if that does not exist continue
      // searching for sol000, sol001, etc.
      // Only one solset of an H5Parm can be opened at once; this object only
      // provides info about one SolSet (even though the file can contain more).
      H5Parm(const std::string& filename, bool forceNew=false, 
             bool forceNewSolSet=false, const std::string& solSetName="");


      H5Parm();

      virtual ~H5Parm();

      // Add metadata (directions on the sky in J2000) about named sources
      void addSources (const std::vector<std::string>& names,
                       const std::vector<std::pair<double, double> >& dirs);

      // Add metadata (positions on earth in ITRF) about antennas
      void addAntennas(const std::vector<std::string>& names,
                       const std::vector<std::vector<double> >& positions);
      
      // Add metadata about polarizations
      void addPolarizations(const std::vector<std::string>& polarizations);

      // Add a version stamp in the attributes of the group
      static void addVersionStamp(H5::Group &node);

      // Create and return a new soltab. Type is the type as used in BBS
      SolTab& createSolTab(const std::string& name,
                           const std::string& type,
                           const std::vector<AxisInfo> axes);

      SolTab& getSolTab(const std::string& name);

      // Get the name of the one SolSet used in this H5Parm
      std::string getSolSetName() const ;

      // Get the number of SolTabs in the active solset of this h5parm
      size_t nSolTabs() { return _solTabs.size(); }

      // Is the given soltab resent in the active solset of this h5parm
      bool hasSolTab(const std::string& solTabName) const;
    private:

      static double takeAbs(std::complex<double> c) {
        return std::abs(c);
      }
      static double takeArg(std::complex<double> c) {
        return std::arg(c);
      }

      std::map<std::string, SolTab> _solTabs;
      H5::Group _solSet;
  };
}

#endif
