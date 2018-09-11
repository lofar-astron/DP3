//# PythonStep.h: A DPStep executed in some python module
//# Copyright (C) 2015
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
//# $Id: DPStep.h 30718 2015-01-19 15:31:51Z diepen $
//#
//# @author Ger van Diepen

#ifndef DPPP_PYTHONSTEP_H
#define DPPP_PYTHONSTEP_H

// @file
// @brief Class to execute a DPStep in some Python module

#include "../DPPP/DPInput.h"
#include "../DPPP/DPBuffer.h"
#include "../DPPP/DPInfo.h"

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"

#include <boost/python.hpp>

#include <casacore/casa/Containers/ValueHolder.h>
#include <casacore/casa/Containers/Record.h>

namespace DP3 {
  namespace DPPP {

    // @ingroup NDPPP

    // This class defines a step in the DPPP pipeline to be executed in
    // Python.
    // the PythonStep functionality is not directly part of the DPPP program,
    // but is automatically loaded on demand from a shared library.

    // The control flow of a python step is as follows:
    // <ol>
    //  <li> PythonStep derives from DPStep as all other DPPP steps do.
    //   Its constructor initializes Python, opens the module given in
    //   the parset key 'python.module' and creates an instance of the
    //   class given in 'python.class'.
    //   That class must derive from the Python class DPStepBase defined
    //   in __init.py__. The C++ class DPStepBase is the interface for
    //   callbacks from Python to C++.
    //  <li> The updateInfo function calls its Python counterpart. The
    //   DPInfo object is encoded in a dict. Currently, Measure objects
    //   in DPInfo are not handled yet.
    //  <li> The process function calls its Python counterpart with the
    //   time and exposure in the DPBuffer object. The Python class can
    //   obtain the other dats (visibilities, flags, weights, UVW, modeldata)
    //   by means of explicit callbacks. This was done for 2 reasons:
    //   <br>- it makes it possible to directly fill a numpy array.
    //   <br>- only data really needed is sent.
    //  <li> When the Python step has output ready, it needs to call the
    //   processNext function with the data that has changed. Note that
    //   (as in e.g. class Averager) it is possible that only every N
    //   input buffers result in an output buffer.
    //  <li> The finish, show, showCounts, showTimings, and addToMS
    //   functions call their Python counterparts. The Python base class
    //   offers default implementations, so they do not need to be
    //   implemented in a derived Python step class.
    // </ol>
    // The class also contains several functions that are basically the
    // callback function for Python.

    class PythonStep: public DPStep
    {
    public:
      PythonStep (DPInput* input,
                  const ParameterSet& parset,
                  const string& prefix);

      virtual ~PythonStep();

      // The 'constructor' for dynamically loaded steps.
      // It is registered when the shared library is loaded.
      static DPStep::ShPtr makeStep (DPInput*, const ParameterSet&,
                                     const std::string&);

      // Process the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Add some data to the MeasurementSet written/updated.
      virtual void addToMS (const string& msName);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the flag counts if needed.
      virtual void showCounts (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

      // Callback functions for Python.
      // <group>
      // Tell that the Python step needs the visibility data.
      void setNeedVisData();
      // Tell that the Python step needs data to be written.
      void setNeedWrite();
      // Get the data into the given Complex array.
      void getData (const casacore::ValueHolder&);
      // Get the flags into the given bool array.
      void getFlags (const casacore::ValueHolder&);
      // Get the weights into the given float array.
      void getWeights (const casacore::ValueHolder&);
      // Get the UVWs into the given double array.
      void getUVW (const casacore::ValueHolder&);
      // Get the model data into the given Complex array.
      void getModelData (const casacore::ValueHolder&);
      // Execute the process function of the next step.
      // The record should contain the changed buffer fields.
      bool processNext (const casacore::Record&);
      // </group>

    private:
      //# Data members.
      DPInput*     itsInput;
      std::string  itsName;
      ParameterSet itsParset;
      DPBuffer     itsBufIn;
      DPBuffer     itsBufTmp;
      DPBuffer     itsBufOut;
      bool         itsNChanChg;
      bool         itsNBlChg;
      std::string  itsPythonClass;
      std::string  itsPythonModule;
      NSTimer      itsTimer;
      boost::python::object itsPyObject;
    };

  } //# end namespace
}

// Define the function (without name mangling) to register the 'constructor'.
extern "C"
{
  void register_pythondppp();
}

#endif
