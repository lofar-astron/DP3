//# PythonStep.cc: A DPStep executed in some python module
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
//# $Id: ApplyCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <PythonDPPP/PythonStep.h>
#include <PythonDPPP/DPStepBase.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/DPRun.h>
#include <Common/ParameterSet.h>

#if defined(casacore)
#include <python/Converters/PycExcp.h>
#include <python/Converters/PycBasicData.h>
#include <python/Converters/PycValueHolder.h>
#include <python/Converters/PycRecord.h>
#include <python/Converters/PycArray.h>
#define pyrap python
#else
#include <pyrap/Converters/PycExcp.h>
#include <pyrap/Converters/PycBasicData.h>
#include <pyrap/Converters/PycValueHolder.h>
#include <pyrap/Converters/PycRecord.h>
#include <pyrap/Converters/PycArray.h>
#endif

#include <casa/OS/Path.h>
#include <unistd.h>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    PythonStep::PythonStep (DPInput* input,
                            const ParameterSet& parset,
                            const string& prefix)
      : itsInput        (input),
        itsName         (prefix),
        itsParset       (parset.makeSubset (prefix)),
        itsNChanChg     (false),
        itsNBlChg       (false),
        itsPythonClass  (itsParset.getString ("python.class")),
        itsPythonModule (itsParset.getString ("python.module", itsPythonClass))
    {
      // Initialize Python interpreter.
      // Note: a second call is a no-op.
      Py_Initialize();
      // Insert the current working directory into the python path.
      // (from http://stackoverflow.com/questions/9285384/
      //   how-does-import-work-with-boost-python-from-inside-python-files)
      string workingDir = Path(".").absoluteName();
      char path[] = "path";    // needed to avoid warning if "path" used below
      PyObject* sysPath = PySys_GetObject(path);
      PyList_Insert (sysPath, 0, PyString_FromString(workingDir.c_str()));
      // Register converters for casa types from/to python types
      casa::pyrap::register_convert_excp();
      casa::pyrap::register_convert_basicdata();
      casa::pyrap::register_convert_casa_valueholder();
      casa::pyrap::register_convert_casa_record();
      try {
        // First import main
        boost::python::object mainModule = boost::python::import
          ("__main__");
        // Import the given module
        boost::python::object dpppModule = boost::python::import
          (itsPythonModule.c_str());
        // Get the python class object from the imported module
        boost::python::object dpppAttr = dpppModule.attr
          (itsPythonClass.c_str());

        // Convert the ParameterSet to a Record (using its string values).
        Record rec;
        for (ParameterSet::const_iterator iter=itsParset.begin();
             iter!=itsParset.end(); ++iter) {
          rec.define (iter->first, iter->second.get());
        }
        // Create an instance of the python class passing the record.
        itsPyObject = dpppAttr(rec);
        // Set the pointer to this object in the DPStepBase object.
        DPStepBase::theirPtr->setStep (this);
      } catch (boost::python::error_already_set const &) {
        // handle the exception in some way
        PyErr_Print();
        throw;
      }
    }

    PythonStep::~PythonStep()
    {
    }

    DPStep::ShPtr PythonStep::makeStep (DPInput* input,
                                        const ParameterSet& pset,
                                        const std::string& prefix)
    {
      return DPStep::ShPtr(new PythonStep(input, pset, prefix));
    }

    void PythonStep::updateInfo (const DPInfo& infoIn)
    {
      try {
        boost::python::object result =
          itsPyObject.attr("_updateInfo")(infoIn.toRecord());
        Record rec = boost::python::extract<Record>(result);
        info() = infoIn;
        // Merge possible result back in DPInfo object.
        info().fromRecord (rec);
        // Check if nr of channels has changed.
        // Note: currently a change in antennae/baselines is not supported
        // because the antenna positions are not passed yet.
        if (infoIn.nchan() != info().nchan()) {
          itsNChanChg = true;
        }
        ASSERT (infoIn.getAnt1().size() == info().getAnt1().size()  &&
                infoIn.getAnt2().size() == info().getAnt2().size()  &&
                allEQ(infoIn.getAnt1(), info().getAnt1())  &&
                allEQ(infoIn.getAnt2(), info().getAnt2()));
        boost::python::object res1 = itsPyObject.attr("needVisData")();
        if (boost::python::extract<bool>(res1)) {
          info().setNeedVisData();
        }
        boost::python::object res2 = itsPyObject.attr("needWrite")();
        if (boost::python::extract<bool>(res2)) {
          info().setWriteData();
        }
      } catch (boost::python::error_already_set const &) {
        // handle the exception in some way
        PyErr_Print();
        throw;
      }
    }

    bool PythonStep::process (const DPBuffer& buf)
    {
      try {
        itsTimer.start();
        itsBufIn.referenceFilled (buf);
        boost::python::object result =
          itsPyObject.attr("process")(itsBufIn.getTime(),
                                      itsBufIn.getExposure());
        bool res = boost::python::extract<bool> (result);
        itsTimer.stop();
        return res;
      } catch (boost::python::error_already_set const &) {
        // handle the exception in some way
        PyErr_Print();
        throw;
      }
    }

    void PythonStep::finish()
    {
      try {
        itsPyObject.attr("finish")();
        getNextStep()->finish();
      } catch (boost::python::error_already_set const &) {
        // handle the exception in some way
        PyErr_Print();
        throw;
      }
    }

    void PythonStep::addToMS (const string& msName)
    {
      try {
        itsPyObject.attr("addToMS")(msName);
      } catch (boost::python::error_already_set const &) {
        // handle the exception in some way
        PyErr_Print();
        throw;
      }
    }

    void PythonStep::show (std::ostream& os) const
    {
      try {
        boost::python::object result = itsPyObject.attr("show")();
        string str = boost::python::extract<string>(result);
        os << "PythonStep " << itsName << " class=" << itsPythonClass << endl;
        if (! str.empty()) {
          os << str;
        }
      } catch (boost::python::error_already_set const &) {
        // handle the exception in some way
        PyErr_Print();
        throw;
      }
    }

    void PythonStep::showCounts (std::ostream& os) const
    {
      try {
        boost::python::object result = itsPyObject.attr("showCounts")();
        string str = boost::python::extract<string>(result);
        if (! str.empty()) {
          os << str;
        }
      } catch (boost::python::error_already_set const &) {
        // handle the exception in some way
        PyErr_Print();
        throw;
      }
    }

    void PythonStep::showTimings (std::ostream& os, double duration) const
    {
      try {
        os << "  ";
        FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
        os << " PythonStep " << itsName << " class=" << itsPythonClass << endl;
        boost::python::object result =
          itsPyObject.attr("showTimings")(itsTimer.getElapsed());
        string str = boost::python::extract<string>(result);
        if (! str.empty()) {
          os << str;
        }
      } catch (boost::python::error_already_set const &) {
        // handle the exception in some way
        PyErr_Print();
        throw;
      }
    }
      
    // Implement all functions to communicate with python.
    void PythonStep::setNeedVisData()
    {
      info().setNeedVisData();
    }
    void PythonStep::setNeedWrite()
    {
      info().setWriteData();
    }

    void PythonStep::getData (const ValueHolder& vh)
    {
      ASSERT (vh.dataType() == TpArrayComplex);
      Array<Complex> arr(vh.asArrayComplex());
      arr = itsBufIn.getData();   // operator= checks if shape matches
    }

    void PythonStep::getFlags (const ValueHolder& vh)
    {
      ASSERT (vh.dataType() == TpArrayBool);
      Array<Bool> arr(vh.asArrayBool());
      arr = itsBufIn.getFlags();
    }

    void PythonStep::getWeights (const ValueHolder& vh)
    {
      ASSERT (vh.dataType() == TpArrayFloat);
      Array<Float> arr(vh.asArrayFloat());
      arr = itsInput->fetchWeights (itsBufIn, itsBufTmp, itsTimer);
    }

    void PythonStep::getUVW (const ValueHolder& vh)
    {
      ASSERT (vh.dataType() == TpArrayDouble);
      Array<Double> arr(vh.asArrayDouble());
      arr = itsInput->fetchUVW (itsBufIn, itsBufTmp, itsTimer);
    }

    void PythonStep::getModelData (const ValueHolder& vh)
    {
      ASSERT (vh.dataType() == TpArrayComplex);
      Cube<Complex> arr(vh.asArrayComplex());
      itsInput->getModelData (RefRows(itsBufIn.getRowNrs()), arr);
    }

    bool PythonStep::processNext (const Record& rec)
    {
      itsTimer.stop();
      uint nproc = 0;
      uint narr  = 0;
      if (rec.isDefined("TIME")) {
        itsBufOut.setTime (rec.asDouble("TIME"));
        nproc++;
      } else {
        itsBufOut.setTime (itsBufIn.getTime());
      }
      if (rec.isDefined("EXPOSURE")) {
        itsBufOut.setExposure (rec.asDouble("EXPOSURE"));
        nproc++;
      } else {
        itsBufOut.setExposure (itsBufIn.getExposure());
      }
      itsBufOut.setRowNrs (itsBufIn.getRowNrs());
      if (rec.isDefined("DATA")) {
        itsBufOut.getData().assign (rec.toArrayComplex("DATA"));
        narr++;
      } else if (! itsNChanChg  &&  ! itsNBlChg) {
        itsBufOut.getData().assign (itsBufIn.getData());
      }
      if (rec.isDefined("FLAGS")) {
        itsBufOut.getFlags().assign (rec.toArrayBool("FLAGS"));
        narr++;
      } else if (! itsNChanChg  &&  ! itsNBlChg) {
        itsBufOut.getFlags().assign (itsBufIn.getFlags());
      }
      if (rec.isDefined("WEIGHTS")) {
        itsBufOut.getWeights().assign (rec.toArrayFloat("WEIGHTS"));
        narr++;
      } else if (! itsNChanChg  &&  ! itsNBlChg) {
        if (! itsBufIn.getWeights().empty()) {
          itsBufOut.getWeights().assign (itsBufIn.getWeights());
        }
      }
      if (rec.isDefined("UVW")) {
        itsBufOut.getUVW().assign (rec.toArrayDouble("UVW"));
        narr++;
      } else if (! itsNChanChg  &&  ! itsNBlChg) {
        if (! itsBufIn.getUVW().empty()) {
          itsBufOut.getUVW().assign (itsBufIn.getUVW());
        }
      }
      nproc += narr;
      if (nproc != rec.nfields()) {
        THROW (Exception,
               "Record/dict given to processNext() contains unknown fields");
      }
      if ((itsNChanChg  ||  itsNBlChg)  &&  narr != 4) {
        THROW (Exception,
               "Record/dict given to processNext() must contain DATA, FLAGS, "
               "WEIGHTS, and UVW if the nr of channels or baselines changes");
      }
      bool res = getNextStep()->process (itsBufOut);
      itsTimer.start();
      return res;
    }

  } //# end namespace
}

// Define the function to make the PythonStep 'constructor' known.
void register_pythondppp()
{
  LOFAR::DPPP::DPRun::registerStepCtor ("pythondppp",
                                        LOFAR::DPPP::PythonStep::makeStep);
}
