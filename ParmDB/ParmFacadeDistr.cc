//# ParmFacadeDistr.cc: Object access the parameter database
//#
//# Copyright (C) 2006
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
//# $Id: ParmFacadeDistr.cc 27747 2013-12-12 11:09:58Z diepen $

#include "ParmFacadeDistr.h"

#include "VdsDesc.h"

#include "../Blob/BlobArray.h"
#include "../Blob/BlobAipsIO.h"

#include <casacore/casa/IO/AipsIO.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/Utilities/GenSort.h>

using namespace DP3::CEP;
using namespace casacore;


namespace DP3 {
  namespace BBS {

    // Initialize statics.
    int ParmFacadeDistr::theirNextPort = 4195;
    vector<string> ParmFacadeDistr::theirFreePorts;


    ParmFacadeDistr::ParmFacadeDistr (const string& tableName)
      : itsPort (getPort()),
        itsConn (itsPort)
    {
      // Get info from VDS. It is automatically closed thereafter.
      int nparts;
      string cdescName;
      VdsDesc vds(tableName);
      nparts    = vds.getParts().size();
      cdescName = vds.getDesc().getClusterDescName();
      // Start all remote processes.
      string command("startparmdbdistr " + itsPort + ' ' +
                     cdescName + ' ' + tableName);
      ASSERT (system(command.c_str()) == 0);
      // Accept a connection from the clients and check if they are
      // initialized correctly.
      // Note that the order in which connections are accepted does not
      // need to be the same as the order of the parts in the VDS.
      // Chances determine the connection order.
      itsConn.addConnections (nparts);
      itsPartNames.reserve (nparts);
      BlobString buf;
      string fname;
      for (int i=0; i<itsConn.size(); ++i) {
        itsConn.read (i, buf);
        MWBlobIn bbi(buf);
        ASSERT (bbi.getOperation() == 1);    // ensure successful init
        bbi.blobStream() >> fname;           // get part name
        itsPartNames.push_back (fname);
        if (i == 0) {
          bbi.blobStream() >> itsParmNames;
          getRecord (bbi.blobStream(), itsDefValues);
        } else {
          vector<string> names;
          Record defValues;
          bbi.blobStream() >> names;
          getRecord (bbi.blobStream(), defValues);
          checkNames (itsParmNames, names, i);
          if (defValues.size() != itsDefValues.size()) {
            LOG_WARN_STR ("DEFAULTVALUES sizes of parts " << itsPartNames[0]
                          << " and " << itsPartNames[i] << " differ");
          }
        }
        bbi.finish();
      }
    }

    ParmFacadeDistr::~ParmFacadeDistr()
    {
      quit();
    }

    void ParmFacadeDistr::quit()
    {
      // Send all parts an end message.
      BlobString buf;
      MWBlobOut bbo(buf, Quit, 0);
      bbo.finish();
      itsConn.writeAll (buf);
      freePort();
    }

    string ParmFacadeDistr::getPort()
    {
      if (! theirFreePorts.empty()) {
        string port = theirFreePorts.back();
        theirFreePorts.pop_back();
        return port;
      }
      ostringstream ostr;
      ostr << theirNextPort;
      theirNextPort++;
      return ostr.str();
    }

    void ParmFacadeDistr::freePort()
    {
      theirFreePorts.push_back (itsPort);
    }

    bool ParmFacadeDistr::checkStatus (MWBlobIn& bbi, int i) const
    {
      if (bbi.getOperation() != 1) {
        string msg;
        bbi.blobStream() >> msg;
        LOG_WARN_STR ("Error '" << msg << "' in part " << itsPartNames[i]);
        return false;
      }
      return true;
    }

    void ParmFacadeDistr::checkStatusAll() const
    {
      BlobString buf;
      for (int i=0; i<itsConn.size(); ++i) {
        itsConn.read (i, buf);
        MWBlobIn bbi(buf);
        checkStatus (bbi, i);
        bbi.finish();
      }
    }

    vector<double> ParmFacadeDistr::getRange (const string& parmNamePattern) const
    {
      BlobString buf;
      MWBlobOut bbo(buf, GetRange, 0);
      bbo.blobStream() << parmNamePattern;
      bbo.finish();
      itsConn.writeAll (buf);
      vector<double> range, result;
      for (int i=0; i<itsConn.size(); ++i) {
        itsConn.read (i, buf);
        MWBlobIn bbi(buf);
        if (checkStatus (bbi, i)) {
          bbi.blobStream() >> range;
          if (result.empty()) {
            result = range;
          } else {
            if (range[0] < result[0]) result[0] = range[0];
            if (range[1] > result[1]) result[1] = range[1];
            if (range[2] < result[2]) result[2] = range[2];
            if (range[3] > result[3]) result[3] = range[3];
          }
        }
        bbi.finish();
      }
      return result;
    }

    // Get all parameter names in the table.
    vector<string> ParmFacadeDistr::getNames (const string& parmNamePattern,
                                              bool includeDefaults) const
    {
      vector<string> result;
      if (parmNamePattern.empty()  ||  parmNamePattern == "*") {
        if (!includeDefaults) {
          return itsParmNames;
        } else {
          result = itsParmNames;
        }
      } else {
        Regex regex(Regex::fromPattern(parmNamePattern));
        for (vector<string>::const_iterator iter=itsParmNames.begin();
             iter!=itsParmNames.end(); ++iter) {
          if (String(*iter).matches (regex)) {
            result.push_back (*iter);
          }
        }
      }
      if (includeDefaults) {
        vector<string> defNames (getDefNames (parmNamePattern));
        result.insert (result.end(), defNames.begin(), defNames.end());
      }
      return result;
    }

    vector<string> ParmFacadeDistr::getDefNames (const string& parmNamePattern) const
    {
      Regex regex(".*");
      if (! parmNamePattern.empty()) {
        regex = Regex(Regex::fromPattern(parmNamePattern));
      }
      vector<string> result;
      result.reserve (itsDefValues.size());
      for (uInt i=0; i<itsDefValues.size(); ++i) {
        const String& name = itsDefValues.name(i);
        if (name.matches (regex)) {
          result.push_back (name);
        }
      }
      return result;
    }

    Record ParmFacadeDistr::getDefValues (const string& parmNamePattern) const
    {
      if (parmNamePattern.empty()  ||  parmNamePattern == "*") {
        return itsDefValues;
      }
      Regex regex(Regex::fromPattern(parmNamePattern));
      Record result;
      for (uInt i=0; i<itsDefValues.size(); ++i) {
        const String& name = itsDefValues.name(i);
        if (name.matches (regex)) {
          result.define (name, itsDefValues.asArrayDouble(i));
        }
      }
      return result;
    }

    void ParmFacadeDistr::addDefValues (const Record& values,
                                        bool check)
    {
      BlobString buf;
      MWBlobOut bbo(buf, AddDefValues, 0);
      putRecord (bbo.blobStream(), values);
      bbo.blobStream() << check;
      bbo.finish();
      itsConn.writeAll (buf);
      checkStatusAll();
    }

    void ParmFacadeDistr::deleteDefValues (const string& parmNamePattern)
    {
      BlobString buf;
      MWBlobOut bbo(buf, DeleteDefValues, 0);
      bbo.blobStream() << parmNamePattern;
      bbo.finish();
      itsConn.writeAll (buf);
      checkStatusAll();
    }

    void ParmFacadeDistr::checkNames (const vector<string>& firstNames,
                                      const vector<string>& names,
                                      uint inx) const
    {
      bool same = (names.size() == firstNames.size());
      uint i=0;
      while (same  &&  i<names.size()) {
        same = (names[i] == firstNames[i]);
        ++i;
      }
      if (!same) {
        LOG_WARN_STR ("names sizes of parts " << itsPartNames[0] << " and "
                      << itsPartNames[inx] << " differ");
      }
    }

    Record ParmFacadeDistr::getValues (const string& parmNamePattern,
                                       double freqv1, double freqv2,
                                       double freqStep,
                                       double timev1, double timev2,
                                       double timeStep,
                                       bool asStartEnd,
                                       bool includeDefaults)
    {
      if (!asStartEnd) {
        double width = freqv2;
        freqv1 -= width/2;
        freqv2 = freqv1 + width;
        width = timev2;
        timev1 -= width/2;
        timev2 = timev1 + width;
      }
      BlobString buf;
      MWBlobOut bbo(buf, GetValues, 0);
      bbo.blobStream() << parmNamePattern
                       << freqv1 << freqv2 << freqStep
                       << timev1 << timev2 << timeStep
                       << includeDefaults;
      bbo.finish();
      itsConn.writeAll (buf);
      // Read the replies.
      vector<Record> recs(itsConn.size());
      for (int i=0; i<itsConn.size(); ++i) {
        itsConn.read (i, buf);
        MWBlobIn bbi(buf);
        if (checkStatus (bbi, i)) {
          getRecord (bbi.blobStream(), recs[i]);
          bbi.finish();
        }
      }
      return combineRemote (recs);
    }

    Record ParmFacadeDistr::getValues (const string& parmNamePattern,
                                       const vector<double>& freqv1,
                                       const vector<double>& freqv2,
                                       const vector<double>& timev1,
                                       const vector<double>& timev2,
                                       bool asStartEnd,
                                       bool includeDefaults)
    {
      BlobString buf;
      MWBlobOut bbo(buf, GetValuesVec, 0);
      bbo.blobStream() << parmNamePattern
                       << freqv1 << freqv2 << timev1 << timev2
                       << asStartEnd << includeDefaults;
      bbo.finish();
      itsConn.writeAll (buf);
      vector<Record> recs(itsConn.size());
      for (int i=0; i<itsConn.size(); ++i) {
        itsConn.read (i, buf);
        MWBlobIn bbi(buf);
        if (checkStatus (bbi, i)) {
          getRecord (bbi.blobStream(), recs[i]);
          bbi.finish();
        }
      }
      return combineRemote (recs);
    }

    Record ParmFacadeDistr::getValuesGrid (const string& parmNamePattern,
                                           double freqv1, double freqv2,
                                           double timev1, double timev2,
                                           bool asStartEnd)
    {
      Box box(freqv1, freqv2, timev1, timev2, asStartEnd);
      BlobString buf;
      MWBlobOut bbo(buf, GetValuesGrid, 0);
      bbo.blobStream() << parmNamePattern
                       << box.lowerX() << box.upperX()
                       << box.lowerY() << box.upperY();
      bbo.finish();
      itsConn.writeAll (buf);
      vector<Record> recs(itsConn.size());
      for (int i=0; i<itsConn.size(); ++i) {
        itsConn.read (i, buf);
        MWBlobIn bbi(buf);
        if (checkStatus (bbi, i)) {
          getRecord (bbi.blobStream(), recs[i]);
          bbi.finish();
        }
      }
      return combineRemote (recs);
    }

    Record ParmFacadeDistr::getCoeff (const string& parmNamePattern,
                                      double freqv1, double freqv2,
                                      double timev1, double timev2,
                                      bool asStartEnd)
    {
      Box box(freqv1, freqv2, timev1, timev2, asStartEnd);
      BlobString buf;
      MWBlobOut bbo(buf, GetCoeff, 0);
      bbo.blobStream() << parmNamePattern
                       << box.lowerX() << box.upperX()
                       << box.lowerY() << box.upperY();
      bbo.finish();
      itsConn.writeAll (buf);
      vector<Record> recs(itsConn.size());
      for (int i=0; i<itsConn.size(); ++i) {
        itsConn.read (i, buf);
        MWBlobIn bbi(buf);
        if (checkStatus (bbi, i)) {
          getRecord (bbi.blobStream(), recs[i]);
          bbi.finish();
        }
      }
      return combineRemote (recs);
    }

    vector<double> ParmFacadeDistr::getDefaultSteps() const
    {
      BlobString buf;
      MWBlobOut bbo(buf, GetDefaultSteps, 0);
      bbo.finish();
      itsConn.writeAll (buf);
      vector<double> steps, defSteps;
      for (int i=0; i<itsConn.size(); ++i) {
        itsConn.read (i, buf);
        MWBlobIn bbi(buf);
        if (checkStatus (bbi, i)) {
          bbi.blobStream() >> steps;
          bbi.finish();
          if (! steps.empty()) {
            if (defSteps.empty()) {
              defSteps = steps;
            } else {
              // The default steps of all parts should match.
              bool ok = true;
              if (steps.size() != defSteps.size()) {
                ok = false;
              } else {
                for (size_t i=0; i<steps.size(); ++i) {
                  if (steps[i] != defSteps[i]) {
                    ok = false;
                    break;
                  }
                }
              }
              if (!ok) {
                LOG_WARN_STR ("Default step values of part " << itsPartNames[i]
                              << " differ from the other parts");
              }
            }
          }
        }
      }
      return defSteps;
    }

    void ParmFacadeDistr::flush (bool fsync)
    {
      BlobString buf;
      MWBlobOut bbo(buf, Flush, 0);
      bbo.blobStream() << fsync;
      bbo.finish();
      itsConn.writeAll (buf);
      checkStatusAll();
    }

    void ParmFacadeDistr::lock (bool lockForWrite)
    {
      BlobString buf;
      MWBlobOut bbo(buf, Lock, 0);
      bbo.blobStream() << lockForWrite;
      bbo.finish();
      itsConn.writeAll (buf);
      checkStatusAll();
    }

    void ParmFacadeDistr::unlock()
    {
      BlobString buf;
      MWBlobOut bbo(buf, Unlock, 0);
      bbo.finish();
      itsConn.writeAll (buf);
      checkStatusAll();
    }

    void ParmFacadeDistr::clearTables()
    {
      BlobString buf;
      MWBlobOut bbo(buf, ClearTables, 0);
      bbo.finish();
      itsConn.writeAll (buf);
      checkStatusAll();
    }

    void ParmFacadeDistr::setDefaultSteps (const vector<double>& steps)
    {
      BlobString buf;
      MWBlobOut bbo(buf, SetDefaultSteps, 0);
      bbo.blobStream() << steps;
      bbo.finish();
      itsConn.writeAll (buf);
      checkStatusAll();
    }

    void ParmFacadeDistr::addValues (const Record&)
      { throw Exception ("ParmFacadeDistr::addValues not implemented"); }

    void ParmFacadeDistr::deleteValues (const string& parmNamePattern,
                                        double freqv1, double freqv2,
                                        double timev1, double timev2,
                                        bool asStartEnd)
    {
      BlobString buf;
      MWBlobOut bbo(buf, DeleteValues, 0);
      bbo.blobStream() << parmNamePattern
                       << freqv1 << freqv2 << timev1 << timev2 << asStartEnd;
      bbo.finish();
      itsConn.writeAll (buf);
      checkStatusAll();
    }


    Record ParmFacadeDistr::combineRemote (const vector<Record>& recs) const
    {
      // Find all possible parm names.
      set<String> names;
      findParmNames (recs, names);
      // Combine the data for each parameter.
      Record result;
      // Step through all parameters.
      for (set<String>::const_iterator iter=names.begin();
           iter!=names.end(); ++iter) {
        combineInfo (*iter, recs, result);
      }
      return result;
    }

    void ParmFacadeDistr::findParmNames (const vector<Record>& recs,
                                         set<String>& names) const
    {
      for (vector<Record>::const_iterator iter=recs.begin();
           iter!=recs.end(); ++iter) {
        for (uint i=0; i<iter->nfields(); ++i) {
          names.insert (iter->name(i));
        }
      }
    }

    void ParmFacadeDistr::combineInfo (const String& name,
                                       const vector<Record>& recs,
                                       Record& result) const
    {
      // Get references to all data arrays.
      vector<Array<double> > values;
      vector<Array<double> > errors;
      vector<Array<double> > fcenters;
      vector<Array<double> > fwidths;
      Array<double> tcenters;
      Array<double> twidths;
      vector<double> freqs;
      values.reserve (recs.size());
      errors.reserve (recs.size());
      fcenters.reserve (recs.size());
      fwidths.reserve (recs.size());
      uint ntime=0;
      bool haveErrors = false;
      for (vector<Record>::const_iterator iter=recs.begin();
           iter!=recs.end(); ++iter) {
        if (iter->isDefined (name)) {
          Record rec = iter->subRecord(name);
          if (values.empty()) {
            // First time.
            tcenters.reference (rec.asArrayDouble("times"));
            twidths.reference  (rec.asArrayDouble("timewidths"));
            ntime = tcenters.size();
            haveErrors = rec.isDefined ("errors");
          } else {
            // Check if the time axis is the same.
            // That should always be the case.
            ASSERT (allNear (tcenters,
                             rec.asArrayDouble("times"), 1e-10));
          }
          values.push_back   (rec.asArrayDouble ("values"));
          if (haveErrors) {
            errors.push_back (rec.asArrayDouble ("errors"));
          }
          fcenters.push_back (rec.asArrayDouble ("freqs"));
          fwidths.push_back  (rec.asArrayDouble ("freqwidths"));
          freqs.push_back (fcenters[freqs.size()].data()[0]);
        }
      }
      // Exit if no matching values found.
      if (values.empty()) {
        return;
      }
      // Now sort (indirectly) the parts in freq order.
      Vector<uInt> indexf;
      GenSortIndirect<double>::sort (indexf, &(freqs[0]), freqs.size());
      // Get the unique frequency domains.
      // This is needed because sometimes BBS solves (partially) globally.
      vector<uint> index;
      index.reserve (indexf.size());
      // The first one always has to be used.
      index.push_back (indexf[0]);
      const Array<double>* lastUsed = &fcenters[indexf[0]];
      uint nfreq = lastUsed->size();
      for (uint i=1; i<indexf.size(); ++i) {
        uint inx = indexf[i];
        const Array<double>& arr = fcenters[inx];
        if (arr.size() != lastUsed->size()  ||
            !allNear (arr, *lastUsed, 1e-10)) {
          index.push_back (inx);
          lastUsed = &fcenters[inx];
          nfreq += lastUsed->size();
        }
      }
      Record rec;
      // Times are the same for all parts, so define them.
      rec.define ("times", tcenters);
      rec.define ("timewidths", twidths);
      // If only one part left, take that one (which is the first one).
      if (index.size() == 1) {
        rec.define ("values", values[0]);
        if (haveErrors) {
          rec.define ("errors", errors[0]);
        }
        rec.define ("freqs", fcenters[0]);
        rec.define ("freqwidths", fwidths[0]);
        result.defineRecord (name, rec);
        return;
      }
      // Combine the values and freq grid.
      // Usually values forms a 2D array, but for funklet coefficients
      // it can be 4D. In that case the first 2 dims are the funklet shape.
      // Get shape from first array.
      IPosition shape = values[index[0]].shape();
      uint ndim = shape.size();
      uint nelem = (ndim <= 2  ?  1 : shape[0]*shape[1]);
      shape[ndim-2] = nfreq;
      shape[ndim-1] = ntime;
      Array<double> data(shape);
      Array<double> errdata;
      if (haveErrors) {
        errdata.resize (shape);
      }
      Array<double> fcenter(IPosition(1, nfreq));
      Array<double> fwidth(IPosition(1, nfreq));
      double* dtp = data.data();
      double* etp = errdata.data();
      double* fcp = fcenter.data();
      double* fwp = fwidth.data();
      // We do the copying ourselves instead of using Array sectioning.
      // It is faster and about as easy to write as we know that all Arrays
      // are contiguous.
      for (uint i=0; i<index.size(); ++i) {
        uint inx = index[i];
        uint nf = fcenters[inx].size();
        memcpy (fcp, fcenters[inx].data(), nf*sizeof(double));
        memcpy (fwp, fwidths[inx].data(),  nf*sizeof(double));
        fcp += nf;
        fwp += nf;
        // Copy values.
        double* to = dtp;
        const double* from = values[inx].data();
        for (uint j=0; j<ntime; ++j) {
          memcpy (to, from, nf*nelem*sizeof(double));
          to += nfreq*nelem;
          from += nf*nelem;
        }
        dtp += nf*nelem;
        if (haveErrors) {
          // Copy errors.
          double* to = etp;
          const double* from = errors[inx].data();
          for (uint j=0; j<ntime; ++j) {
            memcpy (to, from, nf*nelem*sizeof(double));
            to += nfreq*nelem;
            from += nf*nelem;
          }
          etp += nf*nelem;
        }
      }
      rec.define ("values", data);
      if (haveErrors) {
        rec.define ("errors", errdata);
      }
      rec.define ("freqs", fcenter);
      rec.define ("freqwidths", fwidth);
      result.defineRecord (name, rec);
    }

    void ParmFacadeDistr::getRecord (BlobIStream& bis, Record& rec)
    {
      BlobAipsIO baio(bis);
      casacore::AipsIO aio(&baio);
      aio >> rec;
    }

    void ParmFacadeDistr::putRecord (BlobOStream& bos, const Record& rec)
    {
      BlobAipsIO baio(bos);
      casacore::AipsIO aio(&baio);
      aio << rec;
    }

  } // namespace ParmDB
} // namespace LOFAR
