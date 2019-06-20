//# PreFlagger.cc: DPPP step class to (un)flag data on channel, baseline, time
//# Copyright (C) 2010
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
//# $Id$
//#
//# @author Ger van Diepen

#include "PreFlagger.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Exceptions.h"

#include "../Common/ParameterSet.h"
#include "../Common/StreamUtil.h"

#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCEpoch.h>

#include <boost/algorithm/string/case_conv.hpp>

#include <iostream>
#include <stack>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    PreFlagger::PreFlagger (DPInput* input,
                            const ParameterSet& parset,
                            const string& prefix)
      : itsName        (prefix),
        itsInput       (input),
        itsMode        (SetFlag),
        itsPSet        (input, parset, prefix),
        itsCount       (0),
        itsFlagCounter (input->msName(), parset, prefix+"count.")
    {
      string mode = boost::to_lower_copy(parset.getString(prefix+"mode", "set"));
      if (mode == "clear") {
        itsMode = ClearFlag;
      } else if (mode == "setcomplement"  ||  mode == "setother") {
        itsMode = SetComp;
      } else if (mode == "clearcomplement"  ||  mode == "clearother") {
        itsMode = ClearComp;
      } else {
        if (mode!="set")
          throw Exception(
                   "invalid preflagger mode: "
                   "only set, clear, and set/clearcomplement/other are valid");
      }
    }

    PreFlagger::~PreFlagger()
    {}

    void PreFlagger::show (std::ostream& os) const
    {
      os << "PreFlagger " << itsName << std::endl;
      os << "  mode:           ";
      switch (itsMode) {
      case SetFlag:
        os << "set";
        break;
      case ClearFlag:
        os << "clear";
        break;
      case SetComp:
        os << "setcomplement";
        break;
      case ClearComp:
        os << "clearcomplement";
        break;
      }
      os << endl;
      itsPSet.show (os, false);
    }

    void PreFlagger::showCounts (std::ostream& os) const
    {
      os << endl << "Flags set by PreFlagger " << itsName;
      os << endl << "=======================" << endl;
      itsFlagCounter.showBaseline (os, itsCount);
      itsFlagCounter.showChannel  (os, itsCount);
    }

    void PreFlagger::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " PreFlagger " << itsName << endl;
    }

    void PreFlagger::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteFlags();
      itsPSet.updateInfo (getInfo());
      // Initialize the flag counters.
      itsFlagCounter.init (getInfo());
    }

    bool PreFlagger::process (const DPBuffer& buf)
    {
      itsTimer.start();
      // Because no buffers are kept, we can reference the filled arrays
      // in the input buffer instead of copying them.
      itsBuffer.referenceFilled (buf);
      // Do the PSet steps and combine the result with the current flags.
      // Only count if the flag changes.
      Cube<bool>* flags = itsPSet.process (buf, itsBuffer, itsCount,
                                           Block<bool>(), itsTimer);
      const IPosition& shape = flags->shape();
      unsigned int nrcorr = shape[0];
      unsigned int nrchan = shape[1];
      unsigned int nrbl   = shape[2];
      const bool* inPtr = flags->data();
      bool* outPtr = itsBuffer.getFlags().data();
      switch (itsMode) {
      case SetFlag:
        setFlags (inPtr, outPtr, nrcorr, nrchan, nrbl, true);
        break;
      case ClearFlag:
        clearFlags (inPtr, outPtr, nrcorr, nrchan, nrbl, true, buf);
        break;
      case SetComp:
        setFlags (inPtr, outPtr, nrcorr, nrchan, nrbl, false);
        break;
      case ClearComp:
        clearFlags (inPtr, outPtr, nrcorr, nrchan, nrbl, false, buf);
        break;
      }
      itsTimer.stop();
      // Let the next step do its processing.
      getNextStep()->process (itsBuffer);
      itsCount++;
      return true;
    }

    void PreFlagger::setFlags (const bool* inPtr, bool* outPtr,
                               unsigned int nrcorr, unsigned int nrchan, unsigned int nrbl,
                               bool mode)
    {
      for (unsigned int i=0; i<nrbl; ++i) {
        for (unsigned int j=0; j<nrchan; ++j) {
          if (*inPtr == mode  &&  !*outPtr) {
            // Only 1st corr is counted.
            itsFlagCounter.incrBaseline(i);
            itsFlagCounter.incrChannel(j);
            for (unsigned int k=0; k<nrcorr; ++k) {
              outPtr[k] = true;
            }
          }
          inPtr  += nrcorr;
          outPtr += nrcorr;
        }
      }
    }

    void PreFlagger::clearFlags (const bool* inPtr, bool* outPtr,
                                 unsigned int nrcorr, unsigned int nrchan, unsigned int nrbl,
                                 bool mode, const DPBuffer& buf)
    {
      const Complex* dataPtr = buf.getData().data();
      Cube<float> weights = itsInput->fetchWeights (buf, itsBuffer,
                                                    itsTimer);
      const float* weightPtr = weights.data();
      for (unsigned int i=0; i<nrbl; ++i) {
        for (unsigned int j=0; j<nrchan; ++j) {
          if (*inPtr == mode) {
            bool flag = false;
            // Flags for invalid data are not cleared.
            for (unsigned int k=0; k<nrcorr; ++k) {
              if (!isFinite(dataPtr[k].real())  ||
                  !isFinite(dataPtr[k].imag())  ||
                  weightPtr[k] == 0) {
                flag = true;
                break;
              }
            }
            if (*outPtr != flag) {
              itsFlagCounter.incrBaseline(i);
              itsFlagCounter.incrChannel(j);
              for (unsigned int k=0; k<nrcorr; ++k) {
                outPtr[k] = flag;
              }
            }
          }
          inPtr     += nrcorr;
          outPtr    += nrcorr;
          dataPtr   += nrcorr;
          weightPtr += nrcorr;
        }
      }
    }

    void PreFlagger::finish()
    {
      // Let the next step finish its processing.
      getNextStep()->finish();
    }


    PreFlagger::PSet::PSet (DPInput* input,
                            const ParameterSet& parset,
                            const string& prefix)
      : itsInput (input),
        itsName  (prefix),
        itsFlagOnUV    (false),
        itsFlagOnBL    (false),
        itsFlagOnAmpl  (false),
        itsFlagOnPhase (false),
        itsFlagOnReal  (false),
        itsFlagOnImag  (false),
        itsFlagOnAzEl  (false),
        itsSelBL       (parset, prefix, true)
    {
      // Read all possible parameters.
      itsStrTime  = parset.getStringVector (prefix+"timeofday",
                                            vector<string>());
      itsStrLST   = parset.getStringVector (prefix+"lst",
                                            vector<string>());
      itsStrATime = parset.getStringVector (prefix+"abstime",
                                            vector<string>());
      itsStrRTime = parset.getStringVector (prefix+"reltime",
                                            vector<string>());
      itsTimeSlot = parset.getUintVector   (prefix+"timeslot",
                                            vector<unsigned int>(), true); // expand ..
      itsStrAzim  = parset.getStringVector (prefix+"azimuth",
                                            vector<string>());
      itsStrElev  = parset.getStringVector (prefix+"elevation",
                                            vector<string>());
      itsMinUV    = parset.getDouble       (prefix+"uvmmin", -1);
      itsMaxUV    = parset.getDouble       (prefix+"uvmmax", -1);
      itsStrFreq  = parset.getStringVector (prefix+"freqrange",
                                            vector<string>());
      itsStrChan  = parset.getStringVector (prefix+"chan", vector<string>());
      itsAmplMin = fillValuePerCorr
        (ParameterValue (parset.getString  (prefix+"amplmin", string())),-1e30,
         itsFlagOnAmpl);
      itsAmplMax  = fillValuePerCorr
        (ParameterValue (parset.getString  (prefix+"amplmax", string())), 1e30,
         itsFlagOnAmpl);
      itsPhaseMin = fillValuePerCorr
        (ParameterValue (parset.getString  (prefix+"phasemin", string())),-1e30,
         itsFlagOnPhase);
      itsPhaseMax  = fillValuePerCorr
        (ParameterValue (parset.getString  (prefix+"phasemax", string())), 1e30,
         itsFlagOnPhase);
      itsRealMin = fillValuePerCorr
        (ParameterValue (parset.getString  (prefix+"realmin", string())),-1e30,
         itsFlagOnReal);
      itsRealMax  = fillValuePerCorr
        (ParameterValue (parset.getString  (prefix+"realmax", string())), 1e30,
         itsFlagOnReal);
      itsImagMin = fillValuePerCorr
        (ParameterValue (parset.getString  (prefix+"imagmin", string())),-1e30,
         itsFlagOnImag);
      itsImagMax  = fillValuePerCorr
        (ParameterValue (parset.getString  (prefix+"imagmax", string())), 1e30,
         itsFlagOnImag);
      itsStrExpr = parset.getString   (prefix+"expr", string());
      // Parse the possible pset expression and convert to RPN form.
      if (! itsStrExpr.empty()) {
        vector<string> names = exprToRpn (itsStrExpr);
        // Create PSet objects for all operands.
        itsPSets.reserve (names.size());
        for (unsigned int i=0; i<names.size(); ++i) {
          itsPSets.push_back
            (PSet::ShPtr(new PSet(itsInput, parset, prefix+names[i]+'.')));
        }
      }
    }

    void PreFlagger::PSet::updateInfo (const DPInfo& info)
    {
      itsInfo = &info;
      // Fill the matrix with the baselines to flag.
      fillBLMatrix();
      // Handle the possible date/time parameters.
      itsTimes  = fillTimes (itsStrTime,  true,  true);
      itsLST    = fillTimes (itsStrLST,   true,  true);
      itsATimes = fillTimes (itsStrATime, false, false);
      itsRTimes = fillTimes (itsStrRTime, true,  false);
      itsFlagOnTime = !(itsTimeSlot.empty() && itsTimes.empty() &&
                        itsLST.empty() && itsATimes.empty() &&
                        itsRTimes.empty());
      // Handle possible azimuth/elevation ranges.
      itsAzimuth    = fillTimes (itsStrAzim, true, true);
      itsElevation  = fillTimes (itsStrElev, true, true);
      itsFlagOnAzEl = !(itsAzimuth.empty() && itsElevation.empty());
      // Determine if to flag on UV distance.
      // If so, square the distances to avoid having to take the sqrt in flagUV.
      if (itsMinUV >= 0) {
        itsFlagOnUV = true;
        itsMinUV   *= itsMinUV;
      }
      if (itsMaxUV > 0) {
        itsFlagOnUV = true;
        itsMaxUV   *= itsMaxUV;
      } else {
        // Make it a very high number.
        itsMaxUV = 1e30;
      }
      if (itsMinUV>=itsMaxUV)
        throw Exception("PreFlagger uvmmin should be < uvmmax");
      // Determine if only flagging on time info is done.
      itsFlagOnTimeOnly = ( !(itsFlagOnUV || itsFlagOnBL || itsFlagOnAzEl ||
                              itsFlagOnAmpl || itsFlagOnPhase ||
                              itsFlagOnReal || itsFlagOnImag) &&
                            itsPSets.empty());
      // Size the object's buffers (used in process) correctly.
      unsigned int nrcorr = info.ncorr();
      unsigned int nrchan = info.nchan();
      itsFlags.resize (nrcorr, nrchan, info.nbaselines());
      itsMatchBL.resize (info.nbaselines());
      // Determine the channels to be flagged.
      if (!(itsStrChan.empty() && itsStrFreq.empty())) {
        fillChannels (info);
        if (! itsChannels.empty()) {
          itsFlagOnTimeOnly = false;
        }
      }
      // Do the same for the child steps.
      for (unsigned int i=0; i<itsPSets.size(); ++i) {
        itsPSets[i]->updateInfo (info);
      }
    }

    void PreFlagger::PSet::fillChannels (const DPInfo& info)
    {
      unsigned int nrcorr = info.ncorr();
      unsigned int nrchan = info.nchan();
      Vector<bool> selChan(nrchan);
      if (itsStrChan.empty()) {
        selChan = true;
      } else {
        // Set selChan for channels not exceeding nr of channels.
        selChan = false;
        Record rec;
        rec.define ("nchan", nrchan);
        double result;
        for (unsigned int i=0; i<itsStrChan.size(); ++i) {
          // Evaluate possible expressions.
          // Split the value if start..end is given.
          unsigned int startch, endch;
          string::size_type pos = itsStrChan[i].find ("..");
          if (pos == string::npos) {
            TableExprNode node (RecordGram::parse(rec, itsStrChan[i]));
            node.get (rec, result);
            startch = (unsigned int)(result+0.001);
            endch   = startch;
          } else {
            if (pos == 0 || pos >= itsStrChan[i].size() - 2)
              throw Exception(
                       "No start or end given in PreFlagger channel range "
                       + itsStrChan[i]);
            TableExprNode node1
              (RecordGram::parse(rec, itsStrChan[i].substr(0,pos)));
            node1.get (rec, result);
            startch = (unsigned int)(result+0.001);
            TableExprNode node2
              (RecordGram::parse(rec, itsStrChan[i].substr(pos+2)));
            node2.get (rec, result);
            endch = (unsigned int)(result+0.001);
            if(startch > endch)
              throw Exception(
                       "Start " + std::to_string(startch) + " must be <= end " + std::to_string(endch)
                       + " in PreFlagger channel range " + itsStrChan[i]);
          }
          if (startch < nrchan) {
            for (unsigned int ch=startch; ch<std::min(endch+1, nrchan); ++ch) {
              selChan[ch] = true;
            }
          }
        }
      }
      // Now determine which channels to use from given frequency ranges.
      // AND it with the channel selection given above.
      if (! itsStrFreq.empty()) {
        selChan = selChan && handleFreqRanges (itsInfo->chanFreqs());
      }
      // Turn the channels into a mask.
      itsChannels.clear();
      itsChanFlags.resize (nrcorr, nrchan);
      itsChanFlags = false;
      for (unsigned int i=0; i<nrchan; ++i) {
        if (selChan[i]) {
          itsChannels.push_back (i);
          for (unsigned int j=0; j<nrcorr; ++j) {
            itsChanFlags(j,i) = true;
          }
        }
      }
    }

    void PreFlagger::PSet::show (std::ostream& os, bool showName) const
    {
      if (showName) {
        os << "  pset " << itsName << endl;
      }
      if (! itsStrExpr.empty()) {
        os << "   expr:          " << itsStrExpr << std::endl;
      }
      if (! itsStrLST.empty()) {
        os << "   lst:           " << itsStrLST << std::endl;
      }
      if (! itsStrTime.empty()) {
        os << "   timeofday:     " << itsStrTime << std::endl;
      }
      if (! itsStrATime.empty()) {
        os << "   abstime:       " << itsStrATime << std::endl;
      }
      if (! itsStrRTime.empty()) {
        os << "   reltime:       " << itsStrRTime << std::endl;
      }
      if (! itsTimeSlot.empty()) {
        os << "   timeslot:      " << itsTimeSlot << std::endl;
      }
      if (itsFlagOnBL) {
        itsSelBL.show (os);
      }
      if (itsFlagOnUV) {
        if (itsMinUV >= 0) {
          os << "   uvmmin:        " << sqrt(itsMinUV) << std::endl;
        } else {
          os << "   uvmmin:        " << itsMinUV << std::endl;
        }
        os << "   uvmmax:        " << sqrt(itsMaxUV) << std::endl;
      }
      if (itsFlagOnAzEl) {
        os << "   azimuth:       " << itsStrAzim << std::endl;
        os << "   elevation:     " << itsStrElev << std::endl;
      }
      if (! itsChannels.empty()) {
        os << "   channel:       " << itsStrChan << std::endl;
        os << "   freqrange:     " << itsStrFreq << std::endl;
        os << "    chan to flag: " << itsChannels << std::endl;
      }
      if (itsFlagOnAmpl) {
        os << "   amplmin:       " << itsAmplMin << std::endl;
        os << "   amplmax:       " << itsAmplMax << std::endl;
      }
      if (itsFlagOnPhase) {
        os << "   phasemin:      " << itsPhaseMin << std::endl;
        os << "   phasemax:      " << itsPhaseMax << std::endl;
      }
      if (itsFlagOnReal) {
        os << "   realmin:       " << itsRealMin << std::endl;
        os << "   realmax:       " << itsRealMax << std::endl;
      }
      if (itsFlagOnImag) {
        os << "   imagmin:       " << itsImagMin << std::endl;
        os << "   imagmax:       " << itsImagMax << std::endl;
      }
      // Do it for the child steps.
      for (unsigned int i=0; i<itsPSets.size(); ++i) {
        itsPSets[i]->show (os, true);
      }
    }

    Cube<bool>* PreFlagger::PSet::process (const DPBuffer& in,
                                           DPBuffer& out,
                                           unsigned int timeSlot,
                                           const Block<bool>& matchBL,
                                           NSTimer& timer)
    {
      // No need to process it if the time mismatches or if only time selection.
      if (itsFlagOnTime) {
        if (! matchTime (out.getTime(), timeSlot)) {
          itsFlags = false;
          return &itsFlags;
        }
      }
      if (itsFlagOnTimeOnly) {
        itsFlags = itsFlagOnTime;
        return &itsFlags;
      }
      // Initialize the flags.
      itsFlags = false;
      const IPosition& shape = out.getFlags().shape();
      unsigned int nr = shape[0] * shape[1];
      // Take over the baseline info from the parent. Default is all.
      if (matchBL.empty()) {
        itsMatchBL = true;
      } else {
        itsMatchBL = matchBL;
      }
      // The PSet tree is a combination of ORs and ANDs.
      // Depth is AND, breadth is OR.
      // In a PSet flagging is done in two stages.
      // First it is determined which baselines are not flagged. It is kept
      // in the itsMatchBL block.
      // This is passed to the children who do their flagging. In this way
      // a child can minimize the amount of work to do.
      // The resulting flags of the children are handled according to the
      // operators in the RPN list.

      // First flag on baseline if necessary. Stop if no matches. 
      if (itsFlagOnBL  &&  !flagBL()) {
        return &itsFlags;
      }
      // Flag on UV distance if necessary.
      if (itsFlagOnUV  &&  !flagUV (itsInput->fetchUVW (in, out,
                                                        timer))) {
        return &itsFlags;
      }
      // Flag on AzEl is necessary.
      if (itsFlagOnAzEl  &&  !flagAzEl (out.getTime())) {
        return &itsFlags;
      }
      // Convert each baseline flag to a flag per correlation/channel.
      bool* flagPtr = itsFlags.data();
      for (unsigned int i=0; i<itsMatchBL.size(); ++i) {
        if (itsMatchBL[i]) {
          std::fill (flagPtr, flagPtr+nr, itsMatchBL[i]);
        }
        flagPtr += nr;
      }
      // Flag on channel if necessary.
      if (! itsChannels.empty()) {
        flagChannels();
      }
      // Flag on amplitude, phase or real/imaginary if necessary.
      if (itsFlagOnAmpl) {
        flagAmpl (amplitude(out.getData()));
      }
      if (itsFlagOnReal) {
        flagReal (out.getData());
      }
      if (itsFlagOnImag) {
        flagImag (out.getData());
      }
      if (itsFlagOnPhase) {
        flagPhase (out.getData());
      }
      // Evaluate the PSet expression.
      // The expression is in RPN notation. A stack of array pointers is used
      // to keep track of intermediate results. The arrays (in the PSet objects)
      // are reused to AND or OR subexpressions. This can be done harmlessly
      // and saves the creation of too many arrays.
      if (! itsPSets.empty()) {
        std::stack<Cube<bool>*> results;
        for (int oper : itsRpn) {
          if (oper >= 0) {
            results.push (itsPSets[oper]->process (in, out, timeSlot,
                                                    itsMatchBL, timer));
          } else if (oper == OpNot) {
            Cube<bool>* left = results.top();
            // No ||= operator exists, so use the transform function.
            transformInPlace (left->cbegin(), left->cend(),
                              std::logical_not<bool>());
          } else if(oper == OpOr || oper == OpAnd) {
            Cube<bool>* right = results.top();
            results.pop();
            Cube<bool>* left = results.top();
            // No ||= operator exists, so use the transform function.
            if (oper == OpOr) {
              transformInPlace (left->cbegin(), left->cend(),
                                right->cbegin(), std::logical_or<bool>());
            } else {
              transformInPlace (left->cbegin(), left->cend(),
                                right->cbegin(), std::logical_and<bool>());
            }
          } else {
            throw std::runtime_error("Expected operation NOT, OR or AND");
          }
        }
        // Finally AND the children's flags with the flags of this pset.
        if (results.size() != 1)
          throw std::runtime_error("Something went wrong while evaluating expression: results.size() != 1");
        Cube<bool>* mflags = results.top();
        transformInPlace (itsFlags.cbegin(), itsFlags.cend(),
                          mflags->cbegin(), std::logical_and<bool>());
      }
      return &itsFlags;
    }

    bool PreFlagger::PSet::matchTime (double time, unsigned int timeSlot) const
    {
      if (!itsATimes.empty()  &&
          !matchRange (time, itsATimes)) {
        return false;
      }
      if (!itsRTimes.empty()  &&
          !matchRange (time-itsInfo->startTime(), itsRTimes)) {
        return false;
      }
      if (!itsTimes.empty()) {
        MVTime mvtime(time/86400);    // needs time in days
        double timeofday = time - int(mvtime.day()) * 86400.;
        if (!matchRange (timeofday, itsTimes)) {
          return false;
        }
      }
      if (!itsTimeSlot.empty()) {
        if (std::find (itsTimeSlot.begin(), itsTimeSlot.end(), timeSlot) ==
            itsTimeSlot.end()) {
          return false;
        }
      }
      if (!itsLST.empty()) {
        // Convert time from UTC to Local Apparent Sidereal Time.
        MeasFrame frame;
        frame.set (itsInfo->arrayPos());
        Quantity qtime(time, "s");
        MEpoch lst = MEpoch::Convert (MEpoch(MVEpoch(qtime), MEpoch::UTC),
                                      MEpoch::Ref(MEpoch::LAST, frame))();
        double lstSec = lst.getValue().get();        // in days
        lstSec -= int(lstSec);                       // time of day
        if (!matchRange (lstSec*86400, itsLST)) {    // use seconds
          return false;
        }
      }
      return true;
    }

    bool PreFlagger::PSet::matchRange (double v,
                                       const vector<double>& ranges) const
    {
      for (unsigned int i=0; i<ranges.size(); i+=2) {
        if (v > ranges[i]  &&  v < ranges[i+1]) {
          return true;
        }
      }
      return false;
    }

    bool PreFlagger::PSet::flagUV (const Matrix<double>& uvw)
    {
      bool match = false;
      unsigned int nrbl = itsMatchBL.size();
      const double* uvwPtr = uvw.data();
      for (unsigned int i=0; i<nrbl; ++i) {
        if (itsMatchBL[i]) {
          // UV-distance is sqrt(u^2 + v^2).
          // The sqrt is not needed because minuv and maxuv are squared.
          double uvdist = uvwPtr[0] * uvwPtr[0] + uvwPtr[1] * uvwPtr[1];
          if (uvdist >= itsMinUV  &&  uvdist <= itsMaxUV) {
            // UV-dist mismatches, so do not flag baseline.
            itsMatchBL[i] = false;
          } else {
            match = true;
          }
        }
        uvwPtr += 3;
      }
      return match;
    }

    bool PreFlagger::PSet::flagBL()
    {
      bool match = false;
      unsigned int nrbl = itsMatchBL.size();
      const Int* ant1Ptr = itsInfo->getAnt1().data();
      const Int* ant2Ptr = itsInfo->getAnt2().data();
      for (unsigned int i=0; i<nrbl; ++i) {
        if (itsMatchBL[i]) {
          if (! itsFlagBL(ant1Ptr[i], ant2Ptr[i])) {
            // do not flag this baseline
            itsMatchBL[i] = false;
          } else {
            match = true;
          }
        }
      }
      return match;
    }

    bool PreFlagger::PSet::flagAzEl (double time)
    {
      bool match = false;
      unsigned int nrbl = itsMatchBL.size();
      const Int* ant1Ptr = itsInfo->getAnt1().data();
      const Int* ant2Ptr = itsInfo->getAnt2().data();
      // Calculate AzEl for each flagged antenna for this time slot.
      MeasFrame frame;
      Quantity qtime(time, "s");
      MEpoch epoch(MVEpoch(qtime), MEpoch::UTC);
      frame.set (epoch);
      MDirection::Convert converter (itsInfo->phaseCenter(),
                                     MDirection::Ref(MDirection::AZEL, frame));
      unsigned int nrant = itsInfo->antennaNames().size();
      Block<bool> done(nrant, false);
      for (unsigned int i=0; i<nrbl; ++i) {
        if (itsMatchBL[i]) {
          // If needed, check if ant1 matches AzEl criterium.
          // If not matching, itsMatchBL is cleared for this baseline and all
          // subsequent baselines with this antenna.
          int a1 = ant1Ptr[i];
          int a2 = ant2Ptr[i];
          if (!done[a1]) {
            frame.set (itsInfo->antennaPos()[a1]);
            testAzEl (converter, i, a1, ant1Ptr, ant2Ptr);
            done[a1]= true;
          }
          // If needed, check if ant2 matches AzEl criterium.
          if (itsMatchBL[i]  &&  !done[a2]) {
            frame.set (itsInfo->antennaPos()[a2]);
            testAzEl (converter, i, a2, ant1Ptr, ant2Ptr);
            done[a2] = true;
          }
          if (itsMatchBL[i]) {
            match = true;
          }
        }
      }
      return match;
    }

    void PreFlagger::PSet::testAzEl (MDirection::Convert& converter,
                                     unsigned int blnr, int ant,
                                     const int* ant1, const int* ant2)
    {
      // Calculate AzEl (n seconds because ranges are in seconds too).
      MVDirection mvAzel (converter().getValue());
      Vector<double> azel = mvAzel.getAngle("s").getValue();
      double az = azel[0];
      double el = azel[1];
      if (az < 0) az += 86400;
      if (el < 0) el += 86400;
      // If outside the ranges, there is no match.
      // Set no match for all baselines containing this antenna.
      // It needs to be done from this baseline on, because the earlier
      // baselines have already been handled.
      bool res = ((itsAzimuth.empty()   || matchRange(az, itsAzimuth)) &&
                  (itsElevation.empty() || matchRange(el, itsElevation)));
      if (!res) {
        for (unsigned int i=blnr; i<itsMatchBL.size(); ++i) {
          if (ant1[i] == ant  ||  ant2[i] == ant) {
            itsMatchBL[i] = false;
          }
        }
      }
    }

    void PreFlagger::PSet::flagAmpl (const Cube<float>& values)
    {
      const IPosition& shape = values.shape();
      unsigned int nrcorr = shape[0];
      unsigned int nr = shape[1] * shape[2];
      const float* valPtr = values.data();
      bool* flagPtr = itsFlags.data();
      for (unsigned int i=0; i<nr; ++i) {
        bool flag = false;
        for (unsigned int j=0; j<nrcorr; ++j) {
          if (valPtr[j] < itsAmplMin[j]  ||  valPtr[j] > itsAmplMax[j]) {
            flag = true;
            break;
          }
        }
        if (!flag) {
          for (unsigned int j=0; j<nrcorr; ++j) {
            flagPtr[j] = false;
          }
        }
        valPtr  += nrcorr;
        flagPtr += nrcorr;
      }
    }

    void PreFlagger::PSet::flagPhase (const Cube<Complex>& values)
    {
      const IPosition& shape = values.shape();
      unsigned int nrcorr = shape[0];
      unsigned int nr = shape[1] * shape[2];
      const Complex* valPtr = values.data();
      bool* flagPtr = itsFlags.data();
      for (unsigned int i=0; i<nr; ++i) {
        bool flag = false;
        for (unsigned int j=0; j<nrcorr; ++j) {
          float phase = arg(valPtr[j]);
          if (phase < itsPhaseMin[j]  ||  phase > itsPhaseMax[j]) {
            flag = true;
            break;
          }
        }
        if (!flag) {
          for (unsigned int j=0; j<nrcorr; ++j) {
            flagPtr[j] = false;
          }
        }
        valPtr  += nrcorr;
        flagPtr += nrcorr;
      }
    }

    void PreFlagger::PSet::flagReal (const Cube<Complex>& values)
    {
      const IPosition& shape = values.shape();
      unsigned int nrcorr = shape[0];
      unsigned int nr = shape[1] * shape[2];
      const Complex* valPtr = values.data();
      bool* flagPtr = itsFlags.data();
      for (unsigned int i=0; i<nr; ++i) {
        bool flag = false;
        for (unsigned int j=0; j<nrcorr; ++j) {
          if (valPtr[j].real() < itsRealMin[j]  ||
               valPtr[j].real() > itsRealMax[j]) {
            flag = true;
            break;
          }
        }
        if (!flag) {
          for (unsigned int j=0; j<nrcorr; ++j) {
            flagPtr[j] = false;
          }
        }
        valPtr  += nrcorr;
        flagPtr += nrcorr;
      }
    }

    void PreFlagger::PSet::flagImag (const Cube<Complex>& values)
    {
      const IPosition& shape = values.shape();
      unsigned int nrcorr = shape[0];
      unsigned int nr = shape[1] * shape[2];
      const Complex* valPtr = values.data();
      bool* flagPtr = itsFlags.data();
      for (unsigned int i=0; i<nr; ++i) {
        bool flag = false;
        for (unsigned int j=0; j<nrcorr; ++j) {
          if (valPtr[j].imag() < itsImagMin[j]  ||
              valPtr[j].imag() > itsImagMax[j]) {
            flag = true;
            break;
          }
        }
        if (!flag) {
          for (unsigned int j=0; j<nrcorr; ++j) {
            flagPtr[j] = false;
          }
        }
        valPtr  += nrcorr;
        flagPtr += nrcorr;
      }
    }

    void PreFlagger::PSet::flagChannels()
    {
      const IPosition& shape = itsFlags.shape();
      unsigned int nr   = shape[0] * shape[1];
      unsigned int nrbl = shape[2];
      bool* flagPtr = itsFlags.data();
      for (unsigned int i=0; i<nrbl; ++i) {
        transformInPlace (flagPtr, flagPtr+nr,
                          itsChanFlags.cbegin(), std::logical_and<bool>());
        flagPtr += nr;
      }
    }

    // See http://montcs.bloomu.edu/~bobmon/Information/RPN/infix2rpn.shtml
    // for the algorithm used here.
    // Some code was added to check if no two subsequent operators or names
    // are given.
    vector<string> PreFlagger::PSet::exprToRpn (const string& origExpr)
    {
      // Operators & (or &&) | (or ||) and , are used as well as parentheses.
      // The operators must have a value in order of precedence, thus &&
      // has a higher precedence than || (as in C).
      string expr = boost::to_upper_copy(origExpr);
      std::stack<int> tokens;
      vector<string> names;
      unsigned int i=0;
      bool hadName = false;    // the last token was a name.
      while (i < expr.size()) {
        int oper = 0;
        // skip whitespace.
        // Look for parenthesis or operator.
        if (expr[i] == ' '  ||  expr[i] == '\t') {
          i++;
        } else if (expr[i] == '(') {
          if(hadName)
            throw Exception("no operator before opening parenthesis at pos. "
                     + std::to_string(i) + " in expression: " + origExpr);
          oper = OpParen;
          tokens.push (oper);
          i++;
        } else if (expr[i] == '|') {
          oper = OpOr;
          i++;
          if (i < expr.size()  &&  expr[i] == '|') i++;
        } else if (expr[i] == ',') {
          oper = OpOr;
          i++;
        } else if (expr[i] == '&') {
          oper = OpAnd;
          i++;
          if (i < expr.size()  &&  expr[i] == '&') i++;
        } else if (expr[i] == '!') {
          oper = OpNot;
          i++;
        } else if (expr.size()-i >= 3  &&  (expr.substr(i,3) == "OR "  ||
                                            expr.substr(i,3) == "OR\t" ||
                                            expr.substr(i,3) == "OR!"  ||
                                            expr.substr(i,3) == "OR(")) {
          oper = OpOr;
          i+=2;
        } else if (expr.size()-i == 2  &&  (expr.substr(i,2) == "OR")) {
          oper = OpOr;
          i+=2;
        } else if (expr.size()-i >= 4  &&  (expr.substr(i,4) == "AND "  ||
                                            expr.substr(i,4) == "AND\t" ||
                                            expr.substr(i,4) == "AND!"  ||
                                            expr.substr(i,4) == "AND(")) {
          oper = OpAnd;
          i+=3;
        } else if (expr.size()-i == 3  &&  (expr.substr(i,3) == "AND")) {
          oper = OpAnd;
          i+=3;
        } else if (expr.size()-i >= 4  &&  (expr.substr(i,4) == "NOT "  ||
                                            expr.substr(i,4) == "NOT\t" ||
                                            expr.substr(i,4) == "NOT(")) {
          oper = OpNot;
          i+=3;
        } else if (expr.size()-i == 3  &&  (expr.substr(i,3) == "NOT")) {
          oper = OpNot;
          i+=3;
        } else if (expr[i] == ')') {
          // Closing parenthesis. Push till opening parenthesis found.
          if (!hadName)
            throw Exception("no set name before closing parenthesis at pos. "
                     + std::to_string(i) + " in expression: " + origExpr);
          while (true) {
            if(tokens.empty())
              throw Exception("mismatched parentheses at pos. "
                       + std::to_string(i) + " in expression: " + origExpr);
            if (tokens.top() == OpParen) {
              tokens.pop();
              break;
            }
            itsRpn.push_back (tokens.top());
            tokens.pop();
          }
          i++;
        } else {
          // No operator, thus it must be an operand (a set name).
          int st=i;
          if (hadName)
            throw Exception("No operator between set names at pos. "
                     + std::to_string(i) + " in expression: " + origExpr);
          while (i < expr.size() && 
                 expr[i] != ' ' && expr[i] != '\t' &&
                 expr[i] != '(' && expr[i] != ')' && expr[i] != '!' &&
                 expr[i] != ',' && expr[i] != '&' && expr[i] != '|') {
            i++;
          }
          hadName = true;
          itsRpn.push_back (names.size());
          String setName (origExpr.substr(st, i-st));
          // Check the name is valid (no special characters).
          if (!setName.matches (RXidentifier))
            throw Exception(
                     "Invalid set name " + std::string(setName)
                     + " used in set expression " + origExpr);
          names.push_back (setName);
        }
        if (oper < OpParen) {
          // Check if an operator was preceeded correctly.
          if (oper == OpNot) {
            if (hadName)
              throw Exception("No set name before operator ! at pos. "
                       + std::to_string(i) + " in expression: " + origExpr);
          } else {
            if (!hadName)
              throw Exception("No set name before operator at pos. "
                       + std::to_string(i) + " in expression: " + origExpr);
          }
          hadName = false;
          // Push till lower precedence found.
          while (!tokens.empty()  &&  tokens.top() < oper) {
            itsRpn.push_back (tokens.top());
            tokens.pop();
          }
          tokens.push (oper);
        }
      }
      if (!hadName)
        throw Exception("no set name after last operator in expression: "
                 + origExpr);
      while (!tokens.empty()) {
        if (tokens.top()>=OpParen)
          throw Exception(
                   "mismatched parentheses in expression: " + origExpr);
        itsRpn.push_back (tokens.top());
        tokens.pop();
      }
      return names;
    }

    vector<double> PreFlagger::PSet::fillTimes (const vector<string>& vec,
                                                bool asTime,
                                                bool canEndBeforeStart)
    {
      vector<double> result;
      result.reserve (2*vec.size());
      // A time range can be given as time..time or time+-value.
      for (vector<string>::const_iterator str = vec.begin();
           str != vec.end(); ++str) {
        // Find the .. or +- token.
        bool usepm = false;
        string::size_type pos = str->find ("..");
        if (pos == string::npos) {
          usepm = true;
          pos = str->find ("+-");
          if (pos == string::npos)
            throw Exception("PreFlagger time range '" + *str
                     + "' should be range using .. or +-");
        }
        // Get the time or datetime in seconds. The values must be positive.
        double v1 = getSeconds (str->substr(0, pos), asTime, false);
        double v2 = getSeconds (str->substr(pos+2),  asTime, usepm);
        if(v1<0 || v2<0) throw Exception("PreFlagger time range " + *str
                   + " must have positive values");
        if (usepm) {
          double pm = v2;
          v2 = v1 + pm;
          v1 -= pm;
        }
        // If time is used, values around midnight can be given.
        // Note there are 86400 seconds in a day.
        // They are split in 2 ranges.
        if (!canEndBeforeStart) {
          if (v1>=v2)
            throw Exception("PreFlagger time range " + *str + " is invalid");
        } else {
          if (v1 < 0) {
            v1 += 86400;
          }
          if (v2 > 86400) {
            v2 -= 86400;
          }
        }
        if (v1 < v2) {
          result.push_back (v1);
          result.push_back (v2);
        } else {
          result.push_back (-1);
          result.push_back (v2);
          result.push_back (v1);
          result.push_back (86401);
        }
      }
      return result;
    }

    double PreFlagger::PSet::getSeconds (const string& str, bool asTime,
                                         bool usepm)
    {
      Quantity q;
      if (asTime || usepm) {
        if (!MVAngle::read(q, str, true))
          throw Exception(
                   "PreFlagger time " + str + " is invalid");
      } else {
        // It should be a proper date/time, so MVAngle::read should fail.
        if (MVAngle::read(q, str, true))
          throw Exception(
                   "PreFlagger datetime " + str
                   + " is not a proper date/time");
        if (!MVTime::read(q, str, true))
          throw Exception(
                   "PreFlagger datetime " + str + " is invalid"
                   + " is not a proper date/time");
      }
      double v = q.getValue ("s");
      if (usepm) {
        if(v<=0) throw Exception("Preflagger time plusminus value " + str
                   + " must be positive");
      }
      return v;
    }

    vector<float> PreFlagger::PSet::fillValuePerCorr
    (const ParameterValue& value, float defVal, bool& doFlag)
    {
      // Initialize with the default value per correlation.
      vector<float> result(4);
      std::fill (result.begin(), result.end(), defVal);
      if (! value.get().empty()) {
        if (value.isVector()) {
          // Defined as a vector, take the values given.
          vector<string> valstr = value.getStringVector();
          unsigned int sz = std::min(valstr.size(), result.size());
          if (sz > 0) {
            // It contains a value, so set that flagging is done.
            doFlag = true;
            for (unsigned int i=0; i<sz; ++i) {
              if (! valstr[i].empty()) {
                result[i] = strToFloat(valstr[i]);
              }
            }
          }
        } else {
          // A single value means use it for all correlations.
          doFlag = true;
          std::fill (result.begin(), result.end(), value.getFloat());
        }
      }
      return result;
    }

    void PreFlagger::PSet::fillBLMatrix()
    {
      itsFlagOnBL = itsSelBL.hasSelection();
      if (itsFlagOnBL) {
        itsFlagBL.reference (itsSelBL.apply (*itsInfo));
      }
    }

    Vector<bool> PreFlagger::PSet::handleFreqRanges
    (const Vector<double>& chanFreqs)
    {
      unsigned int nrchan = chanFreqs.size();
      Vector<bool> selChan(nrchan, false);
      // A frequency range can be given as  value..value or value+-value.
      // Units can be given for each value; if one is given it applies to both.
      // Default unit is MHz.
      for (vector<string>::const_iterator str = itsStrFreq.begin();
           str != itsStrFreq.end(); ++str) {
        // Find the .. or +- token.
        bool usepm = false;
        string::size_type pos = str->find ("..");
        if (pos == string::npos) {
          usepm = true;
          pos = str->find ("+-");
          if (pos == string::npos)
            throw Exception("PreFlagger freqrange '" + *str
                     + "' should be range using .. or +-");
        }
        string str1 = str->substr (0, pos);
        string str2 = str->substr (pos+2);
        String u1, u2;
        double v1, v2;
        getValue (str1, v1, u1);
        // Default unit for 2nd value is that of 1st value.
        u2 = u1;
        getValue (str2, v2, u2);
        // If no unit, use MHz.
        if (u2.empty()) {
          u2 = "MHz";
        }
        // Default unit of 1st value is that of 2nd value.
        if (u1.empty()) {
          u1 = u2;
        }
        v1 = getFreqHz (v1, u1);
        v2 = getFreqHz (v2, u2);
        if (usepm) {
          double pm = v2;
          v2 = v1 + pm;
          v1 -= pm;
        }
        // Add any channel inside this range.
        for (unsigned int i=0; i<chanFreqs.size(); ++i) {
          if (chanFreqs[i] > v1  &&  chanFreqs[i] < v2) {
            selChan[i] = true;
          }
        }
      }
      return selChan;
    }

    void PreFlagger::PSet::getValue (const string& str, double& value,
                                     String& unit)
    {
      // See if a unit is given at the end.
      String v(str);
      // Remove possible trailing blanks.
      rtrim(v);
      Regex regex("[a-zA-Z]+$");
      string::size_type pos = v.index (regex);
      if (pos != String::npos) {
        unit = v.from   (pos);
        v    = v.before (pos);
      }
      // Set value and unit.
      value = strToDouble(v);
    }

    double PreFlagger::PSet::getFreqHz (double value, const String& unit)
    {
      Quantity q(value, unit);
      return q.getValue ("Hz");
    }

  } //# end namespace
}
