// OneApplyCal.cc: DPPP step class to apply a calibration correction to the data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "OneApplyCal.h"

#include "ApplyCal.h"
#include "MSReader.h"

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "../common/ParameterSet.h"
#include "../common/StringTools.h"

#include <aocommon/dynamicfor.h>

#include <casacore/casa/Arrays/ArrayMath.h>

#include <iostream>
#include <limits>
#include <algorithm>
#include <iomanip>

#include <boost/algorithm/string/case_conv.hpp>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using std::vector;

namespace dp3 {
namespace steps {
// Initialize private static
std::mutex OneApplyCal::theirHDF5Mutex;

OneApplyCal::OneApplyCal(const common::ParameterSet& parset,
                         const std::string& prefix,
                         const std::string& defaultPrefix, bool substep,
                         std::string predictDirection)
    : itsName(prefix),
      itsParmDBName(parset.isDefined(prefix + "parmdb")
                        ? parset.getString(prefix + "parmdb")
                        : parset.getString(defaultPrefix + "parmdb", "")),
      itsParmDBOnDisk(!itsParmDBName.empty()),
      itsUseH5Parm(itsParmDBName.find(".h5") != std::string::npos),
      itsSolSetName(parset.isDefined(prefix + "solset")
                        ? parset.getString(prefix + "solset")
                        : parset.getString(defaultPrefix + "solset", "")),
      itsSigmaMMSE(parset.isDefined(prefix + "MMSE.Sigma")
                       ? parset.getDouble(prefix + "MMSE.Sigma")
                       : parset.getDouble(defaultPrefix + "MMSE.Sigma", 0.)),
      itsUpdateWeights(
          parset.isDefined(prefix + "updateweights")
              ? parset.getBool(prefix + "updateweights")
              : parset.getBool(defaultPrefix + "updateweights", false)),
      itsCount(0),
      itsJonesParameters(nullptr),
      itsTimeStep(0),
      itsNCorr(0),
      itsTimeInterval(-1),
      itsLastTime(-1),
      itsUseAP(false) {
  if (substep) {
    itsInvert = false;
  } else {
    itsInvert = (parset.isDefined(prefix + "invert")
                     ? parset.getBool(prefix + "invert")
                     : parset.getBool(defaultPrefix + "invert", true));
  }

  if (itsUseH5Parm) {
    const std::string interpolationStr =
        (parset.isDefined(prefix + "interpolation")
             ? parset.getString(prefix + "interpolation")
             : parset.getString(prefix + "interpolation", "nearest"));
    if (interpolationStr == "nearest") {
      itsInterpolationType = JonesParameters::InterpolationType::NEAREST;
    } else if (interpolationStr == "linear") {
      itsInterpolationType = JonesParameters::InterpolationType::LINEAR;
    } else {
      throw std::runtime_error("Unsupported interpolation mode: " +
                               interpolationStr);
    }
    itsTimeSlotsPerParmUpdate = 0;
    const std::string directionStr =
        (parset.isDefined(prefix + "direction")
             ? parset.getString(prefix + "direction")
             : parset.getString(defaultPrefix + "direction", predictDirection));

    if (itsParmDBOnDisk) {
      itsH5Parm = schaapcommon::h5parm::H5Parm(itsParmDBName, false, false,
                                               itsSolSetName);
      // The correction only makes sense if the parmdb is being read from disk.
      itsSolTabName = (parset.isDefined(prefix + "correction")
                           ? parset.getString(prefix + "correction")
                           : parset.getString(defaultPrefix + "correction"));
    }
    const std::string missingAntennaBehaviorStr =
        (parset.isDefined(prefix + "missingantennabehavior")
             ? parset.getString(prefix + "missingantennabehavior")
             : parset.getString(defaultPrefix + "missingantennabehavior",
                                "error"));
    itsMissingAntennaBehavior = JonesParameters::StringToMissingAntennaBehavior(
        missingAntennaBehaviorStr);

    if (itsParmDBOnDisk) {
      std::vector<string> solTabs = parset.getStringVector(
          prefix + "soltab", std::vector<string>{"amplitude000", "phase000"});
      setCorrectType(solTabs);

      itsDirection = 0;
      if (directionStr.empty()) {
        if (itsSolTab.HasAxis("dir") && itsSolTab.GetAxis("dir").size != 1)
          throw std::runtime_error(
              "If the soltab contains multiple directions, the direction to "
              "be "
              "applied in applycal should be specified");
        // If there is only one direction, silently assume it is the right one
      } else if (itsSolTab.HasAxis("dir") &&
                 itsSolTab.GetAxis("dir").size > 1) {
        itsDirection = itsSolTab.GetDirIndex(directionStr);
      }
    }
  } else {
    itsMissingAntennaBehavior = JonesParameters::MissingAntennaBehavior::kError;
    itsTimeSlotsPerParmUpdate =
        parset.isDefined(prefix + "timeslotsperparmupdate")
            ? parset.getInt(prefix + "timeslotsperparmupdate")
            : parset.getInt(defaultPrefix + "timeslotsperparmupdate", 500);
    const std::string correctTypeStr = boost::to_lower_copy(
        parset.isDefined(prefix + "correction")
            ? parset.getString(prefix + "correction")
            : parset.getString(defaultPrefix + "correction", "gain"));
    // TODO this shouldn't be dependent on the H5 parm type, but use
    // more logical names (e.g. "diagonal" instead of "gain").
    itsCorrectType =
        JonesParameters::H5ParmTypeStringToGainType(correctTypeStr);
  }

  if (itsCorrectType == GainType::kFullJones && itsUpdateWeights) {
    if (!itsInvert)
      throw std::runtime_error(
          "Updating weights has not been implemented for invert=false and "
          "fulljones");
  }
}

void OneApplyCal::setCorrectType(std::vector<std::string>& solTabs) {
  if (itsSolTabName == "fulljones") {
    if (solTabs.size() != 2)
      throw std::runtime_error(
          "The soltab parameter requires two soltabs for fulljones "
          "correction (amplitude and phase)");
    itsSolTab = itsH5Parm.GetSolTab(solTabs[0]);
    itsSolTab2 = itsH5Parm.GetSolTab(solTabs[1]);
    itsSolTabName =
        solTabs[0] + ", " +
        solTabs[1];  // this is only so that show() shows these tables
    itsCorrectType = GainType::kFullJones;
  } else {
    itsSolTab = itsH5Parm.GetSolTab(itsSolTabName);
    itsCorrectType =
        JonesParameters::H5ParmTypeStringToGainType(itsSolTab.GetType());
  }
  if (itsCorrectType == GainType::kDiagonalPhase && nPol("") == 1) {
    itsCorrectType = GainType::kScalarPhase;
  }
  if (itsCorrectType == GainType::kDiagonalAmplitude && nPol("") == 1) {
    itsCorrectType = GainType::kScalarAmplitude;
  }
}

OneApplyCal::~OneApplyCal() {}

void OneApplyCal::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  itsTimeInterval = infoIn.timeInterval();
  itsNCorr = infoIn.ncorr();

  if (itsNCorr != 4)
    throw std::runtime_error("Applycal only works with 4 correlations");

  if (itsParmDBOnDisk) {
    if (itsUseH5Parm) {
      itsTimeSlotsPerParmUpdate = info().ntime();
    } else {  // Use ParmDB
      itsParmDB = std::make_shared<parmdb::ParmFacade>(itsParmDBName);
    }

    // Detect if full jones solutions are present
    if (!itsUseH5Parm &&
        (itsCorrectType == GainType::kDiagonalComplex ||
         itsCorrectType == GainType::kFullJones) &&
        (itsParmDB->getNames("Gain:0:1:*").size() +
             itsParmDB->getDefNames("Gain:0:1:*").size() >
         0)) {
      itsCorrectType = GainType::kFullJones;
    }

    // Detect if solutions are saved as Real/Imag or Ampl/Phase
    if (itsCorrectType == GainType::kDiagonalComplex ||
        itsCorrectType == GainType::kFullJones) {
      if (itsUseH5Parm) {
        // H5Parm uses amplitude / phase by definition
        itsUseAP = true;
      } else {
        // Determine from values present in parmdb what to use
        if (!itsParmDB->getNames("Gain:0:0:Real*").empty()) {
          // Values with :Real present
          itsUseAP = false;
        } else if (!itsParmDB->getNames("Gain:0:0:Ampl*").empty() ||
                   !itsParmDB->getNames("Phase:0:0:Ampl*").empty()) {
          // Values with :Ampl present
          itsUseAP = true;
        } else if (!itsParmDB->getDefNames("Gain:0:0:Real*").empty()) {
          // Defvalues with :Real present
          itsUseAP = false;
        } else if (!itsParmDB->getDefNames("Gain:0:0:Ampl*").empty() ||
                   !itsParmDB->getDefNames("Gain:0:0:Phase*").empty()) {
          // Defvalues with :Ampl present
          itsUseAP = true;
        } else {
          throw std::runtime_error("No gains found in parmdb " + itsParmDBName);
        }
      }
    }

    if (itsCorrectType == GainType::kDiagonalComplex) {
      if (itsUseAP) {
        itsParmExprs.push_back("Gain:0:0:Ampl");
        itsParmExprs.push_back("Gain:0:0:Phase");
        itsParmExprs.push_back("Gain:1:1:Ampl");
        itsParmExprs.push_back("Gain:1:1:Phase");
      } else {
        itsParmExprs.push_back("Gain:0:0:Real");
        itsParmExprs.push_back("Gain:0:0:Imag");
        itsParmExprs.push_back("Gain:1:1:Real");
        itsParmExprs.push_back("Gain:1:1:Imag");
      }
    } else if (itsCorrectType == GainType::kFullJones) {
      if (itsUseAP) {
        itsParmExprs.push_back("Gain:0:0:Ampl");
        itsParmExprs.push_back("Gain:0:0:Phase");
        itsParmExprs.push_back("Gain:0:1:Ampl");
        itsParmExprs.push_back("Gain:0:1:Phase");
        itsParmExprs.push_back("Gain:1:0:Ampl");
        itsParmExprs.push_back("Gain:1:0:Phase");
        itsParmExprs.push_back("Gain:1:1:Ampl");
        itsParmExprs.push_back("Gain:1:1:Phase");
      } else {
        itsParmExprs.push_back("Gain:0:0:Real");
        itsParmExprs.push_back("Gain:0:0:Imag");
        itsParmExprs.push_back("Gain:0:1:Real");
        itsParmExprs.push_back("Gain:0:1:Imag");
        itsParmExprs.push_back("Gain:1:0:Real");
        itsParmExprs.push_back("Gain:1:0:Imag");
        itsParmExprs.push_back("Gain:1:1:Real");
        itsParmExprs.push_back("Gain:1:1:Imag");
      }
    } else if (itsCorrectType == GainType::kTec) {
      if (nPol("TEC") == 1) {
        itsParmExprs.push_back("TEC");
      } else {
        itsParmExprs.push_back("TEC:0");
        itsParmExprs.push_back("TEC:1");
      }
    } else if (itsCorrectType == GainType::kClock) {
      if (nPol("Clock") == 1) {
        itsParmExprs.push_back("Clock");
      } else {
        itsParmExprs.push_back("Clock:0");
        itsParmExprs.push_back("Clock:1");
      }
    } else if (itsCorrectType == GainType::kRotationAngle) {
      itsParmExprs.push_back("{Common,}RotationAngle");
    } else if (itsCorrectType == GainType::kScalarPhase) {
      itsParmExprs.push_back("{Common,}ScalarPhase");
    } else if (itsCorrectType == GainType::kRotationMeasure) {
      itsParmExprs.push_back("RotationMeasure");
    } else if (itsCorrectType == GainType::kScalarAmplitude) {
      itsParmExprs.push_back("{Common,}ScalarAmplitude");
    } else if (itsCorrectType == GainType::kDiagonalPhase) {
      if (!itsUseH5Parm)
        throw std::runtime_error("A H5Parm is required for phase correction");
      itsParmExprs.push_back("Phase:0");
      itsParmExprs.push_back("Phase:1");
    } else if (itsCorrectType == GainType::kDiagonalAmplitude) {
      if (!itsUseH5Parm)
        throw std::runtime_error(
            "A H5Parm is required for amplitude correction");
      itsParmExprs.push_back("Amplitude:0");
      itsParmExprs.push_back("Amplitude:1");
    } else {
      throw std::runtime_error(
          "Correction type " +
          JonesParameters::GainTypeToHumanReadableString(itsCorrectType) +
          " is unknown");
    }
  }

  itsFlagCounter.init(getInfo());

  // Check that channels are evenly spaced
  if (!itsUseH5Parm && !info().channelsAreRegular()) {
    throw std::runtime_error(
        "ApplyCal with parmdb requires evenly spaced channels.");
  }
}

void OneApplyCal::show(std::ostream& os) const {
  os << "ApplyCal " << itsName << '\n';
  if (itsUseH5Parm) {
    os << "  H5Parm:         " << itsParmDBName << '\n';
    os << "    SolSet:       " << itsH5Parm.GetSolSetName() << '\n';
    os << "    SolTab:       " << itsSolTabName << '\n';
    os << "  Direction:      " << itsDirection << '\n';
    os << "  Interpolation:  "
       << (itsInterpolationType == JonesParameters::InterpolationType::NEAREST
               ? "nearest"
               : "linear")
       << '\n';
    os << "  Missing antennas: "
       << JonesParameters::MissingAntennaBehaviorToString(
              itsMissingAntennaBehavior)
       << '\n';
  } else if (itsParmDBOnDisk) {
    os << "  Parmdb:         " << itsParmDBName << '\n';
  } else {
    os << "  Parm solutions read from buffer"
       << "\n";
  }
  os << "  Correction:       "
     << JonesParameters::GainTypeToHumanReadableString(itsCorrectType) << '\n';
  if (itsCorrectType == GainType::kDiagonalComplex ||
      itsCorrectType == GainType::kFullJones) {
    os << "    Ampl/Phase:   " << std::boolalpha << itsUseAP << '\n';
  }
  os << "  Update weights:   " << std::boolalpha << itsUpdateWeights << '\n';
  os << "  Invert:           " << std::boolalpha << itsInvert << '\n';
  if (itsInvert) {
    os << "    SigmaMMSE:    " << itsSigmaMMSE << '\n';
  }
  os << "  TimeSlotsPerParmUpdate: " << itsTimeSlotsPerParmUpdate << '\n';
}

void OneApplyCal::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " OneApplyCal " << itsName << '\n';
}

bool OneApplyCal::process(std::unique_ptr<DPBuffer> buffer) {
  itsTimer.start();

  if (buffer->GetTime() > itsLastTime) {
    if (itsParmDBOnDisk && itsUseH5Parm) {
      updateParmsH5(buffer->GetTime());
    } else if (itsParmDBOnDisk) {
      updateParmsParmDB(buffer->GetTime());
    } else {
      if (buffer->GetSolution().size() == 0) {
        throw std::runtime_error(
            "No buffer stored before OneApplyCal step. Ensure a solution is "
            "computed before this step, or specify a parmdb/h5parm file in "
            "the applycal step.");
      }
      /*
       * Read solutions from buffer instead.
       */
      // Make parameters complex
      GainType gain_type = itsCorrectType;
      if (itsCorrectType == GainType::kDiagonalComplex && !itsUseAP) {
        gain_type = GainType::kDiagonalRealImaginary;
      } else if (itsCorrectType == GainType::kFullJones && !itsUseAP) {
        gain_type = GainType::kFullJonesRealImaginary;
      }

      vector<double> times(info().ntime());
      for (size_t t = 0; t < times.size(); ++t) {
        // time centroids
        times[t] = info().startTime() + (t + 0.5) * info().timeInterval();
      }

      // Validate that the data is in the correct shape
      const size_t n_chan = buffer->GetData().shape(1);
      const size_t n_corrs =
          buffer->GetSolution()[0].size() / info().antennaNames().size();
      if (buffer->GetSolution().size() != n_chan ||
          (n_corrs != 2 && n_corrs != 4)) {
        throw std::runtime_error(
            "The solution is not in the correct shape. Was the solution "
            "computed on different data than it is being applied on?");
      }

      itsJonesParameters = std::make_unique<JonesParameters>(
          info().chanFreqs(), times, info().antennaNames(), gain_type,
          buffer->GetSolution(), itsInvert, itsSigmaMMSE);
    }
    itsTimeStep = 0;
  } else {
    itsTimeStep++;
  }

  // Loop through all baselines in the buffer.
  const size_t n_bl = buffer->GetData().shape(0);
  const size_t n_chan = buffer->GetData().shape(1);
  const casacore::Cube<casacore::Complex>& gains =
      itsJonesParameters->GetParms();
  const size_t n_corr = gains.shape()[0];

  aocommon::DynamicFor<size_t> loop;
  loop.Run(0, n_bl, [&](size_t bl) {
    const unsigned int ant_a = info().getAnt1()[bl];
    const unsigned int ant_b = info().getAnt2()[bl];

    for (size_t chan = 0; chan < n_chan; chan++) {
      const unsigned int time_freq_offset =
          (itsTimeStep * info().nchan()) + chan;
      const std::complex<float>* gain_a = &gains(0, ant_a, time_freq_offset);
      const std::complex<float>* gain_b = &gains(0, ant_b, time_freq_offset);
      if (n_corr > 2) {
        ApplyCal::ApplyFull(aocommon::MatrixComplexFloat2x2{gain_a},
                            aocommon::MatrixComplexFloat2x2{gain_b}, *buffer,
                            bl, chan, itsUpdateWeights, itsFlagCounter);
      } else {
        ApplyCal::ApplyDiag(aocommon::DiagonalMatrixComplexFloat2x2{gain_a},
                            aocommon::DiagonalMatrixComplexFloat2x2{gain_b},
                            *buffer, bl, chan, itsUpdateWeights,
                            itsFlagCounter);
      }
    }
  });

  itsTimer.stop();
  getNextStep()->process(std::move(buffer));

  itsCount++;
  return true;
}

void OneApplyCal::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

void OneApplyCal::updateParmsH5(const double bufStartTime) {
  itsLastTime = bufStartTime - 0.5 * itsTimeInterval +
                itsTimeSlotsPerParmUpdate * itsTimeInterval;

  const double lastMSTime =
      info().startTime() + info().ntime() * itsTimeInterval;
  if (itsLastTime > lastMSTime &&
      !casacore::nearAbs(itsLastTime, lastMSTime, 1.e-3)) {
    itsLastTime = lastMSTime;
  }

  std::lock_guard<std::mutex> lock(theirHDF5Mutex);

  // Figure out whether time or frequency is first axis
  if (itsSolTab.HasAxis("freq") && itsSolTab.HasAxis("time") &&
      itsSolTab.GetAxisIndex("freq") < itsSolTab.GetAxisIndex("time")) {
    throw std::runtime_error("Fastest varying axis should be freq");
  }

  vector<double> times(info().ntime());
  for (size_t t = 0; t < times.size(); ++t) {
    // TODO(AST-1078) for investigating the 0.5 offset.
    // time centroids
    times[t] = info().startTime() + (t + 0.5) * info().timeInterval();
  }

  std::vector<std::string> ant_names;
  for (const std::string& name : info().antennaNames()) {
    ant_names.push_back(name);
  }

  itsJonesParameters = std::make_unique<JonesParameters>(
      info().chanFreqs(), times, ant_names, itsCorrectType,
      itsInterpolationType, itsDirection, &itsSolTab, &itsSolTab2, itsInvert,
      itsSigmaMMSE, itsParmExprs.size(), itsMissingAntennaBehavior);
}

void OneApplyCal::updateParmsParmDB(const double bufStartTime) {
  unsigned int numAnts = info().antennaNames().size();

  vector<vector<vector<double>>> parmvalues;
  parmvalues.resize(itsParmExprs.size());
  for (size_t i = 0; i < parmvalues.size(); ++i) {
    parmvalues[i].resize(numAnts);
  }

  unsigned int numFreqs(info().chanFreqs().size());
  double freqInterval(info().chanWidths()[0]);
  if (numFreqs > 1) {  // Handle data with evenly spaced gaps between channels
    freqInterval = info().chanFreqs()[1] - info().chanFreqs()[0];
  }
  double minFreq(info().chanFreqs()[0] - 0.5 * freqInterval);
  double maxFreq(info().chanFreqs()[numFreqs - 1] + 0.5 * freqInterval);
  itsLastTime = bufStartTime - 0.5 * itsTimeInterval +
                itsTimeSlotsPerParmUpdate * itsTimeInterval;
  unsigned int numTimes = itsTimeSlotsPerParmUpdate;

  double lastMSTime = info().startTime() + info().ntime() * itsTimeInterval;
  if (itsLastTime > lastMSTime &&
      !casacore::nearAbs(itsLastTime, lastMSTime, 1.e-3)) {
    itsLastTime = lastMSTime;
    numTimes = info().ntime() % itsTimeSlotsPerParmUpdate;
  }

  std::vector<double> times(numTimes);
  for (size_t t = 0; t < times.size(); ++t) {
    // time centroids
    times[t] = bufStartTime + (t + 0.5) * itsTimeInterval;
  }

  std::map<std::string, std::vector<double>> parmMap;
  std::map<std::string, std::vector<double>>::iterator parmIt;

  unsigned int tfDomainSize = numTimes * numFreqs;

  for (unsigned int parmExprNum = 0; parmExprNum < itsParmExprs.size();
       ++parmExprNum) {
    // parmMap contains parameter values for all antennas
    parmMap = itsParmDB->getValuesMap(
        itsParmExprs[parmExprNum] + "{:phase_center,}*", minFreq, maxFreq,
        freqInterval, bufStartTime - 0.5 * itsTimeInterval, itsLastTime,
        itsTimeInterval, true);

    std::string parmExpr = itsParmExprs[parmExprNum];

    // Resolve {Common,}Bla to CommonBla or Bla
    if (!parmMap.empty() && parmExpr.find("{Common,}") != std::string::npos) {
      // Take the name of the first parm, e.g. Bla:CS001, and remove the
      // antenna name
      unsigned int colonPos = (parmMap.begin()->first).find(":");
      parmExpr = (parmMap.begin()->first).substr(0, colonPos);
    }

    std::string name_postfix = "";
    // Remove :phase_center postfix
    if (!parmMap.empty()) {
      // Take the name of the first parm, e.g. Bla:CS001, and remove the
      // antenna name If necessary, append :phase_center
      if ((parmMap.begin()->first).find(":phase_center") != std::string::npos) {
        name_postfix = ":phase_center";
      }
    }

    for (unsigned int ant = 0; ant < numAnts; ++ant) {
      parmIt =
          parmMap.find(parmExpr + ":" +
                       std::string(info().antennaNames()[ant]) + name_postfix);

      if (parmIt != parmMap.end()) {
        parmvalues[parmExprNum][ant].swap(parmIt->second);
        if (parmvalues[parmExprNum][ant].size() != tfDomainSize)
          throw std::runtime_error("Size of parmvalue != tfDomainSize");
      } else {  // No value found, try default
        casacore::Array<double> defValues;
        double defValue;

        std::string defParmNameAntenna =
            parmExpr + ":" + std::string(info().antennaNames()[ant]) +
            name_postfix;
        if (itsParmDB->getDefValues(defParmNameAntenna).size() ==
            1) {  // Default for antenna
          itsParmDB->getDefValues(defParmNameAntenna).get(0, defValues);
          if (defValues.size() != 1)
            throw std::runtime_error(
                "Multiple default values found in parmdb for " +
                defParmNameAntenna + ". " +
                "Did you unintentially overwrite an existing parmdb?");
          defValue = defValues.data()[0];
        } else if (itsParmDB->getDefValues(parmExpr).size() ==
                   1) {  // Default value
          itsParmDB->getDefValues(parmExpr).get(0, defValues);
          if (defValues.size() != 1)
            throw std::runtime_error("Size of defValues != 1");
          defValue = defValues.data()[0];
        } else if (parmExpr.substr(0, 5) == "Gain:") {
          defValue = 0.;
        } else {
          throw std::runtime_error(
              "No parameter value found for " + parmExpr + ":" +
              std::string(info().antennaNames()[ant]) + name_postfix);
        }

        parmvalues[parmExprNum][ant].resize(tfDomainSize);
        for (unsigned int tf = 0; tf < tfDomainSize; ++tf) {
          parmvalues[parmExprNum][ant][tf] = defValue;
        }
      }
    }
  }

  if (parmvalues[0][0].size() > tfDomainSize)  // Catches multiple matches
    throw std::runtime_error("Parameter was found multiple times in ParmDB");

  // Make parameters complex
  GainType gain_type = itsCorrectType;
  if (itsCorrectType == GainType::kDiagonalComplex && !itsUseAP) {
    gain_type = GainType::kDiagonalRealImaginary;
  } else if (itsCorrectType == GainType::kFullJones && !itsUseAP) {
    gain_type = GainType::kFullJonesRealImaginary;
  }

  std::vector<std::string> ant_names;
  for (const std::string& name : info().antennaNames()) {
    ant_names.push_back(name);
  }

  itsJonesParameters = std::make_unique<JonesParameters>(
      info().chanFreqs(), times, ant_names, gain_type, itsInterpolationType,
      itsDirection, std::move(parmvalues), itsInvert, itsSigmaMMSE);
}

unsigned int OneApplyCal::nPol(const std::string& parmName) {
  if (itsUseH5Parm) {
    if (!itsSolTab.HasAxis("pol")) {
      return 1;
    } else {
      return itsSolTab.GetAxis("pol").size;
    }
  } else {  // Use ParmDB
    if (itsParmDB->getNames(parmName + ":0:*").empty() &&
        itsParmDB->getDefNames(parmName + ":0:*").empty()) {
      return 1;
    } else {
      return 2;
    }
  }
}

void OneApplyCal::showCounts(std::ostream& os) const {
  os << "\nFlags set by OneApplyCal " << itsName;
  os << "\n=======================\n";
  itsFlagCounter.showBaseline(os, itsCount);
  itsFlagCounter.showChannel(os, itsCount);
}

}  // namespace steps
}  // namespace dp3
