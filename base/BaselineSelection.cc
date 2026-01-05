// BaselineSelection.cc: Class to handle the baseline selection
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "BaselineSelection.h"

#include <vector>

#include <boost/algorithm/string.hpp>

#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Utilities/Regex.h>
#include <casacore/casa/version.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntenna.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/ms/MSSel/MSAntennaGram.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/Table.h>

#include <aocommon/logger.h>

#include "common/ParameterSet.h"
#include "common/ParameterValue.h"
#include "common/StreamUtil.h"

using casacore::IPosition;
using casacore::Matrix;
using casacore::MS;
using casacore::MSAntenna;
using casacore::MSAntennaParse;
using casacore::SetupNewTable;
using casacore::Table;
using dp3::common::operator<<;

using aocommon::Logger;

namespace dp3 {
namespace base {

LogAntennaParseErrors::LogAntennaParseErrors()
    : old_handler_(MSAntennaParse::thisMSAErrorHandler) {
  // The new handler logs all errors as warnings and does not throw exceptions.
  class ErrorHandler : public casacore::MSSelectionErrorHandler {
   public:
    ErrorHandler() = default;
    ~ErrorHandler() override = default;
    void reportError(const char* token,
                     const casacore::String message) override {
      Logger::Warn << message << token << '\n';
    }
  };

  // This syntax works both when ErrorHandlerPointer is a raw pointer
  // (casacore < 3.1.2) and when it's a smart pointer.
  MSAntennaParse::thisMSAErrorHandler = ErrorHandlerPointer(new ErrorHandler());
}

LogAntennaParseErrors::~LogAntennaParseErrors() {
#if CASACORE_MAJOR_VERSION < 3 ||    \
    (CASACORE_MAJOR_VERSION == 3 &&  \
     (CASACORE_MINOR_VERSION == 0 || \
      (CASACORE_MINOR_VERSION == 1 && CASACORE_PATCH_VERSION < 2)))
  // In casacore < 3.1.2 thisMSAErrorHandler is a raw pointer,
  // From casacore 3.1.2. it's a CountedPtr or another smart pointer.
  delete MSAntennaParse::thisMSAErrorHandler;
#endif

  MSAntennaParse::thisMSAErrorHandler = old_handler_;
}

BaselineSelection::BaselineSelection() {}

BaselineSelection::BaselineSelection(const common::ParameterSet& parset,
                                     const std::string& prefix, bool minmax,
                                     const std::string& defaultCorrType,
                                     const std::string& defaultBaseline)
    : itsStrBL(parset.getString(prefix + "baseline", defaultBaseline)),
      itsCorrType(parset.getString(prefix + "corrtype", defaultCorrType)),
      itsRangeBL(
          parset.getDoubleVector(prefix + "blrange", std::vector<double>())) {
  if (minmax) {
    double minbl = parset.getDouble(prefix + "blmin", -1);
    double maxbl = parset.getDouble(prefix + "blmax", -1);
    if (minbl > 0) {
      itsRangeBL.push_back(0.);
      itsRangeBL.push_back(minbl);
    }
    if (maxbl > 0) {
      itsRangeBL.push_back(maxbl);
      itsRangeBL.push_back(1e30);
    }
  }
  if (itsRangeBL.size() % 2 != 0)
    throw std::runtime_error(
        "DP3 error: uneven number of lengths in baseline range");
}

bool BaselineSelection::hasSelection() const {
  return !((itsStrBL.empty() || itsStrBL == "[]") && itsCorrType.empty() &&
           itsRangeBL.empty());
}

void BaselineSelection::show(std::ostream& os,
                             const std::string& blanks) const {
  os << "  Baseline selection:" << '\n';
  os << "    baseline:     " << blanks << itsStrBL << '\n';
  os << "    corrtype:     " << blanks << itsCorrType << '\n';
  os << "    blrange:      " << blanks << itsRangeBL << '\n';
}

Matrix<bool> BaselineSelection::apply(const DPInfo& info) const {
  // Size and initialize the selection matrix.
  const std::size_t n_antennas = info.antennaNames().size();
  Matrix<bool> selection(n_antennas, n_antennas, true);
  // Apply the various parts if given.
  if (!itsStrBL.empty() && itsStrBL != "[]") {
    handleBL(selection, info);
  }
  if (!itsCorrType.empty()) {
    handleCorrType(selection);
  }
  if (!itsRangeBL.empty()) {
    handleLength(selection, info);
  }
  return selection;
}

casacore::Vector<bool> BaselineSelection::applyVec(const DPInfo& info) const {
  Matrix<bool> sel = apply(info);
  casacore::Vector<bool> vec;
  vec.resize(info.nbaselines());
  for (unsigned int i = 0; i < info.nbaselines(); ++i) {
    vec[i] = sel(info.getAnt1()[i], info.getAnt2()[i]);
  }
  return vec;
}

void BaselineSelection::handleBL(Matrix<bool>& selectBL,
                                 const DPInfo& info) const {
  // Handle the value(s) in the baseline selection string.
  common::ParameterValue pvBL(itsStrBL);
  // The value can be a vector or an MSSelection string.
  // Alas the ParameterValue vector test cannot be used, because
  // the first character of a MSSelection string can also be [.
  // So if the first is [ and a ] is found before the end and before
  // another [, it must be a MSSelection string.
  bool mssel = true;
  if (itsStrBL[0] == '[') {
    std::string::size_type rb = itsStrBL.find(']');
    if (rb == std::string::npos)
      throw std::runtime_error("Baseline selection " + itsStrBL +
                               " has no ending ]");
    if (rb == itsStrBL.size() - 1) {
      mssel = false;
    } else {
      std::string::size_type lb = itsStrBL.find('[', 1);
      mssel = (lb == std::string::npos || lb > rb);
    }
  }
  if (!mssel) {
    // Specified as a vector of antenna name patterns.
    selectBL =
        selectBL && handleBLVector(pvBL, casacore::Vector<casacore::String>(
                                             info.antennaNames()));
  } else {
    // Specified in casacore's MSSelection format.
    selectBL = selectBL && HandleMsSelection(info);
  }
}

Matrix<bool> BaselineSelection::HandleMsSelection(const DPInfo& info) const {
  const casacore::String& antenna1_string = MS::columnName(MS::ANTENNA1);
  const casacore::String& antenna2_string = MS::columnName(MS::ANTENNA2);

  // Create a temporary MSAntenna table in memory for parsing purposes.
  SetupNewTable antenna_setup(casacore::String(),
                              MSAntenna::requiredTableDesc(), Table::New);
  Table antenna_table(antenna_setup, Table::Memory, info.antennaNames().size());
  MSAntenna ms_antenna(antenna_table);
  casacore::MSAntennaColumns ms_antenna_columns(ms_antenna);
  for (size_t i = 0; i < info.antennaNames().size(); ++i) {
    ms_antenna_columns.name().put(i, info.antennaNames()[i]);
    ms_antenna_columns.positionMeas().put(i, info.antennaPos()[i]);
  }

  // Create a temporary table holding the antenna numbers of the baselines.
  casacore::TableDesc table_description;
  table_description.addColumn(casacore::ScalarColumnDesc<int>(antenna1_string));
  table_description.addColumn(casacore::ScalarColumnDesc<int>(antenna2_string));
  SetupNewTable setup(casacore::String(), table_description, Table::New);
  Table table(setup, Table::Memory, info.nbaselines());
  casacore::ScalarColumn<int> column1(table, antenna1_string);
  casacore::ScalarColumn<int> column2(table, antenna2_string);
  for (size_t i = 0; i < info.nbaselines(); ++i) {
    column1.put(i, info.getAnt1()[i]);
    column2.put(i, info.getAnt2()[i]);
  }

  // Do the selection using the temporary tables.
  casacore::TableExprNode antenna_node1 = table.col(antenna1_string);
  casacore::TableExprNode antenna_node2 = table.col(antenna2_string);

  casacore::Vector<int> selected_antennas1;
  casacore::Vector<int> selected_antennas2;
  {
    // Overwrite the error handler to ignore errors for unknown antennas.
    dp3::base::LogAntennaParseErrors ignore_antenna_errors;

    // Parse the selection.
    // 'selected_antennas[12]' will contain the selected antenna indices.
    // 'selected_baselines' becomes an n x 2 matrix with the antenna indices
    // for each selected baseline.
    Matrix<int> selected_baselines;
    casacore::TableExprNode selection_node =
        casacore::msAntennaGramParseCommand(
            antenna_table, antenna_node1, antenna_node2, itsStrBL,
            selected_antennas1, selected_antennas2, selected_baselines);

    // msAntennaGramParseCommand may put negative indices (with unknown
    // semantics) into 'selected_antennas[12]' and 'selected_baselines'.
    // -> Apply 'selection_node' and extract the correct antenna indices.
    Table selection_table = table(selection_node);
    selected_antennas1 =
        casacore::ScalarColumn<int>(selection_table, antenna1_string)
            .getColumn();
    selected_antennas2 =
        casacore::ScalarColumn<int>(selection_table, antenna2_string)
            .getColumn();
  }

  // Convert selected_antennas[12] to a selection matrix.
  Matrix<bool> selection(info.nantenna(), info.nantenna(), false);
  for (size_t bl = 0; bl < selected_antennas1.size(); ++bl) {
    const int a1 = selected_antennas1[bl];
    const int a2 = selected_antennas2[bl];
    selection(a1, a2) = true;
    selection(a2, a1) = true;
  }
  return selection;
}

Matrix<bool> BaselineSelection::handleBLVector(
    const common::ParameterValue& pvBL,
    const casacore::Vector<casacore::String>& antNames) const {
  Matrix<bool> sel(antNames.size(), antNames.size());
  sel = false;
  std::vector<common::ParameterValue> pairs = pvBL.getVector();
  // Each ParameterValue can be a single value (antenna) or a pair of
  // values (a baseline).
  // Note that [ant1,ant2] is somewhat ambiguous; it means two antennae,
  // but one might think it means a baseline [[ant1,ant2]].
  if (pairs.size() == 2 && !(pairs[0].isVector() || pairs[1].isVector())) {
    Logger::Warn << "PreFlagger baseline " << pvBL.get()
                 << " means two antennae, but is somewhat ambigious; "
                    "it's more clear to use [[ant1],[ant2]]\n";
  }
  for (unsigned int i = 0; i < pairs.size(); ++i) {
    std::vector<std::string> bl = pairs[i].getStringVector();
    if (bl.size() == 1) {
      // Turn the given antenna name pattern into a regex.
      casacore::Regex regex(casacore::Regex::fromPattern(bl[0]));
      int nmatch = 0;
      // Loop through all antenna names and set matrix for matching ones.
      for (unsigned int i2 = 0; i2 < antNames.size(); ++i2) {
        if (casacore::String(antNames[i2]).matches(regex)) {
          nmatch++;
          // Antenna matches, so set all corresponding flags.
          for (unsigned int j = 0; j < antNames.size(); ++j) {
            sel(i2, j) = true;
            sel(j, i2) = true;
          }
        }
      }
      if (nmatch == 0) {
        Logger::Warn << "PreFlagger: no matches for antenna name pattern ["
                     << bl[0] << "]\n";
      }
    } else {
      if (bl.size() != 2)
        throw std::runtime_error(
            "PreFlagger baseline "
            " should contain 1 or 2 antenna name patterns");
      // Turn the given antenna name pattern into a regex.
      casacore::Regex regex1(casacore::Regex::fromPattern(bl[0]));
      casacore::Regex regex2(casacore::Regex::fromPattern(bl[1]));
      int nmatch = 0;
      // Loop through all antenna names and set matrix for matching ones.
      for (unsigned int i2 = 0; i2 < antNames.size(); ++i2) {
        if (casacore::String(antNames[i2]).matches(regex2)) {
          // Antenna2 matches, now try Antenna1.
          for (unsigned int i1 = 0; i1 < antNames.size(); ++i1) {
            if (casacore::String(antNames[i1]).matches(regex1)) {
              nmatch++;
              sel(i1, i2) = true;
              sel(i2, i1) = true;
            }
          }
        }
      }
      if (nmatch == 0) {
        Logger::Warn << "PreFlagger: no matches for baseline name pattern ["
                     << bl[0] << ',' << bl[1] << "]\n";
      }
    }
  }
  return sel;
}

void BaselineSelection::handleCorrType(Matrix<bool>& selectBL) const {
  // Process corrtype if given.
  string corrType = boost::to_lower_copy(itsCorrType);
  if (corrType != "auto" && corrType != "cross")
    throw std::runtime_error(
        "DP3 corrType " + corrType +
        " is invalid; must be auto, cross, or empty string");
  if (corrType == "auto") {
    casacore::Vector<bool> diag = selectBL.diagonal().copy();
    selectBL = false;
    selectBL.diagonal() = diag;
  } else {
    selectBL.diagonal() = false;
  }
}

void BaselineSelection::handleLength(Matrix<bool>& selectBL,
                                     const DPInfo& info) const {
  // Get baseline lengths.
  const std::vector<double>& blength = info.getBaselineLengths();
  const std::vector<int>& ant1 = info.getAnt1();
  const std::vector<int>& ant2 = info.getAnt2();
  for (unsigned int i = 0; i < ant1.size(); ++i) {
    // Clear selection if no range matches.
    bool match = false;
    for (unsigned int j = 0; j < itsRangeBL.size(); j += 2) {
      if (blength[i] >= itsRangeBL[j] && blength[i] <= itsRangeBL[j + 1]) {
        match = true;
        break;
      }
    }
    if (!match) {
      int a1 = ant1[i];
      int a2 = ant2[i];
      selectBL(a1, a2) = false;
      selectBL(a2, a1) = false;
    }
  }
}

}  // namespace base
}  // namespace dp3
