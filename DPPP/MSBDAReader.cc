// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#include "MSBDAReader.h"
#include "BDABuffer.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Exceptions.h"

#include <EveryBeam/load.h>
#include <EveryBeam/lofarreadutils.h>

#include "../Common/ParameterSet.h"
#include "../Common/BaselineSelect.h"

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSAntennaParse.h>
#include <casacore/ms/MSSel/MSSelectionErrorHandler.h>
#include <casacore/casa/version.h>

#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/OS/Conversion.h>

#include <iostream>

using namespace casacore;

namespace DP3 {
namespace DPPP {

MSBDAReader::MSBDAReader()
    : readVisData_(False), lastMSTime_(0), nread_(0), ninserted_(0) {}

MSBDAReader::MSBDAReader(const string& msName, const ParameterSet& parset,
                         const string& prefix)
    : readVisData_(False), lastMSTime_(0), nread_(0), ninserted_(0) {}

MSBDAReader::~MSBDAReader() {}

void MSBDAReader::updateInfo(const DPInfo& dpInfo) {}

std::string MSBDAReader::msName() const { return ""; }

void MSBDAReader::setReadVisData(bool readVisData) {
  readVisData_ = readVisData;
}

bool MSBDAReader::process(std::unique_ptr<BDABuffer>) { return true; }

void MSBDAReader::finish() { getNextStep()->finish(); }

void MSBDAReader::show(std::ostream& os) const {}

void MSBDAReader::showCounts(std::ostream& os) const {}

void MSBDAReader::showTimings(std::ostream& os, double duration) const {}

void MSBDAReader::getUVW(const RefRows& rowNrs, double time, DPBuffer& buf) {}

void MSBDAReader::getWeights(const RefRows& rowNrs, DPBuffer& buf) {}

bool MSBDAReader::getFullResFlags(const RefRows& rowNrs, DPBuffer& buf) {
  return true;
}

void MSBDAReader::getModelData(const casacore::RefRows& rowNrs,
                               casacore::Cube<casacore::Complex>& arr) {}

void MSBDAReader::fillBeamInfo(
    vector<everybeam::Station::Ptr>& vec,
    const casacore::Vector<casacore::String>& antNames) {}

}  // namespace DPPP
}  // namespace DP3
