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

#include "ColumnReader.h"

#include <casacore/tables/Tables/ArrayColumn.h>

using casacore::ArrayColumn;

namespace DP3 {
namespace DPPP {

ColumnReader::ColumnReader(DPInput& input, const ParameterSet& parset,
                           const string& prefix, const string& column)
    : input_(input),
      name_(prefix),
      column_name_(parset.getString(prefix + "column", column)) {}

bool ColumnReader::process(const DPBuffer& buffer) {
  DPBuffer buffer_model = buffer;
  ArrayColumn<casacore::Complex> model_col(input_.table(), column_name_);
  model_col.getColumnCells(buffer.getRowNrs(), buffer_model.getData());
  getNextStep()->process(buffer_model);

  return false;
}

void ColumnReader::updateInfo(const DPInfo& _info) {
  DPStep::updateInfo(_info);
  info().setNeedVisData();
  info().setWriteData();
}

void ColumnReader::finish() { getNextStep()->finish(); }

void ColumnReader::show(std::ostream& os) const {
  os << "ColumnReader " << name_ << '\n';
  os << "  column:      " << column_name_ << '\n';
}

}  // namespace DPPP
}  // namespace DP3
