// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MsColumnReader.h"

#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/tables/Tables/ArrayColumn.h>

using casacore::ArrayColumn;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

MsColumnReader::MsColumnReader(const common::ParameterSet& parset,
                               const std::string& prefix,
                               const std::string& column)
    : table_(),
      name_(prefix),
      column_name_(parset.getString(prefix + "column", column)) {}

bool MsColumnReader::process(std::unique_ptr<DPBuffer> buffer) {
  buffer->ResizeData(
      {getInfo().nbaselines(), getInfo().nchan(), getInfo().ncorr()});
  const casacore::IPosition shape(3, getInfo().ncorr(), getInfo().nchan(),
                                  getInfo().nbaselines());
  casacore::Cube<casacore::Complex> data(shape, buffer->GetData().data(),
                                         casacore::SHARE);

  ArrayColumn<casacore::Complex> model_column(table_, column_name_);
  model_column.getColumnCells(buffer->getRowNrs(), data);

  getNextStep()->process(std::move(buffer));

  return false;
}

void MsColumnReader::updateInfo(const DPInfo& _info) {
  Step::updateInfo(_info);
  table_ = casacore::Table(_info.msName());
}

void MsColumnReader::finish() { getNextStep()->finish(); }

void MsColumnReader::show(std::ostream& os) const {
  os << "MsColumnReader " << name_ << '\n';
  os << "  column:      " << column_name_ << '\n';
  os << '\n';
}

void MsColumnReader::showTimings(std::ostream& os,
                                 [[maybe_unused]] double duration) const {
  os << " MsColumnReader " << name_ << '\n';
}

base::Direction MsColumnReader::GetFirstDirection() const {
  using casacore::MDirection;
  const MDirection dirJ2000(
      MDirection::Convert(getInfo().phaseCenter(), MDirection::J2000)());
  const casacore::Quantum<casacore::Vector<double>> angles =
      dirJ2000.getAngle();
  return {angles.getBaseValue()[0], angles.getBaseValue()[1]};
}

}  // namespace steps
}  // namespace dp3
