// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
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
      column_name_(parset.getString(prefix + "column", column)),
      operation_(Operation::kReplace),
      buffer_() {
  const std::string operation =
      parset.getString(prefix + "operation", "replace");
  if (operation == "replace") {
    operation_ = Operation::kReplace;
  } else if (operation == "add") {
    operation_ = Operation::kAdd;
  } else if (operation == "subtract") {
    operation_ = Operation::kSubtract;
  } else {
    throw std::invalid_argument("Invalid MsColumnReader operation " +
                                operation);
  }
}

bool MsColumnReader::process(const DPBuffer& buffer) {
  buffer_.copy(buffer);
  ArrayColumn<casacore::Complex> model_col(table_, column_name_);
  model_col.getColumnCells(buffer.getRowNrs(), buffer_.getData());

  if (operation_ == Operation::kAdd) {
    buffer_.setData(buffer.getData() + buffer_.getData());
  } else if (operation_ == Operation::kSubtract) {
    buffer_.setData(buffer.getData() - buffer_.getData());
  }

  getNextStep()->process(buffer_);

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
  os << "  operation:   ";
  switch (operation_) {
    case Operation::kReplace:
      os << "replace";
      break;
    case Operation::kAdd:
      os << "add";
      break;
    case Operation::kSubtract:
      os << "subtract";
      break;
  }
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
