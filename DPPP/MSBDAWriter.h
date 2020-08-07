// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
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

/// @file
/// @brief Step for writing BDA data to an MS.
/// @author Maik Nijhuis

#ifndef DPPP_MSBDAWRITER_H
#define DPPP_MSBDAWRITER_H

#include <casacore/tables/Tables/Table.h>

#include "DPStep.h"
#include "MSReader.h"
#include "../Common/ParameterSet.h"

using casacore::Table;

namespace DP3 {
namespace DPPP {

class MSBDAWriter : public DPStep {
 public:
  MSBDAWriter(MSReader*, const string&, const ParameterSet&,
              const string&);
  virtual ~MSBDAWriter();

  /// Update the general info.
  void updateInfo(const DPInfo&) override;

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Create the MS by cloning all subtables from the input MS.
  /// All output columns in the main table are using normal storage managers.
  void createMS(const DPInfo&);

 private:
  void createMainTable(const DPInfo&);
  void createBDATimeAxis(const DPInfo&);

 private:
  MSReader* reader_;
  const std::string outName_;
  const ParameterSet parset_;
  const std::string prefix_;

  bool overwrite_;
  casacore::Table ms_;
};

}  // namespace DPPP
}  // namespace DP3

#endif
