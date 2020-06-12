#ifndef DPPP_SETBEAM_H
#define DPPP_SETBEAM_H

/// @file
/// @brief DPPP step class to set the beam keywords in a ms

#include "DPInput.h"
#include "DPBuffer.h"
#include "Position.h"

#include <casacore/measures/Measures/MDirection.h>

namespace DP3 {

class ParameterSet;

namespace DPPP {
  
class SetBeam final : public DPStep
{
public:
  /// Parameters are obtained from the parset using the given prefix.
  SetBeam (DPInput* input, const ParameterSet& parameters, const string& prefix);

  bool process(const DPBuffer& buffer) override;

  void finish() override { };

  void updateInfo(const DPInfo& info) override;

  void show(std::ostream&) const override;
  
private:
  DPInput* _input;
  string _name;
  std::vector<string> _directionStr;
  casacore::MDirection _direction;
  BeamCorrectionMode _mode;
};

} } // end namespaces

#endif

