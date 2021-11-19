// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_DPRUN_H
#define DPPP_DPRUN_H

#include "../steps/Step.h"
#include "../steps/MSReader.h"

#include <map>

namespace dp3 {
namespace common {
class ParameterSet;
}
namespace base {

/// @brief Class to run steps like averaging and flagging on an MS

/// This class contains a single static function that creates and executes
/// the steps defined in the parset file.
/// The parset file is documented on the LOFAR wiki.

class DP3 {
 public:
  /// Define the function to create a step from the given parameterset.
  typedef steps::Step::ShPtr StepCtor(steps::InputStep*,
                                      const common::ParameterSet&,
                                      const std::string& prefix);

  /// Add a function creating a Step to the map.
  static void registerStepCtor(const std::string&, StepCtor*);

  /// Create a step object from the given parameters.
  /// It looks up the step type in theirStepMap. If not found, it will
  /// try to load a shared library with that name and execute the
  /// register function in it.
  static StepCtor* findStepCtor(const std::string& type);

  /// Execute the steps defined in the parset file.
  /// Possible parameters given at the command line are taken into account.
  static void execute(const std::string& parsetName, int argc = 0,
                      char* argv[] = 0);

  /// Create a chain of step objects that are connected together.
  /// A writer will be added to the steps if it is not defined,
  /// and a terminating NullStep is added.
  static std::shared_ptr<steps::InputStep> makeMainSteps(
      const common::ParameterSet& parset);

  /// Create a chain of step objects that are connected together.
  /// Unlike makeMainSteps(), this does neither add a writer nor
  /// terminate the chain with a NullStep.
  /// @param step_names_key Parset key for getting the step names.
  static steps::Step::ShPtr makeStepsFromParset(
      const common::ParameterSet& parset, const std::string& prefix,
      const std::string& step_names_key, steps::InputStep& inputStep,
      bool terminateChain, steps::Step::MsType initial_step_output);

  static steps::Step::ShPtr makeSingleStep(const std::string& type,
                                           steps::InputStep* inputStep,
                                           const common::ParameterSet& parset,
                                           const std::string& prefix,
                                           std::string& msName,
                                           steps::Step::MsType inputType);

 private:
  /// Create an output step, either an MSWriter, MSUpdater or an MSBDAWriter
  /// If no data are modified (for example if only count was done),
  /// still an MSUpdater is created, but it will not write anything.
  /// It reads the output name from the parset. If the prefix is "", it
  /// reads msout or msout.name, otherwise it reads name from the output step
  /// Create an updater step if an input MS was given; otherwise a writer.
  /// Create an updater step only if needed (e.g. not if only count is done).
  /// If the user specified an output MS name, a writer or updater is always
  /// created If there is a writer, the reader needs to read the visibility
  /// data. reader should be the original reader
  static steps::Step::ShPtr makeOutputStep(steps::InputStep* reader,
                                           const common::ParameterSet& parset,
                                           const std::string& prefix,
                                           std::string& currentMSName,
                                           steps::Step::MsType inputType);

  /// The map to create a step object from its type name.
  static std::map<std::string, StepCtor*> theirStepMap;
};

}  // namespace base
}  // namespace dp3

#endif
