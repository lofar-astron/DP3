// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_DP3_H_
#define DP3_BASE_DP3_H_

#include "../common/ParameterSet.h"
#include "../steps/InputStep.h"

namespace dp3 {
namespace base {

void ShowUsage();

/// Command-line interface
void ExecuteFromCommandLine(int argc, char* argv[]);

/// Execute the steps defined in the parset file.
/// Possible parameters given at the command line are taken into account.
void Execute(const std::string& parsetName, int argc = 0,
             char* argv[] = nullptr);

/// Create a step
/// @param type Type of the step.
/// @param parset ParameterSet containing the configuration for the step.
/// @param prefix Prefix, including trailing dot ("."), to use for looking up
/// parameters in the ParameterSet.
/// @param input_type Type of input data, BDA or regular.
/// @return Pointer to the newly created step, or a null pointer if the
/// type not recognized.
std::shared_ptr<steps::Step> MakeSingleStep(const std::string& type,
                                            const common::ParameterSet& parset,
                                            const std::string& prefix,
                                            steps::Step::MsType input_type);

/// Create a chain of step objects that are connected together.
/// A writer will be added to the steps if it is not defined,
/// and a terminating NullStep is added.
std::shared_ptr<steps::InputStep> MakeMainSteps(
    const common::ParameterSet& parset);

/// Create a chain of step objects that are connected together.
/// Unlike MakeMainSteps(), this does neither add a writer nor
/// terminate the chain with a NullStep.
/// @param step_names_key Parset key for getting the step names.
/// @return The first step of the step chain or a null pointer if the chain is
/// empty.
std::shared_ptr<steps::Step> MakeStepsFromParset(
    const common::ParameterSet& parset, const std::string& prefix,
    const std::string& step_names_key, const std::string& input_ms_name,
    bool terminateChain, steps::Step::MsType initial_step_output);

/// Go through the steps to define which fields of the measurement set should
/// be read by the input step.
/// @param first_step The first step of a chain of steps.
/// @return The combined required fields of the step chain.
dp3::common::Fields GetChainRequiredFields(
    std::shared_ptr<steps::Step> first_step);

/// Go through a step chain, combine provided fields of non-output steps and
/// call SetFieldsToWrite for all output steps.
/// @param first_step The first step of a chain of steps.
/// @param provided_fields The provided fields before the step chain.
/// Steps that have sub-steps should use this argument.
/// @return The remaining provided fields after the last step in the chain.
dp3::common::Fields SetChainProvidedFields(
    std::shared_ptr<steps::Step> first_step,
    dp3::common::Fields provided_fields = dp3::common::Fields());

/// Get the number of threads used in DP3. This function is needed since
/// each library / application has its own aocommon Threadpool object.
size_t GetNThreads();

/// Set the number of threads used in DP3. This function is needed since
/// each library / application has its own aocommon Threadpool object.
void SetNThreads(size_t n_threads);

}  // namespace base
}  // namespace dp3

#endif  // DP3_BASE_DP3_H_
