// Step.h: Abstract base class for a DP3 step
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Class to hold code for virtual base class for Flaggers in DP3
/// @author Ger van Diepen

#ifndef DP3_STEPS_STEP_H_
#define DP3_STEPS_STEP_H_

#include <iosfwd>
#include <memory>
#include <string>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/BDABuffer.h>
#include <dp3/base/DPInfo.h>
#include <dp3/base/Direction.h>

#include <dp3/common/Fields.h>

namespace dp3 {
namespace steps {

/// @brief Abstract base class for a DP3 step

/// This class defines a step in the DP3 pipeline.
/// It is an abstract class from which all steps should be derived.
/// A few functions can or must be implemented. They are called by
/// the DP3 program in the following order.
/// <ul>
///  <li> 'updateInfo' should update its DPInfo object with the specific
///        step information. For example, in this way it is known
///       in all steps how the data are averaged and what the shape is.
///  <li> 'show' can be used to show the attributes.
///  <li> 'process' is called continuously to process the next time slot.
///        When processed, it should call 'process' of the next step.
///        When done (i.e. at the end of the input), it should return False.
///  <li> 'finish' finishes the processing which could mean that 'process'
///       of the next step has to be called several times. When done,
///       it should call 'finish' of the next step.
///  <li> 'addToMS' is called in the 'finish' step of tasks that write
///       or update a measurement set. It gives a step the opportunity
///       to add some data to the MS written/updated. It is, for example,
///       used by AOFlagger to write its statistics.
///  <li> 'showCounts' can be used to show possible counts of flags, etc.
/// </ul>
/// A Step object contains a DPInfo object telling the data settings for
/// a step (like channel info, baseline info, etc.).

class Step {
 public:
  typedef std::shared_ptr<Step> ShPtr;

  /// To check compatibility between steps before running.
  enum class MsType { kRegular, kBda };

  // Predefined fields values, that can be OR'd together:
  static constexpr dp3::common::Fields kDataField{
      dp3::common::Fields::Single::kData};
  static constexpr dp3::common::Fields kFlagsField{
      dp3::common::Fields::Single::kFlags};
  static constexpr dp3::common::Fields kWeightsField{
      dp3::common::Fields::Single::kWeights};
  static constexpr dp3::common::Fields kUvwField{
      dp3::common::Fields::Single::kUvw};

  /// Constructor to initialize.
  Step() : itsPrevStep(nullptr) {}

  /// Destructor.
  virtual ~Step();

  /// Process the data (legacy interface).
  /// When processed, it invokes the process function of the next step.
  /// It should return False at the end.
  /// A Step should either implement this legacy overload or the new overload
  /// that accepts a unique pointer.
  /// This base implementation tries calling the unique pointer overload,
  /// and throws if a Step does not implement that overload, too.
  virtual bool process(const base::DPBuffer& buffer);

  /// Process the data.
  /// When processed, the step should invoke the process function of the next
  /// step with the same buffer as argument.
  /// This base implementation tries calling the legacy process() overload,
  /// and throws if a Step does not implement that overload, too.
  /// When all Steps use this new overload, we should remove the legacy overload
  /// and this base implementation should always throw (similarly to the BDA
  /// overload below).
  /// @return False at the end of the input. True if there is more input.
  virtual bool process(std::unique_ptr<base::DPBuffer> buffer);

  /// Process the BDA data.
  /// When processed, it invokes the process function of the next step.
  /// It should return False at the end.
  virtual bool process(std::unique_ptr<base::BDABuffer>) {
    throw std::runtime_error("Step does not support BDA data processing.");
  }

  /// Finish the processing of this step and subsequent steps.
  virtual void finish() = 0;

  /// Set the info of this step and its next step.
  /// It calls the virtual function updateInfo to do the real work.
  /// It returns the info of the last step.
  const base::DPInfo& setInfo(const base::DPInfo&);

  /// Get the fields required by the current step.
  virtual dp3::common::Fields getRequiredFields() const = 0;

  /// Get the fields provided (modified and/or created) by the current step.
  /// The returned fields thus should not include (required) fields that are
  /// forwarded without modifications.
  virtual dp3::common::Fields getProvidedFields() const = 0;

  /// Get access to the info of the input.
  const base::DPInfo& getInfoIn() const { return itsInfoIn; }

  /// Get access to the info of the output.
  const base::DPInfo& getInfoOut() const { return itsInfoOut; }

  /// Legacy version of getInfoOut().
  /// New code should use getInfoOut() instead of getInfo()
  const base::DPInfo& getInfo() const { return itsInfoOut; }

  /// Show the step parameters.
  virtual void show(std::ostream&) const = 0;

  /// Show the flag counts if needed.
  /// The default implementation does nothing.
  virtual void showCounts(std::ostream&) const;

  /// Show the timings.
  /// The default implementation does nothing.
  virtual void showTimings(std::ostream&, double duration) const;

  /// Set the previous step.
  void setPrevStep(Step* prevStep) { itsPrevStep = prevStep; }

  /// Get the previous step.
  Step* getPrevStep() const { return itsPrevStep; }

  /// Set the next step.
  virtual void setNextStep(Step::ShPtr nextStep) {
    itsNextStep = nextStep;
    nextStep->setPrevStep(this);
  }

  /// Get the next step.
  const Step::ShPtr& getNextStep() const { return itsNextStep; }

  /// Return which datatype this step outputs.
  virtual MsType outputs() const { return MsType::kRegular; }

  /// Boolean if this step can process this type of data.
  virtual bool accepts(MsType dt) const { return dt == MsType::kRegular; }

 protected:
  /// @return Non-const reference to ouput info.
  base::DPInfo& infoOut() { return itsInfoOut; }

  /// Legacy version of infoOut().
  /// New code should use infoOut() instead of info().
  base::DPInfo& info() { return itsInfoOut; }

  /// Update the general info (called by setInfo).
  /// The default implementation copies the info.
  virtual void updateInfo(const base::DPInfo&);

  /// Add some data to the MeasurementSet written/updated.
  /// The default implementation only calls addToMS from the previous step
  virtual void addToMS(const std::string& msName);

 private:
  Step::ShPtr itsNextStep;
  Step* itsPrevStep;  /// Normal pointer for back links, prevent
                      /// two shared pointers to same object
  base::DPInfo itsInfoIn;
  base::DPInfo itsInfoOut;
  bool recursive_process_ = false;
};

/// Common interface for steps that produce model data.
class ModelDataStep : public Step {
 public:
  common::Fields getProvidedFields() const override { return kDataField; }

  /// @return The direction of the first patch.
  virtual base::Direction GetFirstDirection() const = 0;
};

}  // namespace steps
}  // namespace dp3

#endif  // DP3_STEPS_STEP_H_
