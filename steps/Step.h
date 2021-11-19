// Step.h: Abstract base class for a DPPP step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Class to hold code for virtual base class for Flaggers in DPPP
/// @author Ger van Diepen

#ifndef DPPP_DPSTEP_H
#define DPPP_DPSTEP_H

#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"
#include "../base/Direction.h"

#include "../common/Timer.h"

#include <iosfwd>
#include <memory>

namespace dp3 {
namespace base {
class BDABuffer;
}
namespace steps {

/// @brief Abstract base class for a DPPP step

/// This class defines a step in the DPPP pipeline.
/// It is an abstract class from which all steps should be derived.
/// A few functions can or must be implemented. They are called by
/// the NDPPP program in the following order.
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
///  <li> 'addToMS' is called after 'finish'. It gives a step the opportunity
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

  /// Constructor to initialize.
  Step() : itsPrevStep(0) {}

  /// Destructor.
  virtual ~Step();

  /// Process the data.
  /// When processed, it invokes the process function of the next step.
  /// It should return False at the end.
  virtual bool process(const base::DPBuffer&) {
    throw std::runtime_error("Step does not support regular data processing.");
  }

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

  /// Get access to the info.
  const base::DPInfo& getInfo() const { return itsInfo; }

  /// Add some data to the MeasurementSet written/updated.
  /// The default implementation only calls addToMS from the previous step
  virtual void addToMS(const string& msName);

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

  /// True when the step modifies any data (or flags, meta data, etc)
  /// and therefore requires the step to be followed by a write step.
  virtual bool modifiesData() const { return true; }

 protected:
  base::DPInfo& info() { return itsInfo; }

  /// Update the general info (called by setInfo).
  /// The default implementation copies the info.
  virtual void updateInfo(const base::DPInfo&);

 private:
  Step::ShPtr itsNextStep;
  Step* itsPrevStep;  /// Normal pointer for back links, prevent
                      /// two shared pointers to same object
  base::DPInfo itsInfo;
};

/// @brief This class defines a null step in the DPPP pipeline.
/// It can be used as the last step in the pipeline, so other steps
/// do not need to test if there is a next step.

class NullStep : public Step {
 public:
  virtual ~NullStep();

  /// Process the data. It does nothing.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  /// It does nothing.
  virtual void finish();

  /// Show the step parameters.
  /// It does nothing.
  virtual void show(std::ostream&) const;
};

/// @brief This class defines step in the DPPP pipeline that keeps the result
/// to make it possible to get the result of another step.
/// It keeps the result and calls process of the next step.

class ResultStep : public Step {
 public:
  /// Create the object. By default it sets its next step to the NullStep.
  ResultStep();

  virtual ~ResultStep();

  /// Keep the buffer.
  virtual bool process(const base::DPBuffer&);

  /// Finish does not do anything.
  virtual void finish();

  /// Show the step parameters.
  /// It does nothing.
  virtual void show(std::ostream&) const;

  /// Get the result.
  const base::DPBuffer& get() const { return itsBuffer; }
  base::DPBuffer& get() { return itsBuffer; }

  /// Clear the buffer.
  void clear() { itsBuffer = base::DPBuffer(); }

 private:
  base::DPBuffer itsBuffer;
};

/// @brief This class defines step in the DPPP pipeline that keeps the result
/// to make it possible to get the result of another step.
/// It keeps the result and calls process of the next step.
/// Buffers are accumulated until cleared.

class MultiResultStep : public Step {
 public:
  /// Create the object. By default it sets its next step to the NullStep.
  MultiResultStep(unsigned int size);

  virtual ~MultiResultStep();

  /// Add the buffer to the vector of kept buffers.
  virtual bool process(const base::DPBuffer&);

  /// Finish does not do anything.
  virtual void finish();

  /// Show the step parameters.
  /// It does nothing.
  virtual void show(std::ostream&) const;

  /// Get the result.
  const std::vector<base::DPBuffer>& get() const { return itsBuffers; }
  std::vector<base::DPBuffer>& get() { return itsBuffers; }

  /// Get the size of the result.
  size_t size() const { return itsSize; }

  /// Clear the buffers.
  void clear() { itsSize = 0; }

 private:
  std::vector<base::DPBuffer> itsBuffers;
  size_t itsSize;
};

/// Common interface for steps that produce model data.
class ModelDataStep : public Step {
 public:
  /// @return The direction of the first patch.
  virtual base::Direction GetFirstDirection() const = 0;
};

}  // namespace steps
}  // namespace dp3

#endif
