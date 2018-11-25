//# DPStep.h: Abstract base class for a DPPP step
//# Copyright (C) 2010
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$
//#
//# @author Ger van Diepen

#ifndef DPPP_DPSTEP_H
#define DPPP_DPSTEP_H

// @file
// @brief Class to hold code for virtual base class for Flaggers in DPPP

#include "DPBuffer.h"
#include "DPInfo.h"

#include "../Common/Timer.h"

#include <iosfwd>
#include <memory>

namespace DP3 {
  namespace DPPP {

    // @ingroup NDPPP

    // This class defines a step in the DPPP pipeline.
    // It is an abstract class from which all steps should be derived.
    // A few functions can or must be implemented. They are called by
    // the NDPPP program in the following order.
    // <ul>
    //  <li> 'updateInfo' should update its DPInfo object with the specific
    //        step information. For example, in this way it is known
    //       in all steps how the data are averaged and what the shape is.
    //  <li> 'show' can be used to show the attributes.
    //  <li> 'process' is called continuously to process the next time slot.
    //        When processed, it should call 'process' of the next step.
    //        When done (i.e. at the end of the input), it should return False.
    //  <li> 'finish' finishes the processing which could mean that 'process'
    //       of the next step has to be called several times. When done,
    //       it should call 'finish' of the next step.
    //  <li> 'addToMS' is called after 'finish'. It gives a step the opportunity
    //       to add some data to the MS written/updated. It is, for example,
    //       used by AOFlagger to write its statistics.
    //  <li> 'showCounts' can be used to show possible counts of flags, etc.
    // </ul>
    // A DPStep object contains a DPInfo object telling the data settings for
    // a step (like channel info, baseline info, etc.).

    class DPStep
    {
    public:
      // Define the shared pointer for this type.
      typedef std::shared_ptr<DPStep> ShPtr;

      // Constructor to initialize.
      DPStep()
        : itsPrevStep(0)
      {}

      // Destructor.
      virtual ~DPStep();

      // Process the data.
      // When processed, it invokes the process function of the next step.
      // It should return False at the end.
      virtual bool process (const DPBuffer&) = 0;

      // Finish the processing of this step and subsequent steps.
      virtual void finish() = 0;

      // Set the info of this step and its next step.
      // It calls the virtual function updateInfo to do the real work.
      // It returns the info of the last step.
      const DPInfo& setInfo (const DPInfo&);

      // Get access to the info.
      const DPInfo& getInfo() const
        { return itsInfo; }

      // Add some data to the MeasurementSet written/updated.
      // The default implementation only calls addToMS from the previous step
      virtual void addToMS (const string& msName);

      // Show the step parameters.
      virtual void show (std::ostream&) const = 0;

      // Show the flag counts if needed.
      // The default implementation does nothing.
      virtual void showCounts (std::ostream&) const;

      // Show the timings.
      // The default implementation does nothing.
      virtual void showTimings (std::ostream&, double duration) const;

      // Set the previous step.
      void setPrevStep (DPStep* prevStep)
      { itsPrevStep = prevStep; }

      // Get the previous step.
      DPStep* getPrevStep () const
      { return itsPrevStep; }

      // Set the next step.
      virtual void setNextStep (DPStep::ShPtr nextStep)
        { itsNextStep = nextStep;
          nextStep->setPrevStep(this);
        }

      // Get the next step.
      const DPStep::ShPtr& getNextStep() const
        { return itsNextStep; }
        
    protected:
      DPInfo& info()
        { return itsInfo; }

      // Update the general info (called by setInfo).
      // The default implementation copies the info.
      virtual void updateInfo (const DPInfo&);

    private:
      //# Data members.
      DPStep::ShPtr itsNextStep;
      DPStep* itsPrevStep; // Normal pointer for back links, prevent
                           // two shared pointers to same object
      DPInfo itsInfo;
    };



    // @ingroup NDPPP

    // This class defines a null step in the DPPP pipeline.
    // It can be used as the last step in the pipeline, so other steps
    // do not need to test if there is a next step.

    class NullStep: public DPStep
    {
    public:
      virtual ~NullStep();

      // Process the data. It does nothing.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      // It does nothing.
      virtual void finish();

      // Show the step parameters.
      // It does nothing.
      virtual void show (std::ostream&) const;
    };



    // @ingroup NDPPP

    // This class defines step in the DPPP pipeline that keeps the result
    // to make it possible to get the result of another step.
    // It keeps the result and calls process of the next step.

    class ResultStep: public DPStep
    {
    public:
      typedef std::shared_ptr<ResultStep> ShPtr;
      // Create the object. By default it sets its next step to the NullStep.
      ResultStep();

      virtual ~ResultStep();

      // Keep the buffer.
      virtual bool process (const DPBuffer&);

      // Finish does not do anything.
      virtual void finish();

      // Show the step parameters.
      // It does nothing.
      virtual void show (std::ostream&) const;

      // Get the result.
      const DPBuffer& get() const
        { return itsBuffer; }
      DPBuffer& get()
        { return itsBuffer; }

      // Clear the buffer.
      void clear()
      { itsBuffer = DPBuffer(); }

    private:
      DPBuffer itsBuffer;
    };



    // @ingroup NDPPP

    // This class defines step in the DPPP pipeline that keeps the result
    // to make it possible to get the result of another step.
    // It keeps the result and calls process of the next step.
    // Buffers are accumulated until cleared.

    class MultiResultStep: public DPStep
    {
    public:
      // Define the shared pointer for this type.
      typedef std::shared_ptr<MultiResultStep> ShPtr;

      // Create the object. By default it sets its next step to the NullStep.
      MultiResultStep (uint size);

      virtual ~MultiResultStep();

      // Add the buffer to the vector of kept buffers.
      virtual bool process (const DPBuffer&);

      // Finish does not do anything.
      virtual void finish();

      // Show the step parameters.
      // It does nothing.
      virtual void show (std::ostream&) const;

      // Get the result.
      const std::vector<DPBuffer>& get() const
        { return itsBuffers; }
      std::vector<DPBuffer>& get()
        { return itsBuffers; }

      // Get the size of the result.
      size_t size() const
        { return itsSize; }

      // Clear the buffers.
      void clear()
        { itsSize = 0; }

    private:
      std::vector<DPBuffer> itsBuffers;
      size_t           itsSize;
    };

  } //# end namespace
}

#endif
