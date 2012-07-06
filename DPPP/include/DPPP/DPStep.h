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
// @brief Class to hold code for virtual base class for Flaggers in IDPPP

#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <Common/lofar_smartptr.h>
#include <Common/Timer.h>
#include <iosfwd>

namespace LOFAR {
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
      typedef shared_ptr<DPStep> ShPtr;

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
      // The default implementation does nothing.
      virtual void addToMS (const string& msName);

      // Show the step parameters.
      virtual void show (std::ostream&) const = 0;

      // Show the flag counts if needed.
      // The default implementation does nothing.
      virtual void showCounts (std::ostream&) const;

      // Show the timings.
      // The default implementation does nothing.
      virtual void showTimings (std::ostream&, double duration) const;

      // Set the next step.
      void setNextStep (const DPStep::ShPtr& nextStep)
        { itsNextStep = nextStep; }

      // Get the next step.
      const DPStep::ShPtr& getNextStep() const
        { return itsNextStep; }

    protected:
      DPInfo& info()
        { return itsInfo; }

    private:
      // Update the general info (called by setInfo).
      // The default implementation copies the info.
      virtual void updateInfo (const DPInfo&);

      //# Data mamabers.
      DPStep::ShPtr itsNextStep;
      DPInfo        itsInfo;
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
    // It only keeps the buffer, but does not process it in next steps.

    class ResultStep: public DPStep
    {
    public:
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
        { itsBuffer.clear(); }

    private:
      DPBuffer itsBuffer;
    };



    // @ingroup NDPPP

    // This class defines step in the DPPP pipeline that keeps the result
    // to make it possible to get the result of another step.
    // It only keeps buffers, but does not process them in next steps.
    // Buffers are accumulated until cleared.

    class MultiResultStep: public DPStep
    {
    public:
      // Create the object. By default it sets its next step to the NullStep.
      MultiResultStep (uint reserveSize);

      virtual ~MultiResultStep();

      // Add the buffer to the vector of kept buffers.
      virtual bool process (const DPBuffer&);

      // Finish does not do anything.
      virtual void finish();

      // Show the step parameters.
      // It does nothing.
      virtual void show (std::ostream&) const;

      // Get the result.
      const vector<DPBuffer>& get() const
        { return itsBuffers; }
      vector<DPBuffer>& get()
        { return itsBuffers; }

      // Clear the buffers.
      void clear()
        { itsBuffers.clear(); }

    private:
      vector<DPBuffer> itsBuffers;
    };

  } //# end namespace
}

#endif
