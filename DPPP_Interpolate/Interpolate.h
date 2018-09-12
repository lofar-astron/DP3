#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <deque>
#include <string>

#include "../DPPP/DPInput.h"
#include "../DPPP/DPBuffer.h"

#include "../Common/Lane.h"
#include "../Common/ParameterSet.h"

#include <casacore/casa/Arrays/Cube.h>

extern "C" void register_interpolate();

namespace DP3 { namespace DPPP {
		
	class Interpolate : public DPStep
	{
	public:
		// Construct the object.
		// Parameters are obtained from the parset using the given prefix.
		Interpolate (DPInput*, const ParameterSet&, const string& prefix);

		virtual ~Interpolate() = default;

		// Process the data.
		// It keeps the data.
		// When processed, it invokes the process function of the next step.
		virtual bool process (const DPBuffer&);

		// Finish the processing of this step and subsequent steps.
		virtual void finish();

		// Update the general info.
		virtual void updateInfo (const DPInfo&);

		// Show the step parameters.
		virtual void show (std::ostream&) const;

		// Show the timings.
		virtual void showTimings (std::ostream&, double duration) const;

		static DPStep::ShPtr makeStep(DPInput* input, const ParameterSet& parset, const std::string& prefix);
		
	private:
		void interpolateTimestep(size_t index);
		void interpolateSample(size_t timestep, size_t baseline, size_t channel, size_t pol);
		void sendFrontBufferToNextStep();
		void interpolationThread();
		
		struct Sample {
			Sample() = default;
			Sample(size_t timestep_, size_t baseline_, size_t channel_, size_t pol_) :
				timestep(timestep_), baseline(baseline_), channel(channel_), pol(pol_)
				{ }
			size_t timestep;
			size_t baseline;
			size_t channel;
			size_t pol;
		};
		
		//# Data members.
		std::string _name;
		size_t _interpolatedPos;
		std::deque<DPBuffer> _buffers;
		size_t _windowSize;
		NSTimer _timer;
		ao::lane<Sample> _lane;
		std::vector<float> _kernelLookup;
	};

} } //# end namespace

#endif
