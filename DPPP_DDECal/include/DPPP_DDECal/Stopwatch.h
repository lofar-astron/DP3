#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <string>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/**
	@author A.R. Offringa <offringa@astro.rug.nl>
*/
class Stopwatch{
	public:
		Stopwatch();
		explicit Stopwatch(bool start);
		~Stopwatch();
		void Start();
		void Pause();
		void Reset();

		std::string ToString() const;
		std::string ToShortString() const;
		long double Seconds() const;
		
		std::string ToDaysString() const;
		std::string ToHoursString() const;
		std::string ToMinutesString() const;
		std::string ToSecondsString() const;
		std::string ToMilliSecondsString() const;
		std::string ToMicroSecondsString() const;
		std::string ToNanoSecondsString() const;
	private:
		bool _running;
		boost::posix_time::ptime _startTime;
		boost::posix_time::time_duration _sum;
};

// end of Doxygen skip
#endif

#endif
