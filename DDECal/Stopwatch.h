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
