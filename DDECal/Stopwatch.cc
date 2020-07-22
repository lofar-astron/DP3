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

#include "Stopwatch.h"

#include <cmath>
#include <sstream>

#include <boost/date_time/posix_time/posix_time.hpp>

Stopwatch::Stopwatch() : _running(false), _sum(boost::posix_time::seconds(0)) {}

Stopwatch::Stopwatch(bool start)
    : _running(start),
      _startTime(boost::posix_time::microsec_clock::local_time()),
      _sum(boost::posix_time::seconds(0)) {}

Stopwatch::~Stopwatch() {}

void Stopwatch::Start() {
  if (!_running) {
    _startTime = boost::posix_time::microsec_clock::local_time();
    _running = true;
  }
}

void Stopwatch::Pause() {
  if (_running) {
    _sum += (boost::posix_time::microsec_clock::local_time() - _startTime);
    _running = false;
  }
}

void Stopwatch::Reset() {
  _running = false;
  _sum = boost::posix_time::seconds(0);
}

std::string Stopwatch::ToString() const {
  if (_running) {
    boost::posix_time::time_duration current =
        _sum + (boost::posix_time::microsec_clock::local_time() - _startTime);
    return to_simple_string(current);
  } else {
    return to_simple_string(_sum);
  }
}

std::string Stopwatch::ToShortString() const {
  const long double seconds = Seconds();
  if (seconds >= 60 * 60 * 24)
    return ToDaysString();
  else if (seconds >= 60 * 60)
    return ToHoursString();
  else if (seconds >= 60)
    return ToMinutesString();
  else if (seconds >= 1.0)
    return ToSecondsString();
  else if (seconds >= 0.001)
    return ToMilliSecondsString();
  else if (seconds >= 0.000001)
    return ToMicroSecondsString();
  else
    return ToNanoSecondsString();
}

std::string Stopwatch::ToDaysString() const {
  const long double days = roundl(Seconds() / (60.0 * 60.0)) / 24.0;
  std::stringstream str;
  if (days >= 10.0)
    str << roundl(days) << " days";
  else
    str << floorl(days) << 'd' << (days * 24.0) << 'h';
  return str.str();
}

std::string Stopwatch::ToHoursString() const {
  const long double hours = roundl(Seconds() / 60.0) / 60.0;
  std::stringstream str;
  if (hours >= 10.0)
    str << roundl(hours) << 'h';
  else
    str << floorl(hours) << 'h' << (hours * 60.0) << 'm';
  return str.str();
}

std::string Stopwatch::ToMinutesString() const {
  const long double mins = roundl(Seconds()) / 60.0;
  std::stringstream str;
  if (mins >= 10.0)
    str << roundl(mins) << " min";
  else
    str << floorl(mins) << 'm' << fmod(mins * 60.0, 60.0) << 's';
  return str.str();
}

std::string Stopwatch::ToSecondsString() const {
  const long double seconds = roundl(Seconds() * 10.0) / 10.0;
  std::stringstream str;
  if (seconds >= 10.0)
    str << roundl(Seconds()) << 's';
  else
    str << seconds << 's';
  return str.str();
}

std::string Stopwatch::ToMilliSecondsString() const {
  const long double msec = roundl(Seconds() * 10000.0) / 10.0;
  std::stringstream str;
  if (msec >= 10.0)
    str << roundl(Seconds() * 1000.0) << "ms";
  else
    str << msec << "ms";
  return str.str();
}

std::string Stopwatch::ToMicroSecondsString() const {
  const long double usec = roundl(Seconds() * 10000000.0) / 10.0;
  std::stringstream str;
  if (usec >= 10.0)
    str << roundl(Seconds() * 1000000.0) << "µs";
  else
    str << usec << "µs";
  return str.str();
}

std::string Stopwatch::ToNanoSecondsString() const {
  const long double nsec = roundl(Seconds() * 10000000000.0) / 10.0;
  std::stringstream str;
  if (nsec >= 10.0)
    str << roundl(Seconds() * 1000000000.0) << "ns";
  else
    str << nsec << "ns";
  return str.str();
}

long double Stopwatch::Seconds() const {
  if (_running) {
    boost::posix_time::time_duration current =
        _sum + (boost::posix_time::microsec_clock::local_time() - _startTime);
    return (long double)current.total_microseconds() / 1000000.0;
  } else {
    return (long double)_sum.total_microseconds() / 1000000.0;
  }
}
