// Timer.h: Accurate timer
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_TIMER_H
#define LOFAR_COMMON_TIMER_H

#include "PrettyUnits.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

#if defined __ia64__ && defined __INTEL_COMPILER
#include <ia64regs.h>
#endif

namespace dp3 {
namespace common {

/// \brief Very accurate timer for elapsed times.

/// Put timer.start() and timer.stop() calls around the piece of
/// code to be timed; make sure that start() and stop() calls alternate.
/// A timer can be started and stopped multiple times; both the average and
/// total time, as well as the number of iterations are printed.
/// The measured time is real time (as opposed to user or system time).
/// The timer can be used to measure from 1 nanosecond to a century interval.
/// The accuracy of the timer depends on the platform.
///
/// The internal RAII class StartStop can be used to do the start/stop.
/// The constructor starts the timer, while the destructor stops it. It has
/// the advantage that no explicit stop has to be given.
/// Moreover, it makes the start/stop exception-safe.
///
/// @tparam Clock Clock type, which meets the
/// <a href="https://en.cppreference.com/w/cpp/named_req/Clock">Clock named
/// requirements</a>. This argument allows using a mocked clock.
template <class Clock>
class BaseTimer {
 public:
  /// Constructor.
  /// @param name The name that will be used when printing the timer.
  explicit BaseTimer(const std::string& name = std::string())
      : name_(name), count_(0), duration_(std::chrono::seconds{0}) {}

  /// Starts the timer.
  void start() { start_ = Clock::now(); }
  /// Stops the timer
  void stop() {
    duration_ += Clock::now() - start_;
    ++count_;
  }

  /// Resets the timer to zero.
  void reset() {
    count_ = 0;
    duration_ = std::chrono::seconds{0};
  }

  /// Prints the timer.
  std::ostream& print(std::ostream& str) const {
    if (name_.empty()) {
      str << "timer: ";
    } else {
      str << std::left << std::setw(25) << name_ << ": " << std::right;
    }
    if (count_ == 0) {
      str << "not used";
    } else {
      const double total = getElapsed();
      // clang-format off
      str << "avg = " << PrettyTime(total / count_)
          << ", total = " << PrettyTime(total)
          << ", count = " << std::setw(9) << count_;
      // clang-format on
    }
    return str;
  }

  /// Get the elapsed time (in seconds).
  double getElapsed() const {
    return std::chrono::duration<double>{duration_}.count();
  }

  /// Get the average time (in seconds) between start/stop.
  double getAverage() const { return getElapsed() / getCount(); }

  /// Get the total number of times start/stop is done.
  uint64_t getCount() const { return count_; }

  /// Accumulate timer statistics.
  BaseTimer& operator+=(const BaseTimer& other) {
    duration_ += other.duration_;
    count_ += other.count_;

    return *this;
  }

  /// Internal class to do an automatic start/stop.
  class StartStop {
   public:
    StartStop(BaseTimer<Clock>& timer) : timer_(timer) { timer_.start(); }
    ~StartStop() { timer_.stop(); }

    /// Forbid copy.
    StartStop(const StartStop&) = delete;
    StartStop& operator=(const StartStop&) = delete;

   private:
    BaseTimer<Clock>& timer_;
  };

 private:
  std::string name_;

  /// The number of times the timing has been stopped.
  uint64_t count_;

  /// The total duration of all start/stop cycles.
  typename Clock::duration duration_;

  /// The timestamp of the most recent call to @rer start().
  typename Clock::time_point start_;
};

using NSTimer = BaseTimer<std::chrono::steady_clock>;

/**
 * Helper class to accumulate the execution time of multiple timers.
 *
 * Objects of this class measure duration of their life-time and add that time
 * to @a accumulator. The class is intended to be used in scope based fashion,
 * thus it's neither copyable nor movable.
 *
 * It's intended the class can be used to accumulate one timer from multiple
 * threads, so it should be able to use a @ref std::atomic<T>. Unfortunately
 * @ref std::atomic<T> for floating-point types requires C++20.
 *
 * The resolution of the class is based on the µs resolution of the @ref NSTimer
 * class.
 */
template <class T, class Clock = std::chrono::steady_clock>
#if __cplusplus > 201703L
requires requires(T& t, double d) {
  t += d;
}
#endif
class ScopedMicroSecondAccumulator {
 public:
  explicit ScopedMicroSecondAccumulator(T& accumulator)
      : accumulator_(accumulator) {
    timer_.start();
  }
  ScopedMicroSecondAccumulator(const ScopedMicroSecondAccumulator&) = delete;

  ~ScopedMicroSecondAccumulator() {
    timer_.stop();
    accumulator_ += timer_.getElapsed() * 1e6;
  }

 private:
  BaseTimer<Clock> timer_;
  T& accumulator_;
};

}  // namespace common
}  // namespace dp3

#endif
