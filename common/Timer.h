// Timer.h: Accurate timer
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_TIMER_H
#define LOFAR_COMMON_TIMER_H

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>

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
/// The internal class NSTimer::StartStop can be used to do the start/stop.
/// The constructor starts the timer, while the destructor stops it. It has
/// the advantage that no explicit stop has to be given.
/// Moreover, it makes the start/stop exception-safe.
class NSTimer {
 public:
  /// Construct.
  /// The given name will be used when the timer is printed.
  NSTimer(const std::string& name = std::string(),
          bool print_on_destruction = false, bool log_on_destruction = false);

  /// Destruct.
  /// The time is printed on stderr if print_on_destruction is true.
  ~NSTimer();

  /// Start the timer.
  void start();
  /// Stop the timer
  void stop();

  /// Reset the timer to zero.
  void reset();

  /// Print the timer.
  /// @{
  std::ostream& print(std::ostream&) const;
  friend std::ostream& operator<<(std::ostream&, const NSTimer&);
  /// @}

  /// Get the elapsed time (in seconds).
  double getElapsed() const {
    return std::chrono::duration<double>{duration_}.count();
  }

  /// Get the average time (in seconds) between start/stop.
  double getAverage() const;

  /// Get the total number of times start/stop is done.
  uint64_t getCount() const;

  /// Accumulate timer statistics.
  NSTimer& operator+=(const NSTimer& other);

  /// Internal class to do an automatic start/stop.
  class StartStop {
   public:
    StartStop(NSTimer& timer) : itsTimer(timer) { itsTimer.start(); }
    ~StartStop() { itsTimer.stop(); }

   private:
    /// Forbid copy.
    StartStop(const StartStop&);
    StartStop& operator=(StartStop&);
    NSTimer& itsTimer;
  };

 private:
  std::string name_;
  bool print_on_destruction_;
  bool log_on_destruction_;

  /// The number of times the timing has been stopped.
  uint64_t count_{0};

  /// The total duration of all start/stop cycles.
  std::chrono::steady_clock::duration duration_;

  /// The timestamp of the most recent call to @rer start().
  std::chrono::steady_clock::time_point start_;
};

inline void NSTimer::reset() {
  count_ = 0;
  duration_ = std::chrono::seconds{0};
}

inline double NSTimer::getAverage() const { return getElapsed() / getCount(); }

inline uint64_t NSTimer::getCount() const { return count_; }

inline NSTimer& NSTimer::operator+=(const NSTimer& other) {
  duration_ += other.duration_;
  count_ += other.count_;

  return *this;
}

inline NSTimer::NSTimer(const std::string& name, bool print_on_destruction,
                        bool log_on_destruction)
    : name_(name),
      print_on_destruction_(print_on_destruction),
      log_on_destruction_(log_on_destruction) {
  reset();
}

inline NSTimer::~NSTimer() {
  if (print_on_destruction_) {
    if (log_on_destruction_) {
      std::stringstream logStr;
      print(logStr);
      std::clog << logStr.str() << '\n';
    } else
      std::clog << *this << std::endl;
  }
}

inline std::ostream& operator<<(std::ostream& str, const NSTimer& timer) {
  return timer.print(str);
}

inline void NSTimer::start() { start_ = std::chrono::steady_clock::now(); }

inline void NSTimer::stop() {
  duration_ += std::chrono::steady_clock::now() - start_;
  ++count_;
}

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
template <class T>
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
  NSTimer timer_;
  T& accumulator_;
};

}  // namespace common
}  // namespace dp3

#endif
