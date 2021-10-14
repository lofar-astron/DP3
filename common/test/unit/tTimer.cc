// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Timer.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <atomic>
#include <numeric>
#include <regex>
#include <thread>
#include <type_traits>

BOOST_AUTO_TEST_SUITE(timer)

namespace TestStaticAssert {
using ScopedTime = dp3::common::ScopedMicroSecondAccumulator<int>;
static_assert(!(std::is_copy_constructible<ScopedTime>::value ||
                std::is_move_constructible<ScopedTime>::value ||
                std::is_copy_assignable<ScopedTime>::value ||
                std::is_move_assignable<ScopedTime>::value),
              "The class is designed to be non-copyable and non-movable");
}  // namespace TestStaticAssert

static void Delay(dp3::common::NSTimer& timer, const int delay) {
  const dp3::common::NSTimer::StartStop scoped_time(timer);
  std::this_thread::sleep_for(std::chrono::milliseconds(delay));
}

BOOST_AUTO_TEST_CASE(timer_start_stop, *boost::unit_test::label("timer")) {
  dp3::common::NSTimer timer;
  {
    Delay(timer, 250);

    BOOST_CHECK_CLOSE(timer.getElapsed(), 0.25, 1);
    BOOST_CHECK_EQUAL(timer.getCount(), 1);
    BOOST_CHECK_EQUAL(timer.getElapsed(), timer.getAverage());
  }

  {
    Delay(timer, 50);

    BOOST_CHECK_CLOSE(timer.getElapsed(), 0.3, 1);
    BOOST_CHECK_EQUAL(timer.getCount(), 2);
    BOOST_CHECK_CLOSE(timer.getAverage(), 0.15, 1);
  }
  {
    Delay(timer, 300);

    BOOST_CHECK_CLOSE(timer.getElapsed(), 0.6, 1);
    BOOST_CHECK_EQUAL(timer.getCount(), 3);
    BOOST_CHECK_CLOSE(timer.getAverage(), 0.2, 1);
  }
}

static std::string print(const dp3::common::NSTimer& timer) {
  std::stringstream sstr;
  // Ensure the test doesn't depend on the system's locale settings.
  sstr.imbue(std::locale::classic());
  timer.print(sstr);
  return sstr.str();
}

BOOST_AUTO_TEST_CASE(timer_print, *boost::unit_test::label("timer")) {
  dp3::common::NSTimer timer;
  {
    Delay(timer, 250);
    BOOST_CHECK(std::regex_match(
        print(timer), std::regex{"timer: avg =  25[0-9] ms, total =  "
                                 "25[0-9] ms, count =         1"}));
  }

  {
    Delay(timer, 50);
    BOOST_CHECK(std::regex_match(
        print(timer), std::regex{"timer: avg =  1[5-6][0-9] ms, total =  "
                                 "3[0-1][0-9] ms, count =         2"}));
  }
  {
    Delay(timer, 300);
    BOOST_CHECK(std::regex_match(
        print(timer), std::regex{"timer: avg =  2[0-1][0-9] ms, total =  "
                                 "6[0-2][0-9] ms, count =         3"}));
  }
}

BOOST_AUTO_TEST_CASE(timer_print_unused, *boost::unit_test::label("timer")) {
  dp3::common::NSTimer timer;
  BOOST_CHECK_EQUAL(print(timer), "timer: not used");
}

BOOST_AUTO_TEST_CASE(timer_print_named, *boost::unit_test::label("timer")) {
  dp3::common::NSTimer timer{"foo"};
  Delay(timer, 100);
  BOOST_CHECK(std::regex_match(
      print(timer), std::regex{"foo                      : avg =  10[0-9] ms, "
                               "total =  10[0-9] ms, count =         1"}));
}

template <class T>
static void Delay(T& accumulator, const int delay) {
  using ScopedTime = dp3::common::ScopedMicroSecondAccumulator<T>;
  const ScopedTime scoped_time(accumulator);
  std::this_thread::sleep_for(std::chrono::milliseconds(delay));
}

BOOST_AUTO_TEST_CASE(scoped_micro_second_accumumator_start_zero,
                     *boost::unit_test::label("timer")) {
  int accumulator = 0;
  Delay(accumulator, 1000);

  // Note the timer is known to be not too accurate.
  BOOST_CHECK_CLOSE(accumulator, 1e6, 5);
}

BOOST_AUTO_TEST_CASE(scoped_micro_second_accumumator_start_non_zero,
                     *boost::unit_test::label("timer")) {
  int accumulator = 10000;
  Delay(accumulator, 1);

  // Note the timer is known to be not too accurate.
  BOOST_CHECK_CLOSE(accumulator, 11e3, 5);
}

// Uses const int* since that's the iterator type of a std::aray<int, N>.
static int RunAllDelaysInThreads(const int* first, const int* last) {
  std::atomic<int> accumulator{0};
  std::vector<std::thread> threads;
  std::for_each(first, last, [&](const int delay) {
    threads.emplace_back([&accumulator, delay] { Delay(accumulator, delay); });
  });
  std::for_each(threads.begin(), threads.end(),
                [](std::thread& thread) { thread.join(); });

  return accumulator;
}

BOOST_AUTO_TEST_CASE(scoped_micro_second_accumumator_atomic_integral,
                     *boost::unit_test::label("timer")) {
  constexpr std::array<int, 12> delays{5, 5, 3, 4, 2, 1, 1, 3, 4, 2, 2, 1};
  const double expected =
      1e3 * std::accumulate(delays.begin(), delays.end(), 0);

  const double result = RunAllDelaysInThreads(delays.begin(), delays.end());

  // Note the timer is known to be not too accurate.
  BOOST_CHECK_CLOSE(result, expected, 5);
}

BOOST_AUTO_TEST_SUITE_END()
