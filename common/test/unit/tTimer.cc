// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Timer.h"

#include <boost/test/unit_test.hpp>

#include <deque>
#include <vector>

namespace {

/// Global queue with the time points that MockClock::now() should return.
static std::deque<std::chrono::system_clock::time_point> mock_time_points;

struct MockClock {
  using duration = std::chrono::system_clock::duration;
  using time_point = std::chrono::system_clock::time_point;

  static time_point now() {
    BOOST_REQUIRE(!mock_time_points.empty());
    time_point t = mock_time_points.front();
    mock_time_points.pop_front();
    return t;
  }
};

using TestTimer = dp3::common::BaseTimer<MockClock>;

// Epsilon value for comparing double values in BOOST_CHECK_CLOSE.
const double kEpsilon = 1e-6;

}  // namespace

BOOST_AUTO_TEST_SUITE(timer)

namespace TestStaticAssert {
using ScopedTime = dp3::common::ScopedMicroSecondAccumulator<int>;
static_assert(!(std::is_copy_constructible<ScopedTime>::value ||
                std::is_move_constructible<ScopedTime>::value ||
                std::is_copy_assignable<ScopedTime>::value ||
                std::is_move_assignable<ScopedTime>::value),
              "The class is designed to be non-copyable and non-movable");
}  // namespace TestStaticAssert

static void Delay(TestTimer& timer, const int delay) {
  const MockClock::time_point now = MockClock::time_point::clock::now();
  mock_time_points.push_back(now);
  mock_time_points.push_back(now + std::chrono::milliseconds(delay));
  { const TestTimer::StartStop scoped_time(timer); }
  BOOST_CHECK(mock_time_points.empty());
}

BOOST_AUTO_TEST_CASE(timer_constructor) {
  TestTimer timer;
  BOOST_CHECK_EQUAL(timer.getElapsed(), 0.0);
  BOOST_CHECK_EQUAL(timer.getCount(), 0);
}

BOOST_AUTO_TEST_CASE(timer_start_stop) {
  TestTimer timer;
  {
    Delay(timer, 250);

    BOOST_CHECK_CLOSE(timer.getElapsed(), 0.25, kEpsilon);
    BOOST_CHECK_EQUAL(timer.getCount(), 1);
    BOOST_CHECK_EQUAL(timer.getElapsed(), timer.getAverage());
  }

  {
    Delay(timer, 50);

    BOOST_CHECK_CLOSE(timer.getElapsed(), 0.3, kEpsilon);
    BOOST_CHECK_EQUAL(timer.getCount(), 2);
    BOOST_CHECK_CLOSE(timer.getAverage(), 0.15, kEpsilon);
  }
  {
    Delay(timer, 300);

    BOOST_CHECK_CLOSE(timer.getElapsed(), 0.6, kEpsilon);
    BOOST_CHECK_EQUAL(timer.getCount(), 3);
    BOOST_CHECK_CLOSE(timer.getAverage(), 0.2, kEpsilon);
  }
}

BOOST_AUTO_TEST_CASE(timer_add) {
  TestTimer timer1;
  TestTimer timer2;
  Delay(timer1, 200);
  Delay(timer2, 50);
  timer1 += timer2;
  BOOST_CHECK_CLOSE(timer1.getElapsed(), 0.25, kEpsilon);
  BOOST_CHECK_EQUAL(timer1.getCount(), 2);
  BOOST_CHECK_CLOSE(timer1.getAverage(), 0.125, kEpsilon);
}

static std::string Print(const TestTimer& timer) {
  std::stringstream sstr;
  // Ensure the test doesn't depend on the system's locale settings.
  sstr.imbue(std::locale::classic());
  timer.print(sstr);
  return sstr.str();
}

BOOST_AUTO_TEST_CASE(timer_print) {
  TestTimer timer;
  {
    Delay(timer, 250);
    BOOST_CHECK_EQUAL(
        Print(timer),
        "timer: avg =  250 ms, total =  250 ms, count =         1");
  }
  {
    Delay(timer, 50);
    BOOST_CHECK_EQUAL(
        Print(timer),
        "timer: avg =  150 ms, total =  300 ms, count =         2");
  }
  {
    Delay(timer, 300);
    BOOST_CHECK_EQUAL(
        Print(timer),
        "timer: avg =  200 ms, total =  600 ms, count =         3");
  }
}

BOOST_AUTO_TEST_CASE(timer_print_unused) {
  TestTimer timer;
  BOOST_CHECK_EQUAL(Print(timer), "timer: not used");
}

BOOST_AUTO_TEST_CASE(timer_print_named) {
  TestTimer timer{"foo"};
  Delay(timer, 100);
  BOOST_CHECK_EQUAL(Print(timer),
                    "foo                      : avg =  100 ms, "
                    "total =  100 ms, count =         1");
}

BOOST_AUTO_TEST_CASE(scoped_micro_second_accumulator) {
  struct MockAccumulator {
    void operator+=(double d) { accumulated.push_back(d); }
    std::vector<double> accumulated;
  } accumulator;

  const MockClock::time_point now = MockClock::time_point::clock::now();
  mock_time_points.push_back(now);
  mock_time_points.push_back(now + std::chrono::milliseconds(42));
  {
    const dp3::common::ScopedMicroSecondAccumulator<MockAccumulator, MockClock>
        scoped_time(accumulator);
  }
  BOOST_CHECK(mock_time_points.empty());

  BOOST_REQUIRE_EQUAL(accumulator.accumulated.size(), std::size_t(1));
  BOOST_CHECK_CLOSE(accumulator.accumulated.front(), 42000.0, kEpsilon);
}

BOOST_AUTO_TEST_SUITE_END()
