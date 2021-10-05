// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Timer.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <atomic>
#include <numeric>
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

template <class T>
static void Delay(T& accumulator, const int delay) {
  using ScopedTime = dp3::common::ScopedMicroSecondAccumulator<T>;
  const ScopedTime scoped_time(accumulator);
  std::this_thread::sleep_for(std::chrono::milliseconds(delay));
}

BOOST_AUTO_TEST_CASE(scoped_micro_second_accumumator_start_zero) {
  int accumulator = 0;
  Delay(accumulator, 1000);

  // Note the timer is known to be not too accurate.
  BOOST_CHECK_CLOSE(accumulator, 1e6, 5);
}

BOOST_AUTO_TEST_CASE(scoped_micro_second_accumumator_start_non_zero) {
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

BOOST_AUTO_TEST_CASE(scoped_micro_second_accumumator_atomic_integral) {
  constexpr std::array<int, 12> delays{5, 5, 3, 4, 2, 1, 1, 3, 4, 2, 2, 1};
  const double expected =
      1e3 * std::accumulate(delays.begin(), delays.end(), 0);

  const double result = RunAllDelaysInThreads(delays.begin(), delays.end());

  // Note the timer is known to be not too accurate.
  BOOST_CHECK_CLOSE(result, expected, 5);
}

BOOST_AUTO_TEST_SUITE_END()
