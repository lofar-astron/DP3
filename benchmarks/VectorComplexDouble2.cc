// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/avx256/VectorComplexDouble2.h"

#include <nanobench.h>
#include <fstream>
#include <iostream>

#if defined(__AVX2__)
static void Multiply(ankerl::nanobench::Bench& bench, const char* name) {
  aocommon::Avx256::VectorComplexDouble2 a{{1, 2}, {10, 11}};
  std::complex<double> b{10, 8};

  bench.run(name, [&] {
    a = a * b;
    ankerl::nanobench::doNotOptimizeAway(a);
  });
}

static void Generate(const std::string& extension, char const* output_template,
                     const ankerl::nanobench::Bench& bench) {
  std::ofstream output("MatrixComplextFloat2x2.render." + extension);
  ankerl::nanobench::render(output_template, bench, output);
}

int main() {
  ankerl::nanobench::Bench bench;
  bench.title("Benchmarking Vector Complex 2 AVX");
  bench.minEpochIterations(1'000'000);

  Multiply(bench, "Multiply Vector and Value");

  Generate("html", ankerl::nanobench::templates::htmlBoxplot(), bench);
  Generate("json", ankerl::nanobench::templates::json(), bench);
}

#else
int main() {
  std::cerr << R"(Benchmarking "MC2x2 normal versus AVX" requires AVX support.
Try rebuilding the benchmark after setting -DPORTABLE=OFF in CMake.
)";

  return EXIT_FAILURE;
}
#endif
