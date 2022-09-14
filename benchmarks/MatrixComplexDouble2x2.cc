// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/MatrixComplexDouble2x2.h"

#include <aocommon/matrix2x2.h>

#include <nanobench.h>
#include <fstream>
#include <iostream>

#if defined(__AVX2__)

template <class Matrix>
static void Conjugate(ankerl::nanobench::Bench& bench, const char* name) {
  Matrix a{{1, 2}, {10, 11}, {100, 101}, {1000, 1001}};

  Matrix r = a;
  bench.run(name, [&] {
    // assigment needed to avoid optimization.
    ankerl::nanobench::doNotOptimizeAway(r = a.Conjugate());
  });
}

template <class Matrix>
static void Transpose(ankerl::nanobench::Bench& bench, const char* name) {
  Matrix a{{1, 2}, {10, 11}, {100, 101}, {1000, 1001}};

  Matrix r = a;
  bench.run(name, [&] {
    // assigment needed to avoid optimization.
    ankerl::nanobench::doNotOptimizeAway(r = a.Transpose());
  });
}

template <class Matrix>
static void HermitianTranspose(ankerl::nanobench::Bench& bench,
                               const char* name) {
  Matrix a{{1, 2}, {10, 11}, {100, 101}, {1000, 1001}};

  Matrix r = a;
  bench.run(name, [&] {
    // assigment needed to avoid optimization.
    ankerl::nanobench::doNotOptimizeAway(r = HermTranspose(a));
  });
}

template <class Matrix>
static void Norm(ankerl::nanobench::Bench& bench, const char* name) {
  Matrix a{{1, 2}, {10, 11}, {100, 101}, {1000, 1001}};

  double r = 1.0;
  bench.run(name, [&] {
    // assigment needed to avoid optimization.
    ankerl::nanobench::doNotOptimizeAway(r = Norm(a));
  });
}

template <class Matrix>
static void Multiply(ankerl::nanobench::Bench& bench, const char* name) {
  Matrix a{{1, 2}, {10, 11}, {100, 101}, {1000, 1001}};
  Matrix b{{4, 8}, {40, 44}, {400, 404}, {4000, 4004}};

  Matrix r = a;
  bench.run(name, [&] {
    // assigment needed to avoid optimization.
    ankerl::nanobench::doNotOptimizeAway(r = a * b);
  });
}

static void Generate(const std::string& extension, char const* output_template,
                     const ankerl::nanobench::Bench& bench) {
  std::ofstream output("MatrixComplextFloat2x2.render." + extension);
  ankerl::nanobench::render(output_template, bench, output);
}

int main() {
  ankerl::nanobench::Bench bench;
  bench.title("Benchmarking MC2x2 normal versus AVX");
  bench.minEpochIterations(500'000);

  Conjugate<aocommon::MC2x2>(bench, "    Conjugate");
  Conjugate<aocommon::Avx256::MatrixComplexDouble2x2>(bench, "AVX Conjugate");
  Conjugate<aocommon::Avx256::MatrixComplexDouble2x2>(bench, "FMV Conjugate");
  Transpose<aocommon::MC2x2>(bench, "    Transpose");
  Transpose<aocommon::Avx256::MatrixComplexDouble2x2>(bench, "AVX Transpose");
  Transpose<aocommon::Avx256::MatrixComplexDouble2x2>(bench, "FMV Transpose");
  HermitianTranspose<aocommon::MC2x2>(bench, "    HermitianTranspose");
  HermitianTranspose<aocommon::Avx256::MatrixComplexDouble2x2>(
      bench, "AVX HermitianTranspose");
  HermitianTranspose<aocommon::MatrixComplexDouble2x2>(
      bench, "FMV HermitianTranspose");
  Norm<aocommon::MC2x2>(bench, "    Norm");
  Norm<aocommon::Avx256::MatrixComplexDouble2x2>(bench, "AVX Norm");
  Norm<aocommon::Avx256::MatrixComplexDouble2x2>(bench, "FMV Norm");
  Multiply<aocommon::MC2x2>(bench, "    Multiply");
  Multiply<aocommon::Avx256::MatrixComplexDouble2x2>(bench, "AVX Multiply");
  Multiply<aocommon::MatrixComplexDouble2x2>(bench, "FMV Multiply");

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
