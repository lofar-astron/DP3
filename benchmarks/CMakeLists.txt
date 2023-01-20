# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# https://nanobench.ankerl.com/index.html
FetchContent_Declare(
  nanobench
  GIT_REPOSITORY https://github.com/martinus/nanobench.git
  GIT_TAG v4.3.7
  GIT_SHALLOW TRUE)

FetchContent_MakeAvailable(nanobench)

function(add_benchmark name)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} PRIVATE nanobench)

  add_test(NAME benchmark/${name} COMMAND ${name})
  set_tests_properties(benchmark/${name} PROPERTIES LABELS "benchmark")
endfunction()

add_benchmark(MatrixComplexDouble2x2)
add_benchmark(MatrixComplexFloat2x2)
add_benchmark(VectorComplexDouble2)
add_benchmark(VectorComplexFloat2)
