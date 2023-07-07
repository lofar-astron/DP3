// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_PYTHONDP3_UTILS_H_
#define DP3_PYTHONDP3_UTILS_H_

#include <ostream>
#include <string>

namespace dp3::pythondp3 {

// The purpose of the ostream_wrapper is to make it possible pass a
// C++ std::ostream to python to be used there as text stream.
// A python text stream should implement a write() method that accepts a str.
// The python documentation on the io module explicitly states that passing
// anything else but a str to the write method of a text stream is an error,
// see https://docs.python.org/3/library/io.html#module-io.
// Therefore the wrapper only accepts strings, even though a C++ output stream
// accepts many different types.

class ostream_wrapper {
 public:
  ostream_wrapper(std::ostream& output_stream)
      : output_stream_(output_stream) {}
  void write(std::string& text) { output_stream_ << text; }

 private:
  std::ostream& output_stream_;
};

}  // namespace dp3::pythondp3

#endif  // DP3_PYTHONDP3_UTILS_H_