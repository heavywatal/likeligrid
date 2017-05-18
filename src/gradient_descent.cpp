// -*- mode: c++; coding: utf-8 -*-
/*! @file gradient_descent.cpp
    @brief Implementation of GradientDescent class
*/
#include "gradient_descent.hpp"
#include "util.hpp"

#include <wtl/debug.hpp>
#include <wtl/exception.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>

namespace likeligrid {

GradientDescent::GradientDescent(const std::string& infile,
    const size_t max_sites,
    const unsigned int concurrency)
    : GradientDescent(wtl::izfstream(infile), max_sites, concurrency) {HERE;}

GradientDescent::GradientDescent(
    std::istream& ist,
    const size_t max_sites,
    const unsigned int concurrency)
    : model_(ist, max_sites),
      concurrency_(concurrency) {HERE;

}

void GradientDescent::run() {HERE;

}


void GradientDescent::unit_test() {HERE;
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    GradientDescent searcher(sst, 4);
    searcher.run();
}

} // namespace likeligrid
